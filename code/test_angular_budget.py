"""Independent checks of algebra, bounds, and floating-point implementation."""
import math
from fractions import Fraction as Fr
import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from angular_budget import (symbolic_checks, sigma, log_ramp, log_G,
                            gaussian_envelope, positivity_budget,
                            inherited_response, original_recipe_log_B)


def test_symbolic_identities():
    assert all(symbolic_checks().values())


@pytest.mark.parametrize('z',[.001,.0025,.005,.0075,.01,.015,.019,.02,.2,1.,10.])
def test_log_ramp_independent_quadrature(z):
    # Integral of the explicit smooth step, not the log-quadrature implementation.
    pieces=[0.,min(.01,z)]
    if z>.01: pieces.append(min(.02,z))
    if z>.02: pieces.append(z)
    direct=sum(quad(lambda t:sigma(50*t),a,b,epsabs=1e-150,epsrel=2e-11,limit=200)[0]
               for a,b in zip(pieces[:-1],pieces[1:]))
    calculated=math.exp(log_ramp(z))
    assert abs(calculated-direct)<= max(1e-150,abs(direct)*1e-9)


@pytest.mark.parametrize('lam',[1e-5,1e-4,.001,.005,.008,.02,.08])
def test_global_envelope_samples(lam):
    cap=gaussian_envelope(lam)['G_sup_upper']
    for t in np.geomspace(.01,11/lam,501):
        lg=log_G(t,lam)
        assert lg<=math.log(cap)+1e-13


@pytest.mark.parametrize('lam',[1e-5,1e-4,.001,.008,.08])
def test_inherited_response_maximum(lam):
    c=1-lam; h=lam/100000; eta=.4; S0=1.
    tstar=-math.log1p(-lam)/lam
    exact=lam*math.exp(c/lam*math.log1p(-lam))*(1-2*h*eta**2)
    attained=-float(inherited_response(tstar,lam,eta,h,S0))
    assert math.isclose(exact,attained,rel_tol=2e-12)
    for t in np.linspace(0,20,1001):
        assert abs(inherited_response(t,lam,eta,h,S0)) <= exact*(1+1e-12)


@pytest.mark.parametrize('lam',[1e-5,1e-4,.001,.005,.008])
def test_budget_positive(lam):
    b=positivity_budget(lam,Fr(str(lam))/100000)
    assert b['B_sufficient_upper_limit_lower']>0
    assert b['min_equilibrium_Q']>0


def test_general_source_identity_randomized():
    rng=np.random.default_rng(917)
    for _ in range(1000):
        l=rng.uniform(.0001,.08); h=min(l/3,.009)
        eta=rng.uniform(-1,1); E=math.exp(rng.uniform(-3,3))
        n,ne,R=rng.normal(size=3)
        D=.5-h;sg=2*eta/(1+eta**2)
        X=E*n; Xe=E*ne-sg*X;Y=E*R
        full= l*(1-2*D*eta*E*n-(1-eta**2)*E*(ne-sg*n)) \
             -h*(1-2*eta*E*R)+(D*eta+(1-eta**2)*E*R)*sg
        reduced=l-h+D*eta*sg-l*(2*D*eta*X+(1-eta**2)*Xe) \
                +(2*h*eta+(1-eta**2)*sg)*Y
        assert math.isclose(full,reduced,rel_tol=2e-11,abs_tol=2e-11)


def test_direct_product_ode_arbitrary_input():
    lam=.008;alpha=.5+lam;beta=.5-lam;eta=.37
    amp=1+.1*eta;B=1.7;q0=.3
    init=q0*eta*(1+eta**2)
    def R(t): return amp*(.2+math.sin(t)) # deliberately not the main pulse
    def rhs(t,v):
        E=B/(1+eta**2)*math.exp(-alpha*t)
        return [R(t)-beta*v[0], E*R(t)-v[1]]
    E0=B/(1+eta**2)
    sol=solve_ivp(rhs,(0,12),[init,E0*init],rtol=1e-11,atol=1e-13,max_step=.05)
    actual=B/(1+eta**2)*np.exp(-alpha*sol.t)*sol.y[0]
    assert np.max(np.abs(actual-sol.y[1]))<2e-10


def test_flat_edge_two_quadrature_orders():
    for lam in [1e-5,1e-4,.001]:
        for t in np.linspace(2,400,80):
            log128=log_G(t,lam,128)
            # At extremely flat tails neither formula affects numeric quadrature.
            if math.isfinite(log128) and log128>-500:
                assert abs(log128-log_G(t,lam,256))<1e-7


@pytest.mark.parametrize('lam,h',[(-1,.0),(.2,.0),(.001,.002),(.05,.02)])
def test_invalid_parameters_rejected(lam,h):
    with pytest.raises(ValueError): positivity_budget(lam,h)


def test_finite_example_has_margin():
    lam=.001
    B=math.exp(original_recipe_log_B(6.,lam))
    bound=positivity_budget(lam,1e-8)
    assert 0<B<bound['B_sufficient_upper_limit_lower']
