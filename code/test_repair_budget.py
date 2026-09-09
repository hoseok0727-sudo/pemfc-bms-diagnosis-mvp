from fractions import Fraction as F
import math
import numpy as np
import pytest
from scipy.integrate import quad
from repair_budget import (geometric_certificates, full_certificate, numerical_repair,
    moment, V_matrix, stable_V_solve, bump_shape, cutoff, smooth_transition,
    repair_bounds, log10_interval, S)
from exact_example import certificate, exp_lower, exp_upper


def test_exact_rational_certificate():
    assert certificate()['all_exact_checks_pass']


def test_exact_exp_enclosures():
    for x in [F(0),F(1,10),F(2),F(4307,800)]:
        assert exp_lower(x)<exp_upper(x) if x else exp_lower(x)==exp_upper(x)==1


def test_shape_independent_constants():
    g=geometric_certificates()
    assert g['moment_floor_pass'] and g['inverse_factor_pass'] and g['field_weight_pass']


@pytest.mark.parametrize('power',[1,2])
def test_cutoff_required_properties(power):
    grid=np.linspace(-2,2,2001)
    vals=np.array([cutoff(x,power) for x in grid])
    assert np.all((vals>=0)&(vals<=1))
    assert all(cutoff(x,power)==1 for x in [-.5,-.1,0,.5])
    assert all(cutoff(x,power)==0 for x in [-2,-1,1,2])
    assert max(abs(cutoff(x,power)-cutoff(-x,power)) for x in grid)==0


@pytest.mark.parametrize('l',[.00001,.0001,.001,.005,.008,F(1,120)])
def test_full_certificate_example_domains(l):
    l=F(str(l)) if not isinstance(l,F) else l
    r=full_certificate(l,l/100000,F(1,10))
    assert r['positive_certified']
    assert r['repair_source_log10_upper']<-370
    assert 0<r['full_Q_lower']<=r['min_equilibrium_Q']


@pytest.mark.parametrize('l',[.00001,.001,.008])
def test_inverse_determinant_and_bound(l):
    V=V_matrix(l)
    exact=-math.exp(-2+7*l)*(-math.expm1(-2*l))
    assert math.isclose(np.linalg.det(V),exact,rel_tol=1e-10)
    assert np.linalg.norm(np.linalg.inv(V),np.inf)<8/l
    for rhs in [np.array([1.,0.]),np.array([0.,1.]),np.array([2.,-3.])]:
        z=stable_V_solve(l,rhs)
        assert np.linalg.norm(V@z-rhs,np.inf)<1e-9


@pytest.mark.parametrize('power',[1,2])
def test_moment_matrix_factorization_direct_integral(power):
    l=.007; L=8.; betas=[.5-l,.5-2*l]
    for i,beta in enumerate(betas):
        for a in [3.,1.]:
            direct=quad(lambda y: math.exp(beta*y)*bump_shape(y-L+a,power),
                        L-a-float(S),L-a+float(S),epsabs=1e-10,limit=200)[0]
            factored=math.exp(beta*(L-a))*moment(beta,power)
            assert math.isclose(direct,factored,rel_tol=3e-10)


@pytest.mark.parametrize('l,power',[(.00001,1),(.001,1),(.008,1),(.001,2)])
def test_normalized_moment_solve_and_majorant(l,power):
    r=numerical_repair(l,m=1.,B=1.,power=power)
    assert r['scaled_moment_residual']<1e-9
    assert r['V_solve_vs_80digit_relative_error']<2e-13
    assert r['small_moment_min']>.05
    T0,Tm,_,_,_=repair_bounds(l,1)
    # Each representative's weighted component is below the universal bound.
    import mpmath as mp
    pref=log10_interval(mp.iv.mpf(800)/mp.iv.mpf(str(l))*T0)['upper']
    main=log10_interval(mp.iv.mpf(800)/mp.iv.mpf(str(l))*Tm)['upper']
    assert r['prefix_weighted_log10_sampled']<pref
    assert r['main_weighted_log10_sampled']<main


def test_B_monotonicity_and_threshold():
    l=F(1,1000);h=l/100000
    r1=full_certificate(l,h,1);r2=full_certificate(l,h,10);r3=full_certificate(l,h,40)
    assert r1['full_Q_lower']>r2['full_Q_lower']>r3['full_Q_lower']
    assert r1['positive_certified'] and not r3['positive_certified']
    assert 30<r1['B_strict_sufficient_limit_lower']<31


@pytest.mark.parametrize('l,h,B',[(-1,0,1),(0,0,1),(.1,0,1),(.001,.001,1),(.001,0,-1)])
def test_invalid_parameters(l,h,B):
    with pytest.raises(ValueError):full_certificate(l,h,B)


def test_angular_exact_repair_kernel_identity():
    # Independently checks the exact source response for one representative
    # repair pulse (normalized unit coefficient) against direct coupled ODE.
    from scipy.integrate import solve_ivp
    l=.001;h=1e-8;eta=.3;c=1-l;D=.5-h
    ph=1/(1+eta**2);g=2*eta/(1+eta**2)
    Y=lambda t:ph*math.exp(-(.5+l)*t)*bump_shape(t-1.)
    coeff=2*h*eta+(1-eta**2)*g
    def rhs(t,z):
        y=Y(t);ye=-g*y
        return [y-z[0],ye-z[1],-l*(2*D*eta*z[0]+(1-eta**2)*z[1])+coeff*y-c*z[2]]
    sol=solve_ivp(rhs,(0,5),[0.,0.,0.],method='DOP853',rtol=1e-11,atol=1e-13,max_step=.01,dense_output=True)
    assert sol.success
    for t in [1.,1.2,2.,4.]:
        H=quad(lambda s:math.exp(-c*(t-s))*(-math.expm1(-l*(t-s)))/l*Y(s),0,t,
               points=[x for x in [.8,1.,1.2] if x<t],epsabs=1e-12)[0]
        IG=quad(lambda s:math.exp(-c*(t-s))*Y(s),0,t,
                points=[x for x in [.8,1.,1.2] if x<t],epsabs=1e-12)[0]
        exact=-l*(2*D*eta-(1-eta**2)*g)*H+coeff*IG
        assert abs(exact-sol.sol(t)[2])<1e-11
