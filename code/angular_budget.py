"""v0.14: cancellation-aware angular-lag budget.

Independent variable y is a LOG-RADIUS coordinate, not elapsed physical time.

A local scalar subsystem of the pinned outgoing-profile construction.
NOT a Navier--Stokes solver, full cone certificate, or singularity proof.

The envelope coefficients use exact rational root bracketing and mpmath
interval elementary functions. The ODE checks are ordinary floating point.
"""
from __future__ import annotations
import argparse, csv, json, math, platform
from fractions import Fraction as Fr
from pathlib import Path
from typing import Any
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar
from scipy.special import expit, log_expit, logsumexp
import sympy as sp
import mpmath as mp

mp.iv.dps = 60
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(128)


def _fr(x: Any) -> Fr:
    return x if isinstance(x, Fr) else Fr(str(x))


def _iv(x: Any):
    f = _fr(x)
    return mp.iv.mpf(f.numerator) / mp.iv.mpf(f.denominator)


def _up(x) -> float:
    return math.nextafter(math.nextafter(float(x.b), math.inf), math.inf)


def _down(x) -> float:
    return math.nextafter(math.nextafter(float(x.a), -math.inf), -math.inf)


def sigma(x: float) -> float:
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    return float(expit(1/(1-x)**2 - 1/x**2))


def log_ramp(z: float, order: int = 128) -> float:
    """Log of integral_0^z sigma(50 s) ds; stable near the flat edge."""
    if z <= 0: return -math.inf
    if z >= .02: return math.log(z-.01)
    x = 50*z
    reflect = x > .5
    endpoint = 1-x if reflect else x
    nodes, weights = (_NODES, _WEIGHTS) if order == 128 else np.polynomial.legendre.leggauss(order)
    u = endpoint*(nodes+1)/2
    logs = log_expit(1/(1-u)**2 - 1/u**2)
    value = math.log(endpoint/2) + float(logsumexp(np.log(weights)+logs))
    if reflect:
        value = float(np.logaddexp(math.log(x-.5),value))
    return value-math.log(50)


def log_main(z: float, order: int = 128) -> float:
    if z <= 0 or z >= 11: return -math.inf
    lr = log_ramp(z, order)
    x = z-10
    if x <= 0: return lr
    return lr+float(log_expit(1/x**2 - 1/(1-x)**2))


def log_G(y: float, lam: float, order: int = 128) -> float:
    return -(.5+lam)*y+log_main(lam*y,order)


def gaussian_envelope(lam: Any) -> dict[str, Any]:
    """Rigorous upper envelope of exp(-alpha*y) mainPulse(lambda*y).

    For z<=.01, ramp(z)<=1250*z^3*exp(4-1/(2500*z^2)).
    For z>=.01, ramp(z)<=z-7/800 (the .01 endpoint bound is 1/800).
    The first envelope's maximum is bracketed with exact rational arithmetic.
    Values of log/exp are enclosed by mpmath.iv, at 60 decimal digits.
    """
    l = _fr(lam)
    if not 0 < l < Fr(1,10): raise ValueError('0 < lambda < 0.1 required')
    alpha=Fr(1,2)+l; y0=Fr(1,100)/l; k=1/(2500*l*l)
    poly=lambda y:alpha*y**3-3*y*y-2*k
    if poly(y0) <= 0:
        lo=hi=y0
    else:
        lo=Fr(0); hi=y0
        for _ in range(130):
            mid=(lo+hi)/2
            if poly(mid)>0: hi=mid
            else: lo=mid
    # Each varying factor is bounded in its own monotonic direction.
    capA=1250*_iv(l*hi)**3*mp.iv.exp(4-_iv(alpha*lo)-_iv(k/(hi*hi)))
    yB=max(y0, Fr(7,800)/l+1/alpha)
    capB=_iv(l*yB-Fr(7,800))*mp.iv.exp(-_iv(alpha*yB))
    generic=_iv(l/alpha)/mp.iv.exp(1)
    # Upper endpoints are rounded outward again on conversion to float.
    cap=min(max(_up(capA),_up(capB)),_up(generic))
    return {'lambda':float(l),'G_sup_upper':cap,'log10_G_sup_upper':math.log10(cap),
            'first_region_root_lo':str(lo),'first_region_root_hi':str(hi),
            'interval_decimal_digits':mp.iv.dps,
            'scope':'analytic envelope only; does not certify an ODE or PDE'}


def positivity_budget(lam: Any, h: Any, amp_max: Any = Fr(6,5), amp_deriv_max: Any=1) -> dict[str,float]:
    """Whole pre-repair angular positivity under listed scalar assumptions.

    0<lambda<.1, 0<=h<=.01, h<lambda; original entrance/hold formulas;
    smooth amplitude with |a|<=amp_max, |a'|<=amp_deriv_max.
    No claim that these parameters satisfy the full construction.
    """
    l=_fr(lam); hh=_fr(h); aa=_fr(amp_max); ap=_fr(amp_deriv_max)
    if not (0<l<Fr(1,10) and 0<=hh<=Fr(1,100) and hh<l and aa>=0 and ap>=0):
        raise ValueError('invalid scalar assumptions')
    c=1-l
    envelope=gaussian_envelope(l)
    # Upper bound on entrance angular disequilibrium after the original hold.
    e0=10*mp.iv.exp(_iv(60*c)*mp.iv.log(_iv(l)))
    # Max inherited-mass response, with S0 <= 4 e^-12 lambda^60.
    dip=4*mp.iv.exp(-12)*mp.iv.exp(_iv(61)*mp.iv.log(_iv(l))+
                                _iv(c/l)*mp.iv.log(_iv(c)))
    coeff=l*(3*aa+ap)+3*aa
    available=_iv(l-hh)-_iv(c)*(e0+dip)
    bmax=available/(_iv(coeff)*_iv(envelope['G_sup_upper']))
    return {'lambda':float(l),'h':float(hh),
            'G_sup_upper':envelope['G_sup_upper'],
            'incoming_Q_error_upper':_up(e0),'inherited_Q_dip_upper':_up(dip),
            'source_amplitude_coefficient':float(coeff),
            'B_sufficient_upper_limit_lower':_down(bmax),
            'min_equilibrium_Q':float((l-hh)/c)}


def remaining_repair_budget(lam: Any, h: Any, B: Any) -> dict[str,float]:
    """For Y_repair and its eta derivative, require
    (3+lambda)*||Y_repair|| + lambda*||d_eta Y_repair|| < available.
    The returned number is an acceptance threshold, not a measured repair.
    """
    l=_fr(lam);hh=_fr(h);bb=_fr(B)
    if bb<0: raise ValueError('B must be nonnegative')
    budget=positivity_budget(l,hh)
    c=1-l
    e0=10*mp.iv.exp(_iv(60*c)*mp.iv.log(_iv(l)))
    dip=4*mp.iv.exp(-12)*mp.iv.exp(_iv(61)*mp.iv.log(_iv(l))+_iv(c/l)*mp.iv.log(_iv(c)))
    main=_iv(l*Fr(23,5)+Fr(18,5))*_iv(bb)*_iv(budget['G_sup_upper'])
    available=_iv(l-hh)-_iv(c)*(e0+dip)-main
    return {'lambda':float(l),'h':float(hh),'B':float(bb),
      'remaining_source_budget_lower':_down(available),
      'uniform_Q_lower_bound':_down(available/_iv(c)),
      'fraction_of_min_equilibrium_lower':_down(available/_iv(l-hh)),
      'equal_norm_repair_limit_lower':_down(available/_iv(3+2*l))}


def inherited_response(y: np.ndarray|float, lam: float, eta: float, h: float, S0: float) -> np.ndarray:
    y=np.asarray(y,dtype=float); c=1-lam
    return -S0*(1-2*h*eta*eta)*np.exp(-c*y)*(-np.expm1(-lam*y))


def numerical_kernels(lam: float, points: int=2001) -> dict[str, Any]:
    """Floating-point diagnostic kernels near the effective weighted pulse.

    The time window is finite and recorded. The global sufficient bound is
    analytic and is not inferred from this grid.
    """
    a=.5+lam; k=1/(2500*lam*lam)
    guess=(2*k/a)**(1/3)
    opt=minimize_scalar(lambda t:-log_G(t,lam),bounds=(1e-8,max(8*guess,.04/lam)),
                        method='bounded',options={'xatol':1e-10})
    peak_y=float(opt.x); peak_log=-float(opt.fun)
    t_end=min(11/lam,2*peak_y+100/a)
    def forcing(t):
        lg=log_G(t,lam)
        return 0.0 if not math.isfinite(lg) or lg-peak_log < -100 else math.exp(lg-peak_log)
    c=1-lam
    def rhs(t,x):
        g=forcing(t)
        return [g-x[0],g-c*x[1],x[0]-c*x[2]]
    sol=solve_ivp(rhs,(0,t_end),[0.,0.,0.],method='DOP853',rtol=2e-11,atol=3e-13,
                  dense_output=True,max_step=min(1.,max(peak_y/20,.01)))
    if not sol.success: raise ArithmeticError(sol.message)
    times=np.unique(np.r_[np.linspace(0,t_end,points),peak_y])
    vals=sol.sol(times)
    # Independent scalar quadrature at four times. G is rescaled by its peak.
    maxerr=0.
    for t in [peak_y*.65,peak_y,peak_y+1,peak_y+8]:
        fun=lambda s:forcing(s)
        H=quad(lambda s:math.exp(-(t-s))*fun(s),0,t,epsabs=2e-11,limit=150)[0]
        IG=quad(lambda s:math.exp(-c*(t-s))*fun(s),0,t,epsabs=2e-11,limit=150)[0]
        # Convolution of exp(-t) and exp(-c*t), evaluated without cancellation.
        IH=quad(lambda s:math.exp(-c*(t-s))*(-math.expm1(-lam*(t-s)))/lam*fun(s),
                0,t,epsabs=2e-11,limit=150)[0]
        maxerr=max(maxerr,float(np.max(np.abs(np.array([H,IG,IH])-sol.sol(t)))))
    eta=np.unique(np.r_[np.linspace(0,1,401),np.geomspace(1e-7,1,401)])
    h=lam/100000
    ss=eta*eta; phi=1/(1+ss); sg=2*eta/(1+ss); D=.5-h
    AA=2*D*eta-(1-ss)*sg
    CC=2*h*eta+(1-ss)*sg
    eq=(lam-h+D*eta*sg)/(1-lam)
    M=phi[:,None]*(CC[:,None]*vals[1]-lam*AA[:,None]*vals[2])
    relative=np.abs(M)/eq[:,None]
    index=np.unravel_index(np.argmax(relative),relative.shape)
    gain=float(relative[index])*math.exp(peak_log)
    critical=1/gain
    # The constant-amplitude model has symmetry in eta. The adverse sign is used.
    return {'lambda':lam,'times':times,'values':vals,
      'weighted_pulse_peak_y':peak_y,'weighted_pulse_peak':math.exp(peak_log),
      'sampled_logradial_coordinate_end':t_end,'independent_quadrature_max_scaled_error':maxerr,
      'constant_amp_sampled_B_critical':critical,'adverse_eta_magnitude':float(eta[index[0]]),
      'adverse_y':float(times[index[1]]),'solver_evaluations':sol.nfev,
      'response_identity_max_error':float(np.max(np.abs(vals[1]-vals[0]-lam*vals[2]))),
      'kernel_peak_coordinate':float(times[np.argmax(vals[1])])}


def symbolic_checks() -> dict[str,bool]:
    l,h,eta,X,Xe,Y,n,ne,E=sp.symbols('lambda h eta X Xeta Y n neta E',real=True)
    d=sp.Rational(1,2)-h; g=2*eta/(1+eta**2)
    A=2*d*eta-(1-eta**2)*g; C=2*h*eta+(1-eta**2)*g
    original=-l*(A*E*n+(1-eta**2)*E*ne)+C*Y
    replaced=original.subs({n:X/E,ne:(Xe+g*X)/E})
    compact=-l*(2*d*eta*X+(1-eta**2)*Xe)+C*Y
    S,t=sp.symbols('S t',real=True)
    inherited=compact.subs({X:S*eta*sp.exp(-t),Xe:S*sp.exp(-t),Y:0})
    response=-S*(1-2*h*eta**2)*(sp.exp(-(1-l)*t)-sp.exp(-t))
    return {
      'product_ODE_rate_exactly_one':sp.simplify((sp.Rational(1,2)+l)+(sp.Rational(1,2)-l)-1)==0,
      'full_angular_source_product_identity':sp.simplify(replaced-compact)==0,
      'inherited_source_exact':sp.simplify(inherited+l*S*(1-2*h*eta**2)*sp.exp(-t))==0,
      'inherited_angular_response_ODE':sp.simplify(sp.diff(response,t)+(1-l)*response-inherited)==0,
      'inherited_response_initial_zero':sp.simplify(response.subs(t,0))==0}


def original_recipe_log_B(m: float,lam: float) -> float:
    K=math.exp(m)
    lp=float(np.logaddexp(K+11,0.))
    return lp+.3-lam/2-(K+12)/2+(30+60*lam)*math.log(lam)


def write_csv(p:Path,rows:list[dict]) -> None:
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    lams=[1e-5,1e-4,1e-3,.005,.008]
    envelopes=[]; budgets=[]; numerical=[]; paramrows=[]
    for lam in lams:
        env=gaussian_envelope(lam); envelopes.append(env)
        budget=positivity_budget(lam,Fr(str(lam))/100000);budgets.append(budget)
        n=numerical_kernels(lam)
        np.savez_compressed(args.output/f'kernels_lambda_{lam:g}.npz',
                            logradial_coordinate=n.pop('times'),values=n.pop('values'))
        numerical.append(n)
        assert n['weighted_pulse_peak'] <= env['G_sup_upper']
        assert n['independent_quadrature_max_scaled_error'] < 1e-8
        for m in [.1,1.,3.,5.5,6.,6.2]:
            logB=original_recipe_log_B(m,lam)
            budgetratio=math.exp(min(700,logB-math.log(budget['B_sufficient_upper_limit_lower'])))
            paramrows.append({'lambda':lam,'m':m,'log_B':logB,'B':math.exp(logB),
                'fraction_of_sufficient_B_limit':budgetratio,
                'angular_only_sufficient_condition':budgetratio<1,
                'full_witness_claim':False})
    checks=symbolic_checks();assert all(checks.values())
    repair=[remaining_repair_budget(l,Fr(str(l))/100000,b) for l in lams for b in [0,.1,1.]]
    write_csv(args.output/'repair_acceptance_budgets.csv',repair)
    write_csv(args.output/'angular_budgets.csv',budgets)
    write_csv(args.output/'numerical_kernel_checks.csv',numerical)
    write_csv(args.output/'local_recipe_screen.csv',paramrows)
    (args.output/'envelope_certificates.json').write_text(json.dumps(envelopes,indent=2),encoding='utf-8')
    summary={'symbolic_checks':checks,'lambda_cases':len(lams),'local_parameter_screens':len(paramrows),
      'max_independent_quadrature_scaled_error':max(x['independent_quadrature_max_scaled_error'] for x in numerical),
      'max_kernel_identity_error':max(x['response_identity_max_error'] for x in numerical),
      'conditions_satisfied_in_local_screen':sum(r['angular_only_sufficient_condition'] for r in paramrows),
      'status':{'scalar_envelope_interval_evaluated':True,'interval_digits':mp.iv.dps,
                'ODE_interval_certified':False,'PDE_solved':False,'full_cone_checked':False,
                'repair_region_computed':False,'lean_compiled':False,'novelty_established':False,'physical_time_interpretation':False,
                'floating_ODE_relative_forcing_clip':math.exp(-100)}}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (args.output/'environment.json').write_text(json.dumps({'python':platform.python_version(),
      'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,'sympy':sp.__version__},indent=2))
    print(json.dumps(summary,indent=2))
    for b,n in zip(budgets,numerical):
        print('lambda',b['lambda'],'analytic B_limit',b['B_sufficient_upper_limit_lower'],
              'sampled constant-a critical',n['constant_amp_sampled_B_critical'])

if __name__=='__main__': main()
