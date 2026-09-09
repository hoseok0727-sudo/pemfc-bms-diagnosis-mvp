"""Cancellation-aware audit of a local outgoing-pulse construction.

This is NOT a Navier--Stokes solver, a blow-up verification, or a Lean build.
The source equations are listed in RESEARCH_NOTE.md at a pinned commit.

Run: python pulse_audit.py --output results
Dependencies: numpy, scipy, sympy, mpmath, matplotlib
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
import sympy as sp


def sigma(x: float) -> float:
    """Smooth symmetric step, edge(x)=exp(-1/x**2) for x>0."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    return float(expit(1.0/(1.0-x)**2 - 1.0/x**2))


def drop(m: float, y: float) -> float:
    """4*(1-sigma(log(y)/m)), evaluated without subtractive cancellation."""
    if m <= 0:
        raise ValueError('m must be positive')
    if y <= 1.0:
        return 4.0
    z = math.log(y)/m
    if z >= 1.0:
        return 0.0
    return 4.0 * sigma(1.0-z)


def integral_sigma(x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return x-0.5
    return quad(sigma,0.0,x,epsabs=2e-13,epsrel=2e-13)[0]


def log_radial_over_P(y: float, D: float, lam: float) -> float:
    """Exact primitive representation; no huge amplitudes are exponentiated."""
    if y < 0:
        return y/10.0
    return 0.6*(y-integral_sigma(y))-lam*integral_sigma(y-D-1)-y/2


def slope(y: float,D: float,lam: float) -> float:
    return 0.6*sigma(1.0-y)-lam*sigma(y-D-1.0)


def check_params(m: float,lam: float,logP: float) -> None:
    if not all(math.isfinite(t) for t in (m,lam,logP)):
        raise ValueError('parameters must be finite')
    if m <= 0 or not 0 < lam < 0.1:
        raise ValueError('require m>0 and 0<lambda<0.1')


def numerical_case(m: float, lam: float, logP: float) -> dict[str, float]:
    check_params(m,lam,logP)
    K=math.exp(m);D=K+10.0;H=K+12.0;w=60.0*math.log(1/lam)
    # Direct piecewise numerical integration of the original slope law.
    edges=[0.0,1.0,D+1.0,H,H+w]
    numeric_logA=logP+sum(quad(lambda y:slope(y,D,lam)-0.5,a,b,
                             epsabs=5e-12,epsrel=5e-13)[0]
                        for a,b in zip(edges[:-1],edges[1:]))
    formula_logA=logP+0.3-lam/2-H/2+(30+60*lam)*math.log(lam)
    # Prefix moments normalized by exp(K), without catastrophic overflow.
    # All integrands are nonnegative. Splitting at 1 handles the flat prefix.
    def integrate(fn):
        return sum(quad(fn,a,b,epsabs=1e-13,epsrel=1e-11,limit=200)[0]
                   for a,b in [(0.,1.),(1.,K)])
    Mscaled=4*math.exp(-K)+integrate(lambda y:math.exp(y-K)*drop(m,y))
    Jscaled=2.5*math.exp(-K)+integrate(
        lambda y: math.exp(1.5*y+log_radial_over_P(y,D,lam)-K)*drop(m,y))
    if Mscaled <= 0 or Jscaled <= 0:
        raise ArithmeticError('normalized quadrature returned a nonpositive value')
    # M=exp(K)*Mscaled; J=sqrt(2)*P*exp(K)*Jscaled.
    logM=K+math.log(Mscaled)
    logJ=0.5*math.log(2)+logP+K+math.log(Jscaled)
    yb=H+w
    logq0=logM-yb-numeric_logA
    logq1=logJ-0.5*math.log(2)-1.5*yb-2*numeric_logA
    # Two independent ways to calculate A_pulse*q_i.
    direct_logprod0=numeric_logA+logq0
    direct_logprod1=numeric_logA+logq1
    normalized_logprod0=math.log(Mscaled)-12-w
    normalized_logprod1=math.log(Jscaled)-0.3+lam/2-12-(1-lam)*w
    # Analytic upper bounds, valid without P>=amplitudeThreshold.
    cap0=math.log(4)-12+60*math.log(lam)
    cap1=math.log(4)-12+lam/2+60*(1-lam)*math.log(lam)
    hold_exact=logP+0.3-lam/2-H/2
    old_envelope=logP+H
    return {
      'm':m,'lambda':lam,'logP':logP,
      'pulse_log_amplitude_abs_error':abs(numeric_logA-formula_logA),
      'history0_log_identity_abs_error':abs(direct_logprod0-normalized_logprod0),
      'history1_log_identity_abs_error':abs(direct_logprod1-normalized_logprod1),
      'history0_ratio_to_analytic_cap':math.exp(normalized_logprod0-cap0),
      'history1_ratio_to_analytic_cap':math.exp(normalized_logprod1-cap1),
      'history0_log_product':normalized_logprod0,
      'history1_log_product':normalized_logprod1,
      'old_over_sharp_amplitude_coefficient_log10':(old_envelope-hold_exact)/math.log(10),
    }


def symbolic_checks() -> dict[str,bool]:
    K,H,lam,w,lp=sp.symbols('K H lam w logP',real=True)
    la=lp+sp.Rational(3,10)-lam/2-H/2-(sp.Rational(1,2)+lam)*w
    lm,lj=sp.symbols('logM logJ',real=True)
    lq0=lm-(H+w)-la
    lq1=lj-sp.log(2)/2-sp.Rational(3,2)*(H+w)-2*la
    # The second target still contains the exact J; substitute its scaled form.
    jscaled,mscaled=sp.symbols('logJscaled logMscaled',real=True)
    ident0=sp.expand((la+lq0).subs({lm:K+mscaled,H:K+12}))
    ident1=sp.expand((la+lq1).subs({lj:sp.log(2)/2+lp+K+jscaled,H:K+12}))
    j,m,h=sp.symbols('J m h',real=True)
    rate=h*(j/10+sp.Rational(7,10))-(m*(m+2)+sp.Rational(5,2)+2*h+2*h*m)
    threshold=(10/h)*(m*m+2*m+sp.Rational(5,2))+13+20*m
    return {
      'history0_exact_cancellation':sp.simplify(ident0-(mscaled-12-w))==0,
      'history1_exact_cancellation':sp.simplify(ident1-(jscaled-sp.Rational(3,10)+lam/2-12-(1-lam)*w))==0,
      'published_ledger_threshold_identity':sp.simplify(rate-h/10*(j-threshold))==0,
      'old_sharp_log_ratio':sp.simplify((lp+H)-(lp+sp.Rational(3,10)-lam/2-H/2)
                         -(sp.Rational(3,2)*H-sp.Rational(3,10)+lam/2))==0,
    }


def independent_high_precision() -> dict[str,Any]:
    import mpmath as mp
    mp.mp.dps=60
    def sig(x):
        if x<=0: return mp.mpf(0)
        if x>=1: return mp.mpf(1)
        return 1/(1+mp.exp(1/x**2-1/(1-x)**2))
    I=mp.quad(sig,[0,mp.mpf('0.25'),mp.mpf('0.5'),mp.mpf('0.75'),1])
    l=mp.mpf('0.001');m=mp.mpf(2);H=mp.exp(m)+12;logP=mp.mpf(7)
    # Piecewise integration; independent arbitrary-precision step implementation.
    I1=mp.quad(lambda t:mp.mpf('0.6')*(1-sig(t))-mp.mpf('0.5'),[0,mp.mpf('.5'),1])
    I2=-(H-2)/2
    I3=mp.quad(lambda t:-l*sig(t)-mp.mpf('.5'),[0,mp.mpf('.5'),1])
    actual=logP+I1+I2+I3-(mp.mpf('.5')+l)*60*mp.log(1/l)
    expected=logP+mp.mpf('.3')-l/2-H/2+(30+60*l)*mp.log(l)
    return {'precision_decimal_digits':60,
      'symmetric_step_integral_abs_error':str(abs(I-mp.mpf('.5'))),
      'pulse_log_formula_abs_error':str(abs(actual-expected)),
      'log10_sharp_history0_cap_lambda_0_001':str((mp.log(4)-12+60*mp.log(l))/mp.log(10)),
      'log10_sharp_history1_cap_lambda_0_001':str((mp.log(4)-12+l/2+60*(1-l)*mp.log(l))/mp.log(10))}


def exact_integer_ledger() -> list[dict[str,int]]:
    # Exact rational arithmetic fixes off-by-one float rounding in v0.7.
    from fractions import Fraction
    h=Fraction(1,1000)
    out=[]
    for m in [0,1,2,3,5,10]:
        threshold=10/h*(m*m+2*m+Fraction(5,2))+13+20*m
        first=threshold.numerator//threshold.denominator+1
        out.append({'derivative_order':m,'first_J_with_positive_ledger_rate':first})
    return out


def cancellation_scaling_cases() -> list[dict[str,float]]:
    """Exact scale-only test; no ODE or PDE robustness inference."""
    B=np.array([[1.,1.],[1.,-1.]])/math.sqrt(2)
    target=np.array([1.,0.])
    rows=[]
    for exponent in [0,4,8,12]:
        scales=np.array([10.**exponent,10.**(-exponent)])
        H=B*scales
        coeff=np.linalg.solve(B,target)/scales
        physical_columns=H*coeff
        normed=H/np.linalg.norm(H,axis=0)
        rows.append({'column_exponent':exponent,'raw_condition':float(np.linalg.cond(H)),
                     'normalized_condition':float(np.linalg.cond(normed)),
                     'relative_residual':float(np.linalg.norm(physical_columns.sum(axis=1)-target)),
                     'largest_scaled_contribution':float(np.linalg.norm(physical_columns,axis=0).max())})
    return rows


def save_csv(path:Path,rows:list[dict]) -> None:
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    rows=[]
    # These m values validate local identities. They do not satisfy the full base-witness constraints.
    for m in [0.1,0.25,0.5,1.,2.,3.]:
        for lam in [1e-5,1e-3,5e-3,8e-3]:
            for logP in [0.,math.log(1000),math.exp(m)+11.]:
                rows.append(numerical_case(m,lam,logP))
    save_csv(args.output/'quadrature_validation.csv',rows)
    symbolic=symbolic_checks();precision=independent_high_precision()
    unit=cancellation_scaling_cases();save_csv(args.output/'scaling_counterexample.csv',unit)
    budgets=exact_integer_ledger();save_csv(args.output/'corrected_ledger_thresholds.csv',budgets)
    symmetry=max(abs(sigma(x)+sigma(1-x)-1) for x in np.linspace(-.2,1.2,2001))
    summary={'number_of_local_cases':len(rows), 'symbolic_checks':symbolic,
      'sigma_symmetry_max_abs_error':symmetry,
      'max_log_amplitude_error':max(r['pulse_log_amplitude_abs_error'] for r in rows),
      'max_history_log_identity_error':max(max(r['history0_log_identity_abs_error'],r['history1_log_identity_abs_error']) for r in rows),
      'largest_history_ratio_to_proved_bound':max(max(r['history0_ratio_to_analytic_cap'],r['history1_ratio_to_analytic_cap']) for r in rows),
      'high_precision':precision,
      'all_sampled_bounds_pass':all(r['history0_ratio_to_analytic_cap']<=1+1e-10 and r['history1_ratio_to_analytic_cap']<=1+1e-10 for r in rows),
      'status':{'pde_solved':False,'lean_compiled':False,'interval_certified':False,'novelty_established':False}}
    assert all(symbolic.values())
    assert summary['max_log_amplitude_error']<1e-9
    assert summary['max_history_log_identity_error']<1e-9
    assert summary['all_sampled_bounds_pass']
    (args.output/'validation_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    # A plot of a rigorous local bound, not observed fluid or simulated PDE data.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lam=np.logspace(-5,math.log10(.008),300)
    cap0=(math.log(4)-12+60*np.log(lam))/math.log(10)
    cap1=(math.log(4)-12+lam/2+60*(1-lam)*np.log(lam))/math.log(10)
    plt.figure(figsize=(7.2,4.5))
    plt.plot(lam,cap0,label='Amplitude × normalized mass history')
    plt.plot(lam,cap1,linestyle='--',label='Amplitude × normalized angular history')
    plt.xscale('log');plt.xlabel('Pulse parameter lambda');plt.ylabel('log10 of analytic upper bound')
    plt.title('Local history-product bounds: independent of P and m')
    plt.grid(True,alpha=.25);plt.legend();plt.tight_layout()
    plt.savefig(args.output/'history_product_bounds.png',dpi=160);plt.close()
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
