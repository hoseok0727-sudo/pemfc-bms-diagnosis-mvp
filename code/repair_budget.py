"""v0.15: full-interval scalar angular-positivity bound with two-moment repair.

The independent variable is a logarithmic radius, not physical time.
The analytic bound uses only plateau/support/range facts and applies to the
source's selected cutoff. Numerical experiments use explicit representative
cutoffs, not an unproved identification with Mathlib's Nonempty.some choice.
This is NOT a Navier--Stokes solver or a verification of the full construction.
"""
from __future__ import annotations
import argparse, csv, json, math, platform, sys
from fractions import Fraction as F
from pathlib import Path
from typing import Any
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import expit, logsumexp
import mpmath as mp
from angular_budget import gaussian_envelope, log_main
from pulse_audit import numerical_case

mp.iv.dps = 70
S=F(3,20)
CENTERS=np.array([3.0,1.0])

def rat(x:Any)->F:
    return x if isinstance(x,F) else F(str(x))

def iv(x:Any):
    r=rat(x); return mp.iv.mpf(r.numerator)/mp.iv.mpf(r.denominator)

def upper(x)->float:
    return math.nextafter(math.nextafter(float(x.b),math.inf),math.inf)

def lower(x)->float:
    return math.nextafter(math.nextafter(float(x.a),-math.inf),-math.inf)

def log10_interval(x)->dict[str,float]:
    if float(x.a)<=0 and x.a==0: raise ValueError('positive interval needed')
    a=mp.iv.log(x)/mp.iv.log(iv(10))
    return {'lower':lower(a),'upper':upper(a)}

def validate(lam:Any,h:Any=0,B:Any=1)->tuple[F,F,F]:
    l,hh,b=map(rat,(lam,h,B))
    if not 0<l<=F(1,120): raise ValueError('require 0<lambda<=1/120')
    if not 0<=hh<l: raise ValueError('require 0<=h<lambda')
    if b<0: raise ValueError('B must be nonnegative')
    return l,hh,b

def geometric_certificates()->dict[str,Any]:
    s=iv(S); co=(mp.iv.exp(s)+mp.iv.exp(-s))/2
    si=(mp.iv.exp(s)-mp.iv.exp(-s))/2
    plateau=mp.iv.log((co+si/4)/(co-si/4))
    moment_floor=mp.iv.exp(-s/2)*plateau
    inv_const=mp.iv.exp(iv(2))
    field_const=mp.iv.exp(iv(F(61,120)*F(63,20)))
    return {'moment_floor_lower':lower(moment_floor),'chosen_moment_floor':.05,
            'inverse_factor_upper':upper(inv_const),'chosen_inverse_factor':8.,
            'field_weight_factor_upper':upper(field_const),'chosen_field_weight_factor':5.,
            'moment_floor_pass':bool(moment_floor.a>iv(F(1,20)).b),
            'inverse_factor_pass':bool(inv_const.b<iv(8).a),
            'field_weight_pass':bool(field_const.b<iv(5).a),
            'interval_decimal_digits':mp.iv.dps}

def repair_bounds(lam:Any,B:Any=1,a0:Any=F(6,5),a1:Any=1):
    l,_,b=validate(lam,0,B); aa,ad=map(rat,(a0,a1))
    if aa<0 or ad<0: raise ValueError('amplitude bounds must be nonnegative')
    c=1-l; L=13/l
    # Debt families AFTER multiplication by the true background amplitude.
    T0=4*mp.iv.exp(-12+iv(l/2)+iv(60*c)*mp.iv.log(iv(l))-iv(c*L))
    Tm=23*iv(b)*mp.iv.exp(-iv(F(15,2)/l)-9)
    K=iv(800/l)
    f0=K*(T0+iv(aa)*Tm)
    f1=K*(T0+iv(ad+2*aa)*Tm)
    source=iv(3+l)*f0+iv(l)*f1
    return T0,Tm,f0,f1,source

def full_certificate(lam:Any,h:Any,B:Any=1,a0:Any=F(6,5),a1:Any=1)->dict[str,Any]:
    l,hh,b=validate(lam,h,B); aa,ad=map(rat,(a0,a1)); c=1-l
    env=gaussian_envelope(l)
    # Preserve the binary float upper endpoint exactly when importing it.
    G=iv(F.from_float(env['G_sup_upper']))
    e0=10*mp.iv.exp(iv(60*c)*mp.iv.log(iv(l)))
    dip=4*mp.iv.exp(-12+iv(61)*mp.iv.log(iv(l))+iv(c/l)*mp.iv.log(iv(c)))
    main=iv(b)*G*iv(3*aa+l*(3*aa+ad))
    _,_,f0,f1,rep=repair_bounds(l,b,aa,ad)
    available=iv(l-hh)-iv(c)*(e0+dip)-main-rep
    qlow=available/iv(c);eq=iv((l-hh)/c)
    # Solve the affine bound in B for a strict sufficient threshold.
    _,_,_,_,rep0=repair_bounds(l,0,aa,ad)
    _,_,_,_,rep1=repair_bounds(l,1,aa,ad)
    # No cancellation of exponentially small interval quantities:
    rep_slope=iv(800/l)*23*mp.iv.exp(-iv(F(15,2)/l)-9)*iv((3+l)*aa+l*(ad+2*aa))
    limit=(iv(l-hh)-iv(c)*(e0+dip)-rep0)/(G*iv(3*aa+l*(3*aa+ad))+rep_slope)
    return {'lambda':float(l),'h':float(hh),'B':float(b),
      'main_G_sup_upper':env['G_sup_upper'],
      'repair_Y_log10_upper':log10_interval(f0)['upper'],
      'repair_detaY_log10_upper':log10_interval(f1)['upper'],
      'repair_source_log10_upper':log10_interval(rep)['upper'],
      'main_source_upper':upper(main),
      'full_Q_lower':lower(qlow),'min_equilibrium_Q':float((l-hh)/c),
      'fraction_min_equilibrium_lower':lower(qlow/eq),
      'B_strict_sufficient_limit_lower':lower(limit),
      'positive_certified':bool(available.a>0),
      'scope':'scalar angular lag on the entire pulse only; conditional input bounds',
      'interval_decimal_digits':mp.iv.dps}

def smooth_transition(x:float,power:int=1)->float:
    if x<=0:return 0.
    if x>=1:return 1.
    return float(expit(1/(1-x)**power-1/x**power))

def cutoff(x:float,power:int=1)->float:
    return smooth_transition(2*(1-abs(x)),power)

def bump_shape(delta:float,power:int=1)->float:
    x=2*(math.exp(delta)-math.cosh(float(S)))/math.sinh(float(S))
    return cutoff(x,power)

def bump_breakpoints()->list[float]:
    s=float(S);co=math.cosh(s);si=math.sinh(s)
    return [math.log(co+q*si) for q in [-.5,-.25,.25,.5]]

def moment(beta:float,power:int=1,tol:float=2e-13)->float:
    points=[-float(S)]+bump_breakpoints()+[float(S)]
    return sum(quad(lambda t:math.exp(beta*t)*bump_shape(t,power),x,y,
                    epsabs=tol,epsrel=tol,limit=150)[0] for x,y in zip(points,points[1:]))

def V_matrix(lam:float)->np.ndarray:
    betas=np.array([.5-lam,.5-2*lam])
    return np.exp(-betas[:,None]*CENTERS[None,:])

def stable_V_solve(lam:float,rhs:np.ndarray)->np.ndarray:
    t0,t1=map(float,rhs)
    denominator=math.exp(lam)*math.expm1(2*lam)
    d0=(t1-math.exp(lam)*t0)/denominator
    return np.exp((.5-lam)*CENTERS)*np.array([d0,t0-d0])

def log_main_moment(lam:float,beta:float,tol:float=2e-12)->float:
    """Log integral exp(beta*y) f(lambda*y), without huge exponentiation."""
    fun=lambda z: beta*(z-11)/lam+log_main(z)
    opt=minimize_scalar(lambda z:-fun(z),bounds=(10.,11.-1e-10),method='bounded',
                         options={'xatol':1e-12})
    peak=float(opt.x);shift=fun(peak)
    pts=sorted(set([0.,.02,10.,max(10.,peak-(11-peak)*3),peak,11.]))
    q=sum(quad(lambda z:0. if not 0<z<11 else math.exp(fun(z)-shift),x,y,
               epsabs=tol,epsrel=tol,limit=200)[0] for x,y in zip(pts,pts[1:]))
    if q<=0:raise ArithmeticError('nonpositive moment')
    return beta*11/lam+shift+math.log(q)-math.log(lam)

def numerical_repair(lam:float,m:float=1.,B:float=1.,power:int=1)->dict[str,Any]:
    validate(lam,0,B)
    if B<=0:raise ValueError('numerical sample requires B>0')
    H=math.exp(m)+12;logP=math.log(B)-.3+lam/2+H/2-(30+60*lam)*math.log(lam)
    old=numerical_case(m,lam,logP)
    bq=np.array([old['history0_log_product'],old['history1_log_product']])
    beta=np.array([.5-lam,.5-2*lam]);gamma=.5+lam+beta;L=13/lam
    M=np.array([moment(t,power) for t in beta])
    hm=np.array([log_main_moment(lam,t) for t in beta])
    mats=V_matrix(lam);logsP=bq-gamma*L-np.log(M)
    logsA=math.log(B)+hm-gamma*L-np.log(M)
    outputs={}; residual=0.; hp_error=0.
    with mp.workdps(80):
        lv=mp.mpf(str(lam));vv=mp.matrix([[mp.exp(-a*(mp.mpf('.5')-(i+1)*lv)) for a in (3,1)] for i in (0,1)])
        for key,logs in [('prefix',logsP),('main',logsA)]:
            scale=float(max(logs));rhs=-np.exp(logs-scale);sol=stable_V_solve(lam,rhs)
            residual=max(residual,float(np.linalg.norm(mats@sol-rhs,np.inf)/np.linalg.norm(rhs,np.inf)))
            hp=mp.lu_solve(vv,mp.matrix([mp.mpf(str(t)) for t in rhs]))
            hp_float=np.array([float(x) for x in hp]);hp_error=max(hp_error,float(np.linalg.norm(sol-hp_float,np.inf)/np.linalg.norm(hp_float,np.inf)))
            delta=np.linspace(-float(S),float(S),3001)
            weights=np.array([math.exp(-(.5+lam)*t)*bump_shape(t,power) for t in delta])
            effect=np.max(np.abs(sol)*np.exp((.5+lam)*CENTERS))*max(weights)
            outputs[key+'_weighted_log10_sampled']=float((scale+math.log(effect))/math.log(10))
    return {'lambda':lam,'m':m,'B':B,'representative_power':power,
      'small_moment_min':float(min(M)),'V_condition_2':float(np.linalg.cond(mats)),
      'scaled_moment_residual':residual,'V_solve_vs_80digit_relative_error':hp_error,
      **outputs,'cutoff_note':'explicit representative; uniform analytic bounds cover the selected source cutoff'}

def save_csv(path:Path,rows:list[dict]):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run(out:Path):
    out.mkdir(parents=True,exist_ok=True)
    geo=geometric_certificates();assert all(geo[k] for k in ('moment_floor_pass','inverse_factor_pass','field_weight_pass'))
    rows=[full_certificate(l,l/F(100000),1) for l in [F(1,100000),F(1,10000),F(1,1000),F(1,200),F(1,125),F(1,120)]]
    save_csv(out/'full_pulse_certificates.csv',rows)
    numerical=[]
    for l in [.00001,.0001,.001,.005,.008,1/120]:
        for power in (1,2):
            numerical.append(numerical_repair(l,power=power))
    save_csv(out/'normalized_repair_validation.csv',numerical)
    params=[full_certificate(l,l/F(100000),b) for l in [F(1,100000),F(1,10000),F(1,1000),F(1,200),F(1,125)] for b in [F(1,10),F(1,2),1,10,30]]
    save_csv(out/'local_parameter_screen.csv',params)
    summary={'geometry':geo,'full_pulse_scalar_example':full_certificate(F(1,1000),F(1,100000000),1),
      'representative_cutoff_cases':len(numerical),'max_scaled_moment_residual':max(x['scaled_moment_residual'] for x in numerical),
      'max_solve_vs_80digit_error':max(x['V_solve_vs_80digit_relative_error'] for x in numerical),
      'local_parameter_cases':len(params),'sufficient_condition_passes':sum(x['positive_certified'] for x in params),
      'not_verified':['full PDE','full stress cone','finite forcing realizability','Lean kernel','academic novelty'],
      'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__}
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path('results'))
    run(p.parse_args().out)
