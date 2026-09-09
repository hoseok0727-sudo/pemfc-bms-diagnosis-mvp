"""Exact rational audit of a conservative full-pulse scalar example.

This is an arithmetic certificate for inequalities in the handwritten proof,
not a Lean proof or an interval enclosure of the full Navier--Stokes solution.
No floating-point arithmetic is used in the asserted inequalities.
"""
from fractions import Fraction as F
import json

def exp_lower(x:F,n:int=64)->F:
    if x<0 or n<0:raise ValueError('x and n must be nonnegative')
    t=F(1);s=t
    for k in range(1,n+1):t=t*x/k;s+=t
    return s

def exp_upper(x:F,n:int=64)->F:
    if not 0<=x<n+2:raise ValueError('Taylor tail ratio must be below 1')
    t=F(1);s=t
    for k in range(1,n+1):t=t*x/k;s+=t
    nxt=t*x/(n+1)
    return s+nxt/(1-x/F(n+2))

def certificate()->dict:
    l=F(1,1000);h=F(1,100000000);alpha=F(1,2)+l;c=1-l
    y0=F(1,100)/l;k=1/(2500*l*l)
    # First envelope is increasing through y0: its unique positive
    # derivative root is no smaller than y0 (see manuscript).
    assert alpha*y0**3-3*y0*y0-2*k<0
    capA=F(1250)*(l*y0)**3/exp_lower(alpha*y0)
    capB=(l/alpha)/exp_lower(1+alpha*F(7,800)/l)
    G=F(9163,10**9)
    assert capA<G and capB<G
    # Constants used for the shape-independent moment-repair estimate.
    moment_floor=F(111,1840)
    assert moment_floor>F(1,20)
    assert exp_upper(F(2))<8
    assert exp_upper(F(61,120)*F(63,20))<5
    # e>2 -> exp(-N)<2**(-N); all rational upper bounds below.
    e0=F(1,10**176)
    dip=F(4,10**183)
    T0=F(4,10**177*2**12987)
    Tm=F(23,2**7500)
    a0=F(6,5);a1=F(1)
    R0=800/l*(T0+a0*Tm)
    R1=800/l*(T0+(a1+2*a0)*Tm)
    repair=(3+l)*R0+l*R1
    assert repair<F(1,10**2000)
    main=G*(3*a0+l*(3*a0+a1))
    qlo=((l-h)-main-repair)/c-e0-dip
    claim=F(9679,10**7)
    assert qlo>claim
    return {'all_exact_checks_pass':True,'arithmetic':'fractions.Fraction, no floating point in assertions',
            'lambda':'1/1000','h':'1/100000000','B':'1',
            'G_certified_upper':'9163/1000000000',
            'repair_source_certified_upper':'10^(-2000)',
            'Q_strict_lower':'9679/10000000',
            'claim_scope':'conditional scalar angular lag on 0<=y<=13000',
            'not_a_Lean_or_PDE_certificate':True}

if __name__=='__main__':print(json.dumps(certificate(),indent=2))
