"""Smooth near-minimal logarithmic drop under |d k/d log(y)| <= delta.

This only replaces one scalar transition. It does NOT establish a new global
Navier--Stokes profile, compatibility with other lemmas, or lower PDE cost.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
import numpy as np
from scipy.integrate import quad
from pulse_audit import sigma, integral_sigma


def sigma_derivative(x:float)->float:
    if not 0 < x < 1:
        return 0.
    # Both factors are evaluated separately to avoid 1-sigma cancellation.
    return sigma(x)*sigma(1-x)*(2/(1-x)**3+2/x**3)


@dataclass(frozen=True)
class SmoothDrop:
    delta:float
    shoulder:float=1.0
    height:float=4.0

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.delta,self.shoulder,self.height)):
            raise ValueError('all parameters must be finite')
        if min(self.delta,self.shoulder,self.height)<=0:
            raise ValueError('all parameters must be positive')
        if self.shoulder>self.height/self.delta:
            raise ValueError('shoulder must not exceed height/delta')

    @property
    def width(self)->float:
        return self.height/self.delta+self.shoulder

    def rate(self,x:float)->float:
        return self.delta*sigma(x/self.shoulder)*sigma((self.width-x)/self.shoulder)

    def accumulated(self,x:float)->float:
        a=self.shoulder;L=self.width
        if x<=0:return 0.
        if x>=L:return L-a
        if x<a:return a*integral_sigma(x/a)
        if x<=L-a:return x-a/2
        return L-a-a*integral_sigma((L-x)/a)

    def value(self,x:float)->float:
        if x<=0:return self.height
        if x>=self.width:return 0.
        return self.height-self.delta*self.accumulated(x)


def run(output:Path)->dict:
    output.mkdir(parents=True,exist_ok=True)
    xs=np.linspace(0,1,100001)
    peak=max(sigma_derivative(float(x)) for x in xs)
    rows=[]
    for delta in [1/8,1/16,1/32]:
        new=SmoothDrop(delta)
        grid=np.linspace(-1,new.width+1,4001)
        values=np.array([new.value(float(x)) for x in grid])
        rates=np.array([new.rate(float(x)) for x in grid])
        area=sum(quad(new.rate,a,b,epsabs=1e-12,epsrel=1e-12)[0]
                 for a,b in [(0,1),(1,new.width-1),(new.width-1,new.width)])
        old_width=32/delta
        row={'delta':delta,'universal_infimum_log_width':4/delta,
             'old_sigma_family_min_log_width':old_width,
             'new_C_infinity_log_width':new.width,
             'width_reduction_factor':old_width/new.width,
             'new_max_rate':float(rates.max()),'integrated_drop_error':abs(area-4),
             'new_min_value':float(values.min()),'new_max_value':float(values.max()),
             'monotonicity_largest_upward_step':float(np.diff(values).max()),
             'log10_old_clock_end':old_width/math.log(10),
             'log10_new_clock_end':new.width/math.log(10),
             'new_log_coordinate_second_derivative_bound':8*delta/new.shoulder}
        assert row['integrated_drop_error']<1e-10
        assert row['new_max_rate']<=delta+1e-14
        assert row['new_min_value']>=-1e-12 and row['new_max_value']<=4+1e-12
        assert row['monotonicity_largest_upward_step']<=1e-12
        rows.append(row)
    with (output/'drop_profile_validation.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=1/8;new=SmoothDrop(d);old_width=32/d
    x=np.linspace(0,old_width,2001)
    plt.figure(figsize=(7.4,4.5))
    plt.plot(x,[4*sigma(1-z/old_width) for z in x],label='Original step family, width 256')
    plt.plot(x,[new.value(float(z)) for z in x],label='Tapered-rate alternative, width 33')
    plt.xlabel('x = log(y), a mathematical coordinate');plt.ylabel('Drop coefficient k')
    plt.title('Same first-derivative cap: |dk/dx| <= 1/8')
    plt.grid(True,alpha=.25);plt.legend();plt.tight_layout()
    plt.savefig(output/'drop_comparison.png',dpi=160);plt.close()
    summary={'sigma_derivative_max_on_100001_grid':peak,
             'sigma_derivative_at_one_half':sigma_derivative(.5),
             'analytic_global_sigma_derivative_cap':8,
             'new_profiles':rows,
             'limitations':['Only a scalar first-derivative constraint was enforced',
                 'Higher derivatives change and may worsen',
                 'Original profile-dependent lemmas need re-derivation',
                 'No full PDE, cone, energy closure, or physical bandwidth verification']}
    (output/'drop_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=Path('results'))
    run(p.parse_args().output)
