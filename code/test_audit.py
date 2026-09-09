"""Automated checks for local identities, not a theorem-prover certificate."""
import math
from fractions import Fraction
import pytest
import numpy as np
from pulse_audit import (sigma, symbolic_checks, numerical_case, exact_integer_ledger,
                         cancellation_scaling_cases)
from drop_profile import SmoothDrop, sigma_derivative


def test_exact_algebra():
    assert all(symbolic_checks().values())


def test_sigma_symmetry_and_midpoint():
    assert sigma(0.5)==0.5
    assert sigma_derivative(0.5)==8.0
    for x in np.linspace(-1,2,301):
        assert abs(sigma(float(x))+sigma(float(1-x))-1)<1e-13


@pytest.mark.parametrize('m',[0.1,0.5,1.,2.,3.])
@pytest.mark.parametrize('lam',[1e-5,.001,.008])
def test_local_pulse_bounds(m,lam):
    row=numerical_case(m,lam,math.exp(m)+11)
    assert row['pulse_log_amplitude_abs_error']<1e-9
    assert row['history0_log_identity_abs_error']<1e-9
    assert row['history1_log_identity_abs_error']<1e-9
    assert row['history0_ratio_to_analytic_cap']<=1+1e-10
    assert row['history1_ratio_to_analytic_cap']<=1+1e-10


@pytest.mark.parametrize('delta',[1/8,1/16,1/32])
def test_smooth_drop(delta):
    f=SmoothDrop(delta)
    assert f.value(0)==4 and f.value(f.width)==0
    assert f.rate(0)==0 and f.rate(f.width)==0
    assert abs(f.accumulated(f.width)*delta-4)<1e-13
    for x in np.linspace(0,f.width,501):
        assert 0<=f.rate(float(x))<=delta
        assert -1e-12<=f.value(float(x))<=4+1e-12


def test_ledger_strict_integer_threshold():
    rows=exact_integer_ledger()
    assert rows[0]['first_J_with_positive_ledger_rate']==25014
    assert rows[1]['first_J_with_positive_ledger_rate']==55034
    for row in rows:
        h=Fraction(1,1000);m=row['derivative_order'];j=row['first_J_with_positive_ledger_rate']
        def rate(J):
            return h*(Fraction(J,10)+Fraction(7,10))-(m*(m+2)+Fraction(5,2)+2*h+2*h*m)
        assert rate(j)>0 and rate(j-1)<=0


def test_raw_condition_not_scale_invariant():
    rows=cancellation_scaling_cases()
    assert rows[-1]['raw_condition']>1e23
    assert abs(rows[-1]['normalized_condition']-1)<1e-12
    assert rows[-1]['relative_residual']<1e-12


def test_invalid_inputs():
    with pytest.raises(ValueError): SmoothDrop(-1)
    with pytest.raises(ValueError): SmoothDrop(.5,shoulder=100)
    with pytest.raises(ValueError): numerical_case(-1,.001,1)
