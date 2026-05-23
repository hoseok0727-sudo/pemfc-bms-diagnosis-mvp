from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

from .baseline_model import BaselineVoltageModel


@dataclass(frozen=True)
class DiagnosisSummary:
    state: str
    soh_proxy_pct: float
    residual_z_score: float
    degradation_rate_v_per_h: float
    degradation_speed_v_per_h: float
    estimated_rul_h: float | None
    latest_corrected_voltage_v: float
    eol_voltage_v: float
    rapid_drop_v: float
    rapid_drop_detected: bool


def calculate_diagnosis_features(
    df: pd.DataFrame,
    model: BaselineVoltageModel,
    reference_current_a: float | None = None,
    reference_temperature_c: float | None = None,
) -> pd.DataFrame:
    out = df.copy()
    if reference_current_a is None:
        reference_current_a = float(out["current_a_smooth"].median())
    if reference_temperature_c is None:
        reference_temperature_c = float(out["temperature_c_smooth"].median())

    out["power_w"] = out["voltage_v_smooth"] * out["current_a_smooth"]
    out["voltage_hat_v"] = model.predict(out)
    out["voltage_residual_v"] = out["voltage_v_smooth"] - out["voltage_hat_v"]
    out["residual_z_score"] = (
        out["voltage_residual_v"] - model.residual_mean_
    ) / model.residual_std_

    ref_prediction = model.predict_at_reference(
        len(out),
        reference_current_a=reference_current_a,
        reference_temperature_c=reference_temperature_c,
    )
    out["corrected_voltage_v"] = out["voltage_v_smooth"] - (
        out["voltage_hat_v"] - ref_prediction
    )

    initial_reference = float(out["corrected_voltage_v"].head(max(10, len(out) // 20)).median())
    out["soh_proxy_raw_pct"] = out["corrected_voltage_v"] / initial_reference * 100.0
    out["soh_proxy_pct"] = out["soh_proxy_raw_pct"].clip(upper=100.0)
    out["corrected_voltage_delta_v"] = out["corrected_voltage_v"].diff()
    return out


def estimate_degradation_rate(
    df: pd.DataFrame,
    window_fraction: float = 0.4,
    min_points: int = 20,
) -> float:
    n_window = max(min_points, int(len(df) * window_fraction))
    window = df.tail(min(n_window, len(df))).dropna(subset=["time_h", "corrected_voltage_v"])
    if len(window) < 2:
        return 0.0

    x = window[["time_h"]].to_numpy()
    y = window["corrected_voltage_v"].to_numpy()
    reg = LinearRegression().fit(x, y)
    return float(reg.coef_[0])


def estimate_rul_hours(
    latest_corrected_voltage_v: float,
    degradation_speed_v_per_h: float,
    eol_voltage_v: float,
) -> float | None:
    if degradation_speed_v_per_h <= 0:
        return None
    remaining = (latest_corrected_voltage_v - eol_voltage_v) / degradation_speed_v_per_h
    return float(max(0.0, remaining))


def classify_state(
    soh_proxy_pct: float,
    residual_z_score: float,
    degradation_speed_v_per_h: float,
    rapid_drop_detected: bool = False,
) -> str:
    # The MVP status is driven by the voltage-based SOH proxy and rapid-drop
    # rule. The residual z-score is kept as a reference indicator because it can
    # be overly sensitive when the empirical baseline is imperfect.
    _ = residual_z_score, degradation_speed_v_per_h
    if soh_proxy_pct < 90 or rapid_drop_detected:
        return "Critical"
    if soh_proxy_pct < 95:
        return "Check"
    if soh_proxy_pct < 97:
        return "Warning"
    return "Normal"


def summarize_diagnosis(
    df: pd.DataFrame,
    eol_soh_pct: float = 90.0,
    rapid_drop_threshold_v: float = 0.5,
) -> DiagnosisSummary:
    latest = df.iloc[-1]
    initial_reference = float(df["corrected_voltage_v"].head(max(10, len(df) // 20)).median())
    eol_voltage_v = initial_reference * eol_soh_pct / 100.0
    degradation_rate = estimate_degradation_rate(df)
    degradation_speed = max(0.0, -degradation_rate)
    rul_h = estimate_rul_hours(
        latest_corrected_voltage_v=float(latest["corrected_voltage_v"]),
        degradation_speed_v_per_h=degradation_speed,
        eol_voltage_v=eol_voltage_v,
    )
    previous = df["corrected_voltage_v"].iloc[-20] if len(df) >= 20 else np.nan
    rapid_drop_v = 0.0
    rapid_drop_detected = False
    if not np.isnan(previous):
        rapid_drop_v = float(latest["corrected_voltage_v"] - previous)
        rapid_drop_detected = rapid_drop_v < -abs(rapid_drop_threshold_v)
    state = classify_state(
        soh_proxy_pct=float(latest["soh_proxy_pct"]),
        residual_z_score=float(latest["residual_z_score"]),
        degradation_speed_v_per_h=degradation_speed,
        rapid_drop_detected=rapid_drop_detected,
    )
    return DiagnosisSummary(
        state=state,
        soh_proxy_pct=float(latest["soh_proxy_pct"]),
        residual_z_score=float(latest["residual_z_score"]),
        degradation_rate_v_per_h=degradation_rate,
        degradation_speed_v_per_h=degradation_speed,
        estimated_rul_h=rul_h,
        latest_corrected_voltage_v=float(latest["corrected_voltage_v"]),
        eol_voltage_v=eol_voltage_v,
        rapid_drop_v=rapid_drop_v,
        rapid_drop_detected=rapid_drop_detected,
    )
