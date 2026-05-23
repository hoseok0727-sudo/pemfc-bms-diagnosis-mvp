from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


FEATURE_COLUMNS = [
    "current_a_smooth",
    "current_a_sq",
    "temperature_c_smooth",
    "current_temp",
]


class BaselineVoltageModel:
    """Empirical PEMFC normal-voltage baseline model."""

    def __init__(self) -> None:
        self.model = LinearRegression()
        self.residual_mean_: float | None = None
        self.residual_std_: float | None = None
        self.reference_current_a_: float | None = None
        self.reference_temperature_c_: float | None = None
        self.baseline_voltage_v_: float | None = None
        self.use_constant_baseline_ = False
        self.is_fitted = False

    def fit(
        self,
        df: pd.DataFrame,
        normal_fraction: float = 0.1,
        normal_hours: float | None = 24.0,
    ) -> "BaselineVoltageModel":
        train = initial_normal_section(
            df,
            normal_fraction=normal_fraction,
            normal_hours=normal_hours,
        )
        self.reference_current_a_ = float(train["current_a_smooth"].median())
        self.reference_temperature_c_ = float(train["temperature_c_smooth"].median())
        self.baseline_voltage_v_ = float(train["voltage_v_smooth"].median())
        self.use_constant_baseline_ = _has_low_excitation(train)

        if self.use_constant_baseline_:
            x_train = pd.DataFrame({"constant_feature": np.zeros(len(train))})
        else:
            x_train = build_features(
                train,
                reference_current_a=self.reference_current_a_,
                reference_temperature_c=self.reference_temperature_c_,
            )
        y_train = train["voltage_v_smooth"].to_numpy()
        self.model.fit(x_train, y_train)

        residual = y_train - self.model.predict(x_train)
        self.residual_mean_ = float(np.mean(residual))
        residual_std_floor = max(0.01 * abs(self.baseline_voltage_v_), 1e-6)
        self.residual_std_ = float(max(np.std(residual, ddof=1), residual_std_floor))
        self.is_fitted = True
        return self

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        self._check_fitted()
        if self.use_constant_baseline_:
            return np.full(len(df), self.baseline_voltage_v_)
        return self.model.predict(
            build_features(
                df,
                reference_current_a=self.reference_current_a_,
                reference_temperature_c=self.reference_temperature_c_,
            )
        )

    def predict_at_reference(
        self,
        n_rows: int,
        reference_current_a: float,
        reference_temperature_c: float,
    ) -> np.ndarray:
        self._check_fitted()
        if self.use_constant_baseline_:
            return np.full(n_rows, self.baseline_voltage_v_)
        ref = pd.DataFrame(
            {
                "current_a_smooth": np.full(n_rows, reference_current_a),
                "temperature_c_smooth": np.full(n_rows, reference_temperature_c),
            }
        )
        return self.model.predict(
            build_features(
                ref,
                reference_current_a=self.reference_current_a_,
                reference_temperature_c=self.reference_temperature_c_,
            )
        )

    def coefficients(self) -> dict[str, float]:
        self._check_fitted()
        values = {"beta0": float(self.model.intercept_)}
        values.update(
            {
                f"beta{i + 1}": float(coef)
                for i, coef in enumerate(self.model.coef_)
            }
        )
        return values

    def _check_fitted(self) -> None:
        if not self.is_fitted:
            raise RuntimeError("BaselineVoltageModel must be fitted before prediction.")


def initial_normal_section(
    df: pd.DataFrame,
    normal_fraction: float = 0.1,
    normal_hours: float | None = 24.0,
) -> pd.DataFrame:
    if not 0 < normal_fraction <= 1:
        raise ValueError("normal_fraction must be in (0, 1].")
    fraction_rows = max(20, int(len(df) * normal_fraction))
    if normal_hours is not None and "time_h" in df.columns:
        start_time = float(df["time_h"].min())
        hour_section = df[df["time_h"] <= start_time + normal_hours]
        if len(hour_section) >= 20:
            fraction_section = df.head(min(fraction_rows, len(df)))
            if len(hour_section) <= len(fraction_section):
                return hour_section.copy()
    return df.head(min(fraction_rows, len(df))).copy()


def _has_low_excitation(df: pd.DataFrame) -> bool:
    """Detect fixed-load aging data where polynomial fitting is ill-conditioned."""
    current_std = float(df["current_a_smooth"].std(ddof=1) or 0.0)
    temperature_std = float(df["temperature_c_smooth"].std(ddof=1) or 0.0)
    current_span = float(df["current_a_smooth"].max() - df["current_a_smooth"].min())
    temperature_span = float(df["temperature_c_smooth"].max() - df["temperature_c_smooth"].min())
    return current_std < 0.5 and temperature_std < 0.5 and current_span < 2.0 and temperature_span < 2.0


def build_features(
    df: pd.DataFrame,
    reference_current_a: float | None = None,
    reference_temperature_c: float | None = None,
) -> pd.DataFrame:
    current = df["current_a_smooth"].to_numpy()
    temperature = df["temperature_c_smooth"].to_numpy()
    if reference_current_a is None:
        reference_current_a = float(np.median(current))
    if reference_temperature_c is None:
        reference_temperature_c = float(np.median(temperature))
    current_delta = current - reference_current_a
    temperature_delta = temperature - reference_temperature_c
    return pd.DataFrame(
        {
            "current_a_smooth": current_delta,
            "current_a_sq": current_delta**2,
            "temperature_c_smooth": temperature_delta,
            "current_temp": current_delta * temperature_delta,
        },
        index=df.index,
    )
