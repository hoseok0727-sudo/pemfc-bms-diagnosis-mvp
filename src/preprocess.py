from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ["time_h", "voltage_v", "current_a", "temperature_c"]
OPTIONAL_COLUMNS = ["source_file"]


def generate_synthetic_sample(
    output_path: str | Path,
    n_samples: int = 1200,
    seed: int = 42,
) -> Path:
    """Create a simple PEMFC-like aging dataset for first-run demos."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    time_h = np.linspace(0, 1000, n_samples)
    current_a = 46 + 13 * np.sin(time_h / 45) + rng.normal(0, 1.8, n_samples)
    current_a = np.clip(current_a, 18, 75)
    temperature_c = 63 + 3.5 * np.sin(time_h / 85 + 0.4) + rng.normal(0, 0.45, n_samples)

    baseline_voltage = (
        31.8
        - 0.072 * current_a
        - 0.00042 * current_a**2
        + 0.035 * (temperature_c - 60)
        - 0.00032 * current_a * (temperature_c - 60)
    )
    degradation = -0.0041 * time_h
    voltage_v = baseline_voltage + degradation + rng.normal(0, 0.09, n_samples)

    df = pd.DataFrame(
        {
            "time_h": time_h,
            "voltage_v": voltage_v,
            "current_a": current_a,
            "temperature_c": temperature_c,
        }
    )
    df.to_csv(output_path, index=False)
    return output_path


def load_csv(path: str | Path) -> pd.DataFrame:
    """Load a standard CSV file and validate MVP columns."""
    df = pd.read_csv(path)
    return validate_required_columns(df)


def preprocess_data(
    df: pd.DataFrame,
    rolling_window: int = 15,
    clip_quantiles: tuple[float, float] = (0.001, 0.999),
) -> pd.DataFrame:
    """Clean, sort, smooth, and lightly winsorize BMS operating data."""
    df = validate_required_columns(df).copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=REQUIRED_COLUMNS)

    for col in REQUIRED_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=REQUIRED_COLUMNS)
    df = df.sort_values("time_h").drop_duplicates(subset=["time_h"]).reset_index(drop=True)

    for col in ["voltage_v", "current_a", "temperature_c"]:
        low, high = df[col].quantile(list(clip_quantiles))
        df[col] = df[col].clip(low, high)
        smooth_col = f"{col}_smooth"
        df[smooth_col] = (
            df[col]
            .rolling(window=rolling_window, min_periods=1, center=True)
            .mean()
        )

    return df


def validate_required_columns(df: pd.DataFrame) -> pd.DataFrame:
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}. "
            f"Required schema is {REQUIRED_COLUMNS}."
        )
    return df


def load_ieee_phm2014_zip(
    zip_path: str | Path,
    fuel_cell: str = "FC1",
    max_rows_per_file: int | None = None,
) -> pd.DataFrame:
    """Load IEEE PHM 2014 PEMFC aging Excel files from the released zip.

    Only Ageing workbooks are used. EIS, polarization, HFR, Rct, and cell-level
    diagnostic files are intentionally ignored for the BMS-level MVP scope.
    """
    zip_path = Path(zip_path)
    if not zip_path.exists():
        raise FileNotFoundError(zip_path)

    fuel_cell = fuel_cell.upper()
    rows: list[pd.DataFrame] = []
    with zipfile.ZipFile(zip_path) as archive:
        names = _ageing_workbook_names(archive.namelist(), fuel_cell)
        if not names:
            raise ValueError(f"No {fuel_cell} Ageing workbooks found in {zip_path}.")

        for name in names:
            with archive.open(name) as file_obj:
                content = file_obj.read()
            frames = _read_excel_workbook(content, source_name=name)
            for frame in frames:
                mapped = map_bms_columns(frame, source_file=name)
                if mapped is None:
                    continue
                if max_rows_per_file is not None:
                    mapped = mapped.head(max_rows_per_file)
                rows.append(mapped)

    if not rows:
        raise ValueError(
            "Could not map IEEE workbook columns to time_h, voltage_v, "
            "current_a, temperature_c. Inspect the workbook headers and update "
            "COLUMN_ALIASES in preprocess.py."
        )

    df = pd.concat(rows, ignore_index=True)
    df = _normalize_time_axis(df)
    return validate_required_columns(df)


def map_bms_columns(df: pd.DataFrame, source_file: str = "") -> pd.DataFrame | None:
    """Map common IEEE/PEMFC header variants to the MVP schema."""
    if df.empty:
        return None

    cleaned = df.copy()
    cleaned.columns = [str(col).strip() for col in cleaned.columns]
    normalized_to_original = {_normalize_header(col): col for col in cleaned.columns}

    mapping = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        match = _find_column(normalized_to_original, aliases)
        if match is None:
            return None
        mapping[canonical] = match

    out = pd.DataFrame({target: cleaned[source] for target, source in mapping.items()})
    for col in REQUIRED_COLUMNS:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=REQUIRED_COLUMNS)
    if out.empty:
        return None
    out["source_file"] = source_file
    return out


COLUMN_ALIASES: dict[str, Iterable[str]] = {
    "time_h": [
        "time_h",
        "timeh",
        "time",
        "time_s",
        "times",
        "t",
        "hour",
        "hours",
        "operating_time",
        "operatingtime",
    ],
    "voltage_v": [
        "voltage_v",
        "voltage",
        "vstack",
        "stackvoltage",
        "stack_voltage",
        "u",
        "utot",
        "u_tot",
        "fuelcellvoltage",
        "fcvoltage",
    ],
    "current_a": [
        "current_a",
        "current",
        "istack",
        "stackcurrent",
        "stack_current",
        "ia",
        "i",
        "loadcurrent",
        "fc_current",
    ],
    "temperature_c": [
        "temperature_c",
        "temperature",
        "temp",
        "tstack",
        "stacktemperature",
        "stack_temperature",
        "tinwat",
        "toutwat",
        "tinwatc",
        "toutwatc",
        "coolantin",
        "coolantout",
        "waterin",
        "waterout",
        "tin",
        "tout",
        "tc",
        "coolanttemperature",
    ],
}


def _ageing_workbook_names(names: list[str], fuel_cell: str) -> list[str]:
    return sorted(
        name
        for name in names
        if Path(name).name.upper().startswith(f"{fuel_cell}_AGEING")
        and "AGEING" in name.upper()
        and name.lower().endswith((".xlsx", ".xls"))
    )


def _read_excel_workbook(content: bytes, source_name: str) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    workbook = pd.ExcelFile(io.BytesIO(content))
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(workbook, sheet_name=sheet)
        if not frame.empty:
            frame["__sheet_name"] = sheet
            frame["__source_name"] = source_name
            frames.append(frame)
    return frames


def _find_column(normalized_to_original: dict[str, str], aliases: Iterable[str]) -> str | None:
    alias_list = [_normalize_header(alias) for alias in aliases]
    for alias in alias_list:
        for normalized, original in normalized_to_original.items():
            if normalized == alias:
                return original
    for alias in alias_list:
        if len(alias) < 3:
            continue
        for normalized, original in normalized_to_original.items():
            if alias in normalized:
                return original
    for normalized, original in normalized_to_original.items():
        if normalized in alias_list:
            return original
    return None


def _normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _normalize_time_axis(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    time = pd.to_numeric(df["time_h"], errors="coerce")
    if time.dropna().empty:
        return df

    # Many raw files use seconds. Convert to hours when the scale is clearly too large.
    span = time.max() - time.min()
    if span > 10_000:
        time = (time - time.min()) / 3600.0
    else:
        time = time - time.min()
    df["time_h"] = time
    return df.sort_values("time_h").reset_index(drop=True)
