from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATASETS = {
    "normal_data": PROJECT_ROOT / "outputs" / "normal_data_run" / "diagnosis_timeseries.csv",
    "synthetic": PROJECT_ROOT / "outputs" / "synthetic_run" / "diagnosis_timeseries.csv",
    "ieee_fc1": PROJECT_ROOT / "outputs" / "ieee_fc1_run" / "diagnosis_timeseries.csv",
    "ieee_fc2": PROJECT_ROOT / "outputs" / "ieee_fc2_run" / "diagnosis_timeseries.csv",
    "bad_data": PROJECT_ROOT / "outputs" / "bad_data_run" / "diagnosis_timeseries.csv",
}


def export_dashboard_data(
    output_path: str | Path = PROJECT_ROOT / "ui" / "dashboard_data.json",
    max_points: int = 1600,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "generated_from": "PEMFC BMS early anomaly detection MVP outputs",
        "datasets": [],
    }
    for dataset_id, csv_path in DEFAULT_DATASETS.items():
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        report_path = csv_path.with_name("diagnosis_report.txt")
        payload["datasets"].append(
            {
                "id": dataset_id,
                "label": _dataset_label(dataset_id),
                "source": str(csv_path.relative_to(PROJECT_ROOT)),
                "summary": _parse_report(report_path),
                "points": _compact_points(df, max_points=max_points),
            }
        )

    output_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return output_path


def _compact_points(df: pd.DataFrame, max_points: int) -> list[dict[str, float | str]]:
    if len(df) > max_points:
        step = max(1, len(df) // max_points)
        df = df.iloc[::step].copy()
    cols = [
        "time_h",
        "voltage_v_smooth",
        "current_a_smooth",
        "temperature_c_smooth",
        "power_w",
        "voltage_hat_v",
        "voltage_residual_v",
        "baseline_gap_pct",
        "residual_z_score",
        "corrected_voltage_v",
        "corrected_voltage_delta_v",
        "soh_proxy_raw_pct",
        "soh_proxy_pct",
    ]
    if "soh_proxy_raw_pct" not in df.columns:
        df["soh_proxy_raw_pct"] = df["soh_proxy_pct"]
    if "baseline_gap_pct" not in df.columns:
        df["baseline_gap_pct"] = df["voltage_residual_v"] / df["voltage_hat_v"] * 100.0
    if "corrected_voltage_delta_v" not in df.columns:
        df["corrected_voltage_delta_v"] = df["corrected_voltage_v"].diff()
    df["corrected_voltage_delta_v"] = df["corrected_voltage_delta_v"].fillna(0.0)
    out = []
    for row in df[cols].itertuples(index=False):
        point = {col: round(float(value), 6) for col, value in zip(cols, row)}
        point["state"] = _classify_point(
            point["soh_proxy_pct"],
            point["baseline_gap_pct"],
            point["residual_z_score"],
            point["corrected_voltage_delta_v"],
        )
        out.append(point)
    return out


def _classify_point(
    soh_proxy_pct: float,
    baseline_gap_pct: float,
    residual_z_score: float,
    corrected_voltage_delta_v: float,
) -> str:
    _ = residual_z_score
    rapid_drop = corrected_voltage_delta_v < -0.5
    if soh_proxy_pct < 90 or baseline_gap_pct <= -10 or rapid_drop:
        return "Critical"
    if soh_proxy_pct < 95 or baseline_gap_pct <= -5:
        return "Check"
    if soh_proxy_pct < 97 or baseline_gap_pct <= -2:
        return "Warning"
    return "Normal"


def _parse_report(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8")
    wanted = {
        "state": r"Current State:\s*(.+)",
        "soh": r"SOH proxy:\s*(.+)",
        "baseline_gap": r"Baseline voltage gap:\s*(.+)",
        "z": r"(?:Baseline deviation z-score|Voltage residual z-score):\s*(.+)",
        "degradation": r"Degradation speed:\s*(.+)",
        "rul": r"Estimated RUL:\s*(.+)",
    }
    return {
        key: match.group(1).strip()
        for key, pattern in wanted.items()
        if (match := re.search(pattern, text))
    }


def _dataset_label(dataset_id: str) -> str:
    return {
        "normal_data": "Normal operation demo",
        "synthetic": "Synthetic demo",
        "ieee_fc1": "IEEE PHM 2014 FC1",
        "ieee_fc2": "IEEE PHM 2014 FC2",
        "bad_data": "Injected anomaly test",
    }.get(dataset_id, dataset_id)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export compact JSON for the dashboard UI.")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "ui" / "dashboard_data.json")
    parser.add_argument("--max-points", type=int, default=1600)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = export_dashboard_data(args.output, max_points=args.max_points)
    print(f"Dashboard data saved to: {output_path}")


if __name__ == "__main__":
    main()
