from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from .baseline_model import BaselineVoltageModel
from .diagnosis import calculate_diagnosis_features, summarize_diagnosis
from .preprocess import (
    generate_synthetic_sample,
    load_csv,
    load_ieee_phm2014_zip,
    preprocess_data,
)
from .report import generate_report, save_report


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SAMPLE_PATH = PROJECT_ROOT / "data" / "sample_pemfc_data.csv"
DEFAULT_REPORT_PATH = PROJECT_ROOT / "outputs" / "diagnosis_report.txt"


def run_pipeline(
    data_path: str | Path | None = None,
    ieee_zip_path: str | Path | None = None,
    fuel_cell: str = "FC1",
    normal_fraction: float = 0.2,
    output_dir: str | Path = PROJECT_ROOT / "outputs",
) -> tuple[Path, list[Path]]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if ieee_zip_path is not None:
        raw = load_ieee_phm2014_zip(ieee_zip_path, fuel_cell=fuel_cell)
    else:
        data_path = Path(data_path or DEFAULT_SAMPLE_PATH)
        if not data_path.exists():
            generate_synthetic_sample(data_path)
        raw = load_csv(data_path)

    processed = preprocess_data(raw)
    model = BaselineVoltageModel().fit(processed, normal_fraction=normal_fraction)
    diagnosed = calculate_diagnosis_features(processed, model)
    summary = summarize_diagnosis(diagnosed)

    report = generate_report(summary)
    report_path = save_report(report, output_dir / "diagnosis_report.txt")
    plot_paths = create_plots(diagnosed, output_dir)
    diagnosed.to_csv(output_dir / "diagnosis_timeseries.csv", index=False)
    return report_path, plot_paths


def create_plots(df, output_dir: str | Path) -> list[Path]:
    output_dir = Path(output_dir)
    plot_specs = [
        ("voltage_vs_time.png", "voltage_v_smooth", "Stack voltage (V)", "Voltage vs Time"),
        (
            "corrected_voltage_vs_time.png",
            "corrected_voltage_v",
            "Corrected voltage (V)",
            "Corrected Voltage vs Time",
        ),
        ("soh_proxy_vs_time.png", "soh_proxy_pct", "SOH proxy (%)", "SOH Proxy vs Time"),
        (
            "residual_z_score_vs_time.png",
            "residual_z_score",
            "Residual z-score",
            "Voltage Residual Z-score vs Time",
        ),
    ]

    paths: list[Path] = []
    for filename, y_col, ylabel, title in plot_specs:
        fig, ax = plt.subplots(figsize=(9, 4.8))
        ax.plot(df["time_h"], df[y_col], linewidth=1.4)
        ax.set_title(title)
        ax.set_xlabel("Operating time (h)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = output_dir / filename
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PEMFC BMS diagnosis MVP.")
    parser.add_argument("--data", type=Path, default=DEFAULT_SAMPLE_PATH, help="Standard CSV input.")
    parser.add_argument("--ieee-zip", type=Path, default=None, help="IEEE PHM 2014 FC1/FC2 Excel zip.")
    parser.add_argument("--fuel-cell", choices=["FC1", "FC2"], default="FC1")
    parser.add_argument("--normal-fraction", type=float, default=0.2)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report_path, plot_paths = run_pipeline(
        data_path=args.data,
        ieee_zip_path=args.ieee_zip,
        fuel_cell=args.fuel_cell,
        normal_fraction=args.normal_fraction,
        output_dir=args.output_dir,
    )
    print(f"Report saved to: {report_path}")
    for path in plot_paths:
        print(f"Plot saved to: {path}")


if __name__ == "__main__":
    main()
