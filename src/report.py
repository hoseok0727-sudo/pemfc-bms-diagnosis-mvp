from __future__ import annotations

from pathlib import Path

from .diagnosis import DiagnosisSummary


def generate_report(summary: DiagnosisSummary) -> str:
    rul_text = (
        f"{summary.estimated_rul_h:.1f} h"
        if summary.estimated_rul_h is not None
        else "N/A (degradation trend is not negative)"
    )
    interpretation = _interpretation(summary)
    recommended_action = _recommended_action(summary.state)

    return f"""PEMFC Diagnosis Report

Current State: {summary.state}

SOH proxy: {summary.soh_proxy_pct:.2f} %
Voltage residual z-score: {summary.residual_z_score:.2f}
Degradation rate: {summary.degradation_rate_v_per_h:.6f} V/h
Latest corrected voltage: {summary.latest_corrected_voltage_v:.3f} V
EOL corrected voltage threshold: {summary.eol_voltage_v:.3f} V
Estimated RUL: {rul_text}

Interpretation:
{interpretation}

Recommended Action:
{recommended_action}

Limitations:
- This MVP uses only BMS-level stack voltage, current, temperature, and operating time.
- It does not use EIS, HFR, charge transfer resistance, or cell-level voltage diagnosis.
- It cannot identify a specific electrochemical failure mode.
- RUL is a simple trend-based operating-time indicator, not a precise lifetime prediction.
"""


def save_report(report: str, output_path: str | Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report, encoding="utf-8")
    return output_path


def _interpretation(summary: DiagnosisSummary) -> str:
    lines = []
    if summary.residual_z_score <= -2:
        lines.append("- Stack voltage is lower than the learned normal baseline.")
    else:
        lines.append("- Stack voltage is close to the learned normal baseline.")
    if summary.soh_proxy_pct < 95:
        lines.append("- Corrected voltage shows measurable degradation versus the initial reference.")
    else:
        lines.append("- Corrected voltage remains near the initial reference level.")
    lines.append("- The result is an early warning signal based on BMS data, not a root-cause diagnosis.")
    return "\n".join(lines)


def _recommended_action(state: str) -> str:
    if state == "Normal":
        return "- Continue normal monitoring."
    if state == "Warning":
        return "\n".join(
            [
                "- Review cooling and temperature control stability.",
                "- Check hydrogen and air supply trends during similar load conditions.",
                "- Continue trend monitoring with a shorter review interval.",
            ]
        )
    if state == "Check":
        return "\n".join(
            [
                "- Inspect cooling and humidification-related operating conditions.",
                "- Verify hydrogen and air supply stability.",
                "- Schedule detailed inspection if voltage decline persists.",
            ]
        )
    return "\n".join(
        [
            "- Restrict operation if the trend continues.",
            "- Perform detailed inspection with additional diagnostic equipment.",
            "- Review stack safety limits and manufacturer maintenance criteria.",
        ]
    )
