# BMS 데이터 기반 PEMFC 조기 이상감지 및 잔여수명 추정 MVP

This project is a minimum viable product for PEMFC stack early anomaly detection and maintenance support using BMS-level operating data only. It uses stack voltage, current, temperature, and operating time to estimate a normal voltage baseline, voltage residual, voltage-based SOH proxy, degradation speed, and simple RUL.

본 MVP는 EIS, HFR, Rct, 셀별 전압, 유량, 압력, 습도 데이터를 사용하지 않고, BMS에서 비교적 쉽게 얻을 수 있는 스택 전압, 전류, 온도, 운전시간만 사용한다. 따라서 세부 고장 원인을 확정하지는 않으며, 정상 전압 기준선 대비 전압 잔차, 전압 기반 SOH proxy, 열화속도, 간이 RUL을 이용해 Normal, Warning, Check, Critical 상태를 판단하는 조기 이상감지 및 유지보수 보조 시스템으로 설계한다.

## Live Demo Site

This repository is prepared for GitHub Pages deployment from the `docs/` directory.

Expected GitHub Pages URL after deployment:

```text
https://hoseok0727-sudo.github.io/pemfc-bms-diagnosis-mvp/
```

Public demo entry points:

- `docs/index.html`: landing page
- `docs/dashboard.html`: realtime-style dashboard
- `docs/explanation_report.html`: calculation explanation report
- `docs/pemfc_bms_diagnosis_easy_explanation.pdf`: PDF report
- `docs/assets/*.gif`: normal and fault dashboard demos

## Project Scope

Included:

| Item | Description |
| --- | --- |
| BMS data input | Stack voltage, current, temperature, operating time |
| Sensor quality check | Missing values, invalid ranges, spikes, and smoothing |
| Power calculation | `P = V * I` |
| Normal baseline model | Empirical voltage baseline using current and temperature |
| Voltage residual | Difference between measured voltage and predicted normal voltage |
| SOH proxy | Corrected-voltage ratio versus initial reference |
| Degradation speed | Positive voltage decrease speed, `s_deg = -k_deg` |
| Simple RUL | Remaining hours assuming the current voltage decrease speed continues |
| State classification | Normal / Warning / Check / Critical |
| Rapid-drop rule | Escalates status when corrected voltage drops sharply |

Excluded:

| Excluded item | Reason |
| --- | --- |
| EIS analysis | Not continuously available in basic BMS data |
| HFR and Rct estimation | Requires separate diagnostic equipment or experiments |
| Cell-level voltage diagnosis | Requires CVM or cell-level sensing |
| Definite root-cause identification | BMS-level data alone is insufficient |
| Flooding/drying or hydrogen-shortage classification | Requires additional physical signals |
| Generalization to all fuel-cell types | This MVP is scoped to PEMFC |

## Input Data

The standard CSV schema is:

| Column | Unit | Description |
| --- | --- | --- |
| `time_h` | h | Cumulative operating time |
| `voltage_v` | V | Stack voltage |
| `current_a` | A | Stack/load current |
| `temperature_c` | deg C | Stack or coolant temperature |

The code can also read the IEEE PHM 2014 PEMFC aging Excel zip (`FC1_FC2_Excel.zip`) and uses only `Ageing` workbooks. EIS and polarization files are intentionally ignored to keep the MVP within BMS-level data.

## Preprocessing

BMS raw data can include missing values, sensor noise, and spikes. The MVP:

- removes missing values in required columns,
- excludes clearly invalid ranges such as negative voltage or negative current,
- removes impossible temperature values outside `-40` to `120 deg C`,
- applies a rolling median and rolling average to reduce spikes and noise.

## Normal Training Section

The normal baseline and residual statistics are learned only from an initial normal section:

```text
t in [0, t0]
```

Default implementation:

- first `10%` of the dataset by default,
- optionally a larger initial hour window when `--normal-hours` is supplied.

From that normal section, the residual mean and standard deviation are computed:

```text
mu_r = mean(r_V)
sigma_r = std(r_V)
```

This split is important because the later diagnosis section must be compared against a fixed normal reference.

## Equations

Input vector:

```text
X(t) = [V_stack(t), I(t), T(t), t]
```

Power:

```text
P(t) = V_stack(t) * I(t)
```

Empirical normal voltage baseline:

```text
V_hat_normal(t) = beta0 + beta1*I(t) + beta2*I(t)^2 + beta3*T(t) + beta4*I(t)*T(t)
```

전류가 커지면 전압은 자연스럽게 낮아지고, 온도도 전압에 영향을 준다. 따라서 단순히 전압이 낮다고 이상으로 판단하지 않고, 전류와 온도를 반영한 정상 전압 기준선을 먼저 만든다.

This baseline is not a physical electrochemical model. It is an empirical baseline model learned from the initial normal operating data.

In fixed-load aging data where the initial normal section has too little current and temperature variation, the implementation falls back to a constant initial-voltage baseline. This avoids unstable polynomial coefficients while preserving the MVP goal.

Voltage residual:

```text
r_V(t) = V_stack(t) - V_hat_normal(t)
```

Standardized residual:

```text
z_V(t) = (r_V(t) - mu_r) / sigma_r
```

Operating-condition corrected voltage:

```text
V_corr(t) = V_stack(t) - [V_hat_normal(I(t), T(t)) - V_hat_normal(I_ref, T_ref)]
```

SOH proxy:

```text
SOH_proxy(t) = V_corr(t) / V_corr(t0) * 100
```

For dashboard and report display, `SOH_proxy` is capped at `100%`. The uncapped
raw value is kept as `soh_proxy_raw_pct` in the time-series CSV for audit. This
prevents early noise or minor operating-condition correction error from making
the health score look better than new.

Corrected-voltage slope and degradation speed:

```text
V_corr(t) ~= a + k_deg*t
k_deg = dV_corr / dt
s_deg = -k_deg
```

The dashboard and report use `s_deg` as a positive degradation speed. A larger `s_deg` means the corrected voltage is decreasing faster.

EOL threshold:

```text
V_EOL = 0.9 * V_corr(t0)
SOH_proxy,EOL = 90%
```

Simple RUL:

```text
RUL(t) = (V_corr(t) - V_EOL) / s_deg, only when s_deg > 0
```

RUL is not a precise prediction of the actual failure time. It is a maintenance-support indicator calculated under the assumption that the current voltage decrease speed continues.

Rapid corrected-voltage drop:

```text
Delta V_corr(t) = V_corr(t) - V_corr(t - Delta t)
```

If:

```text
Delta V_corr(t) < -Delta V_crit
```

the status is escalated because SOH can still look acceptable even when voltage suddenly drops.

The residual z-score is kept as a reference indicator, not as the main state
decision rule. It is useful for seeing baseline deviation, but it can be too
sensitive when the empirical baseline is imperfect.

The dashboard uses a simpler baseline-gap rule:

```text
baseline_gap_pct = (V_stack - V_hat_normal) / V_hat_normal * 100
```

Degradation speed and RUL are included as trend-support indicators. They do not
override the two main anomaly rules. `s_deg` is displayed as a positive voltage
decrease speed, and RUL is calculated only when recent corrected voltage is
meaningfully decreasing:

```text
RUL(t) = (V_corr(t) - V_EOL) / s_deg
```

## State Logic

| State | Example condition | Meaning | Response |
| --- | --- | --- | --- |
| Normal | `SOH_proxy >= 97%` and `baseline_gap_pct > -2%` | Similar to initial corrected-voltage reference and baseline | Continue monitoring |
| Warning | `95% <= SOH_proxy < 97%` or `baseline_gap_pct <= -2%` | Voltage decline signal starts | Strengthen trend monitoring |
| Check | `90% <= SOH_proxy < 95%` or `baseline_gap_pct <= -5%` | Meaningful performance decline | Check cooling/load/operating conditions |
| Critical | `SOH_proxy < 90%`, `baseline_gap_pct <= -10%`, or rapid corrected-voltage drop | Operation restriction may be needed | Detailed inspection and operation review |

Thresholds are temporary MVP defaults. Real deployment should recalibrate them using stack manufacturer criteria and operating data.

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Run With Synthetic Sample Data

```bash
python -m src.main
```

If `data/sample_pemfc_data.csv` does not exist, it is generated automatically.

Outputs:

- `outputs/diagnosis_report.txt`
- `outputs/diagnosis_timeseries.csv`
- `outputs/voltage_vs_time.png`
- `outputs/corrected_voltage_vs_time.png`
- `outputs/soh_proxy_vs_time.png`
- `outputs/residual_z_score_vs_time.png`

## Run With IEEE PHM 2014 Excel Zip

Example:

```bash
python -m src.main --ieee-zip "%USERPROFILE%\Downloads\FC1_FC2_Excel.zip" --fuel-cell FC1
python -m src.main --ieee-zip "%USERPROFILE%\Downloads\FC1_FC2_Excel.zip" --fuel-cell FC2
```

The loader searches for `Ageing` workbooks and maps common column names to:

- `time_h`
- `voltage_v`
- `current_a`
- `temperature_c`

If a workbook uses unexpected headers, update `COLUMN_ALIASES` in `src/preprocess.py`.

## Realtime-style Dashboard UI

After running the pipelines, export compact dashboard data and start a local static server:

```bash
python -m src.export_dashboard_data
python -m http.server 8765 -d ui
```

Open:

```text
http://localhost:8765/dashboard.html
```

The dashboard replays each time series as a live-looking stream. It is not connected to a real BMS; it simulates realtime updates from saved diagnosis outputs for MVP demonstration.

The calculation explanation report is also available at:

```text
http://localhost:8765/explanation_report.html
```

PDF version:

```text
docs/pemfc_bms_diagnosis_easy_explanation.pdf
```

Dashboard GIF demos:

```text
docs/assets/normal_operation_dashboard.gif
docs/assets/fault_critical_dashboard.gif
```

## GitHub Pages Deployment

1. Create a public GitHub repository named `pemfc-bms-diagnosis-mvp`.
2. Push this project to the repository.
3. In GitHub, open `Settings` -> `Pages`.
4. Set `Source` to `Deploy from a branch`.
5. Select branch `main` and folder `/docs`.
6. Save. The site will be published at:

```text
https://hoseok0727-sudo.github.io/pemfc-bms-diagnosis-mvp/
```

## Limitations

- This MVP is limited to PEMFC.
- It uses only BMS-level stack voltage, current, temperature, and operating time.
- It does not use EIS, HFR, Rct, cell-level voltage, flow, pressure, or humidity data.
- It does not identify specific electrochemical failure modes such as catalyst degradation, membrane drying, flooding, or gas starvation.
- The RUL value is a simple trend-based estimate, not a precise lifetime prediction.
- Applying this workflow to PAFC, SOFC, MCFC, or other fuel-cell types requires retraining the baseline model and recalibrating thresholds.
