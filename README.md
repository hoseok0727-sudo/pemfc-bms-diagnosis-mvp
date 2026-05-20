# PEMFC BMS Diagnosis MVP

본 프로젝트는 PEMFC 스택을 대상으로 BMS 수준의 전압, 전류, 온도, 운전시간 데이터를 이용하여 상태진단을 수행하는 MVP이다. 정상 전압 기준선, 전압 잔차, SOH proxy, 열화율, 간이 RUL을 계산하여 조기 이상감지와 점검 필요성 판단을 지원한다. 본 시스템은 세부 전기화학 고장 원인을 확정하는 정밀 진단기가 아니라, BMS 데이터 기반 유지보수 보조 시스템이다.

This project is a minimum viable product for PEMFC stack health diagnosis using BMS-level operating data. The system uses stack voltage, current, temperature, and operating time to estimate a normal voltage baseline, voltage residual, SOH proxy, degradation rate, and simple RUL. It is intended as an early warning and maintenance-support tool, not as a precise electrochemical failure diagnosis system.

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
| Power calculation | `P = V * I` |
| Normal baseline model | Empirical voltage model using current and temperature |
| Voltage residual | Difference between measured voltage and normal predicted voltage |
| SOH proxy | Corrected-voltage ratio versus initial reference |
| Degradation rate | Linear trend of corrected voltage versus operating time |
| Simple RUL | Remaining hours assuming the current degradation trend continues |
| State classification | Normal / Warning / Check / Critical |
| Report output | Text report with interpretation and recommended action |

Excluded:

| Excluded item | Reason |
| --- | --- |
| EIS analysis | Not continuously available in basic BMS data |
| HFR and Rct estimation | Requires separate diagnostic equipment or experiments |
| Cell-level voltage diagnosis | Requires CVM or cell-level sensing |
| Definite root-cause identification | BMS-level data alone is insufficient |
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

## Equations

Power:

```text
P(t) = V_stack(t) * I(t)
```

Empirical normal voltage baseline:

```text
V_hat_normal(t) = beta0 + beta1*I(t) + beta2*I(t)^2 + beta3*T(t) + beta4*I(t)*T(t)
```

In fixed-load aging data where the initial normal section has too little current and temperature variation, the implementation falls back to a constant initial-voltage baseline. This avoids unstable polynomial coefficients while preserving the MVP goal: detecting voltage decline from a normal BMS operating reference.

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

Degradation rate:

```text
V_corr(t) ~= a + k_deg*t
k_deg = dV_corr / dt
```

Simple RUL:

```text
RUL = (V_EOL - V_corr(t)) / k_deg, only when k_deg < 0
```

본 연구의 RUL은 실제 고장 시점을 정확히 예측하는 값이 아니라, 현재 전압 저하 추세가 유지된다는 가정하에 계산한 간이 잔여 운전시간 지표이다.

## State Logic

| State | Example condition | Meaning |
| --- | --- | --- |
| Normal | `SOH > 95%` and `z_V > -2` | Normal operating trend |
| Warning | `90% < SOH <= 95%` or `-3 < z_V <= -2` | Trend monitoring required |
| Check | `80% < SOH <= 90%` or `-4 < z_V <= -3` | Inspection recommended |
| Critical | `SOH <= 80%`, `z_V <= -4`, or rapid voltage drop | Detailed inspection or operation restriction required |

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

The dashboard replays each diagnosis time series as a live-looking stream. It is not connected to a real BMS; it simulates realtime updates from saved diagnosis outputs for MVP demonstration.

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
- It does not use EIS, HFR, Rct, or cell-level voltage diagnosis.
- It does not identify specific electrochemical failure modes such as catalyst degradation, membrane drying, flooding, or gas starvation.
- The RUL value is a simple trend-based estimate, not a precise lifetime prediction.
- Applying this workflow to PAFC, SOFC, MCFC, or other fuel-cell types requires retraining the baseline model and recalibrating thresholds.
