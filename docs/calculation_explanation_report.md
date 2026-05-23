# PEMFC BMS 조기 이상감지 및 잔여수명 추정 MVP 계산 설명 보고서

## 1. 아주 쉽게 말하면

이 시스템은 "현재 전류와 온도라면 정상 PEMFC는 이 정도 전압이 나와야 한다"는 기준선을 만들고, 실제 전압이 그 기준보다 얼마나 낮아졌는지 보고 상태를 판단한다.

정확한 이름은 **BMS 데이터 기반 PEMFC 조기 이상감지 및 잔여수명 추정 MVP**이다. 세부 고장 원인을 확정하는 시스템이 아니다.

사용 데이터는 네 가지이다.

| 데이터 | 의미 |
| --- | --- |
| 전압 | 연료전지가 실제로 내는 스택 전압 |
| 전류 | 부하가 얼마나 걸렸는지 |
| 온도 | 스택 또는 냉각수 조건 |
| 운전시간 | 시간이 지나며 전압이 떨어지는지 보는 기준 |

## 2. 계산 흐름

1. 전압, 전류, 온도, 운전시간 데이터를 읽는다.
2. 결측값, 비정상 범위, 스파이크를 처리하고 이동평균/median filter로 노이즈를 줄인다.
3. 전체 데이터의 처음 10% 또는 초기 24시간을 정상 기준 구간으로 설정한다.
4. 초기 정상 구간에서 정상 전압 기준선과 잔차 통계 `mu_r`, `sigma_r`를 만든다.
5. 실제 전압과 기준선의 차이, SOH proxy, 열화속도, 급격한 전압 하락 여부를 계산한다.
6. Normal, Warning, Check, Critical 중 하나로 표시한다.

## 3. 정상 전압 기준선

전류가 커지면 전압은 자연스럽게 낮아지고, 온도도 전압에 영향을 준다. 따라서 단순히 전압이 낮다고 이상으로 판단하지 않고, 전류와 온도를 반영한 정상 전압 기준선을 먼저 만든다.

```latex
\hat{V}_{\mathrm{normal}}(t)=
\beta_0+\beta_1 I(t)+\beta_2 I(t)^2+\beta_3 T(t)+\beta_4 I(t)T(t)
```

이 기준선은 물리 모델이 아니라 초기 정상 운전 데이터로 학습한 경험식이다. 그래서 발표에서는 "BMS 데이터 기반 baseline model"이라고 설명하면 된다.

## 4. 전압 잔차

실제 전압에서 정상 기준 전압을 뺀 값이다.

```latex
r_V(t)=V_{\mathrm{stack}}(t)-\hat{V}_{\mathrm{normal}}(t)
```

0 근처면 정상 기준과 비슷하다. 음수로 커질수록 같은 전류와 온도 조건에서도 정상보다 전압이 낮다는 의미다.

## 5. 잔차 z-score

전압 잔차가 정상 구간에서 보이던 흔한 변동인지, 보기 드문 큰 하락인지 나타낸다.

```latex
z_V(t)=\frac{r_V(t)-\mu_r}{\sigma_r}
```

여기서 `mu_r`와 `sigma_r`는 초기 정상 기준 구간에서 계산한다.

```latex
\mu_r=\mathrm{mean}(r_V), \qquad \sigma_r=\mathrm{std}(r_V)
```

`z = -2` 정도면 Warning, `z = -3` 정도면 Check, `z = -4` 이하면 Critical 신호로 볼 수 있다.

## 6. 운전조건 보정 전압

전류와 온도 조건 차이 때문에 생긴 전압 변화를 제거한 값이다.

```latex
V_{\mathrm{corr}}(t)=V_{\mathrm{stack}}(t)-
\left[
\hat{V}_{\mathrm{normal}}(I(t),T(t))-
\hat{V}_{\mathrm{normal}}(I_{\mathrm{ref}},T_{\mathrm{ref}})
\right]
```

서로 다른 운전 조건의 전압을 같은 기준 조건처럼 비교하기 위해 사용한다.

## 7. SOH proxy

초기 정상 상태의 보정 전압을 100%로 놓고, 현재 보정 전압이 몇 % 수준인지 계산한다.

```latex
SOH_{\mathrm{proxy}}(t)=
\frac{V_{\mathrm{corr}}(t)}{V_{\mathrm{corr}}(t_0)}\times100
```

이 값은 실제 완전한 SOH가 아니라 전압 기반 SOH proxy이다.

## 8. 열화속도

최근 구간에서 보정 전압이 시간에 따라 얼마나 빨리 떨어지는지 보는 값이다.

```latex
V_{\mathrm{corr}}(t)\approx a+k_{\mathrm{deg}}t
```

전압이 감소하면 `k_deg`는 음수가 된다. 발표와 대시보드에서는 더 이해하기 쉽게 양수 열화속도 `s_deg`를 사용한다.

```latex
s_{\mathrm{deg}}=-k_{\mathrm{deg}}
```

`s_deg`가 클수록 보정 전압이 빠르게 낮아진다는 뜻이다.

## 9. EOL 기준과 간이 RUL

BMS-only MVP에서는 EOL 기준을 전압 기반으로 둔다.

```latex
V_{\mathrm{EOL}}=0.9V_{\mathrm{corr}}(t_0)
```

같은 말로 쓰면:

```latex
SOH_{\mathrm{proxy,EOL}}=90\%
```

현재 열화속도가 그대로 유지된다고 가정했을 때, EOL 기준 전압까지 남은 시간을 계산한다.

```latex
RUL(t)=\frac{V_{\mathrm{corr}}(t)-V_{\mathrm{EOL}}}{s_{\mathrm{deg}}}
\qquad (s_{\mathrm{deg}}>0)
```

RUL은 실제 고장 시점을 정확히 예측하는 값이 아니라, 현재 전압 저하 추세가 유지된다는 가정하의 유지보수 참고 지표이다.

## 10. 급격한 전압 하락 rule

장기 SOH가 아직 괜찮아 보여도 전압이 갑자기 크게 떨어지면 조기 이상으로 봐야 한다.

```latex
\Delta V_{\mathrm{corr}}(t)=
V_{\mathrm{corr}}(t)-V_{\mathrm{corr}}(t-\Delta t)
```

```latex
\Delta V_{\mathrm{corr}}(t)<-\Delta V_{\mathrm{crit}}
```

이 조건이 만족되면 Check 또는 Critical로 상태를 올린다.

## 11. 상태 판정

| 상태 | 기준 예시 | 의미 | 대응 |
| --- | --- | --- | --- |
| Normal | `SOH >= 97%`, `z > -2` | 정상 기준과 유사 | 일반 모니터링 |
| Warning | `95% <= SOH < 97%` 또는 `z <= -2` | 전압 저하 신호 시작 | 추세 확인 강화 |
| Check | `90% <= SOH < 95%` 또는 `z <= -3` | 의미 있는 성능 저하 | 냉각/부하/운전조건 점검 |
| Critical | `SOH < 90%` 또는 `z <= -4` 또는 급격한 전압 하락 | 운전 제한 검토 필요 | 정밀 점검/운전 조건 재검토 |

## 12. 한계

이 시스템은 EIS, HFR, Rct, 셀별 전압, 유량, 압력, 습도 데이터를 사용하지 않는다. 따라서 촉매 열화, 막 건조, flooding, 수소 부족 같은 구체적인 고장 원인을 확정하지 않는다.

정확한 표현은 **BMS 데이터 기반 조기 이상감지 및 유지보수 보조 시스템**이다.
