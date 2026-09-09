# Two-moment pulse repair: local scalar positivity — v0.15

## 공개 검토용 연구 초안 / Unreviewed research draft

**두 모멘트 펄스 보정의 상쇄 보존 추정과 전 구간 각방향 양성**  
*제안된 Navier–Stokes 구성의 국소 스칼라 하위문제*

작성일: 2026-09-09. 이 브랜치는 독립적인 연구 노트의 공개 보관용이며, 저장소의 PEMFC 프로젝트와 별개의 자료다.

이 자료는 AI 도구를 사용해 유도·구현·작성한 **독립 검토 전 초안**이다. 인간 저자·소속은 확정하지 않았다. 저널 투고, 동료심사, Lean 커널 검증 또는 전체 Navier–Stokes 문제의 검증을 완료한 자료가 아니다. 학술적 신규성도 확인되지 않았다.

## 읽을 자료

- [한국어 LaTeX 원고 전체](manuscript/paper_ko.tex): 정의, 정리, 증명, 수치 사례, 한계와 참고문헌.
- [결과 설명](RESULT_SUMMARY_KO.md): 수학적 결과의 의미와 적용 범위.
- [주장 점검표](CLAIM_AUDIT_KO.md): 주장별 근거, 미검증 사항, 이전 해석의 정정.
- [출처 명세](SOURCE_MANIFEST.json): 사용한 원문 커밋과 함수, 관련 문헌.
- [정확한 유리수 검산 코드](code/exact_example.py): 표준 라이브러리만으로 한 조건부 사례의 마지막 수치 부등식을 검사한다.

## 초록

로그 반경 좌표에서 주 펄스와 두 개의 분리된 모멘트 보정으로 이루어진 스칼라 이력 문제를 분석한다. 배경 진폭과 모멘트 계수를 곱으로 유지하고, 모멘트 행렬의 지수 인자를 분리하여 보정의 최종 가중 기여를 상계한다. 명시된 국소 입력 제한 아래에서 후반 보정까지 포함한 전체 펄스 구간의 각방향 이력 Q가 양수라는 충분조건을 제시한다. 원고의 조건부 사례에서는 Q > 0.0009679라는 하계를 유리수 산술로 검산한다.

Q는 수학적 보조량이고 y는 로그 반경이다. 위 숫자는 실제 유체의 안전율·속도·성공 확률이 아니다. 원고의 입력이 원래 구성의 모든 입구·매칭·고차 미분·압력·축방향·응력 조건과 동시에 실현되는지는 미검증이다.

## 재현

이 공개본은 원고와 실행 가능한 소스, 검산 자료를 제공한다. PDF와 그림은 아래 명령으로 생성할 수 있으며, 글꼴 파일은 배포하지 않는다.

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
cd code
python repair_budget.py --out ../results
python exact_example.py
python -m pytest -q
python make_figure.py
```

제공된 검산 패키지의 87개 자동 테스트는 게시 준비 과정에서 다시 실행하여 통과했다. 테스트의 수는 수학적 정리의 개수, 신규성 또는 독립적인 학술 검증을 뜻하지 않는다.

XeLaTeX, kotex와 원고에 지정된 글꼴을 갖춘 환경에서:

```bash
cd ../manuscript
xelatex -interaction=nonstopmode -halt-on-error paper_ko.tex
xelatex -interaction=nonstopmode -halt-on-error paper_ko.tex
```

## 검토 요청

모멘트 지수와 Jacobian의 원문 대응, 조건 누락, 가중 상계의 부호와 지수, 전체 구성에 대한 실제 유용성, 선행연구 대비 신규성에 관한 검토가 필요하다. 일반적인 적분인자법이나 행렬 스케일링을 새 기법으로 주장하지 않는다.

출처의 전체 특이점 주장과 이 노트의 조건부 국소 정리는 별개로 검토해야 한다. 출처 URL과 버전은 원고 및 명세에 기재되어 있다. 이 공개 게시 자체는 학술지 게재 또는 원문 저자들의 승인을 의미하지 않는다.
