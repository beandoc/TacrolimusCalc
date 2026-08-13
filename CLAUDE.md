# TacrolimusCalc — working notes for future sessions

This file exists because the same wrong turn keeps being available: **tuning the
model against a 92-prediction cohort and reading sub-1-point MAPE movements as
progress.** Everything below was measured, not assumed. Re-measure before
trusting it, but do not re-derive it from scratch.

---

## 1. What this app is

One-compartment oral PK (Bateman), apparent CL/F and V/F, **MAP-Bayesian on two
parameters (ηCL, ηV)** with a log-normal prior. KA is population-fixed at 4.53
and never estimated. C/D ratio is **not** the predictive model — it only feeds
the metabolizer label and IPV flags.

- Engine: [pk-engine.js](pk-engine.js) — pure, DOM-free, loaded as a classic script.
- UI: [index.html](index.html) — one 2700-line inline `<script>`, no framework.
- Tests: `npm run verify` (= `check` + `test`). 171 + 13 assertions.
- Backtest: `npm run backtest` → [tools/backtest.js](tools/backtest.js).

Sampling convention: troughs are drawn at **06:45**, 15 min before the 07:00
dose, morning only. Dose rows at 07:00/19:00 are change-points, not
administrations. A 07:15 draw instead of 06:45 is **+30%** on the level.

---

## 2. The evaluation cohort — know these numbers before proposing anything

The eval set is `tacrolimus.db` (gitignored, PHI). After dropping the 5 fixture
MRNs, **10 patients / 122 levels / 92 held-out predictions** at a 3-level warm-up.
Level counts: 23, 20, 15, 14, 10, 10, 9, 8, 7, 6.

Measured on that cohort with `npm run backtest`:

| config | MAPE | 95% CI | bucket acc | missed HIGH |
|---|---|---|---|---|
| stationary | 29.5% | [24.3, 35.4] | 57.6% | 15 |
| adaptive | 25.3% | [20.8, 30.5] | 55.4% | 18 |
| selector (shipped) | 24.6% | [20.4, 29.2] | 58.7% | 17 |
| baseline: carry-last-level | 26.7% | [22.0, 31.9] | 54.3% | 17 |
| baseline: **always predict 8.0** | 31.6% | [25.8, 37.6] | 48.9% | 30 |

**The 95% CI on any absolute MAPE here is ~±5 points, i.e. ~10 points wide.**
A 0.5-point difference is 1/20th of the noise on the measurement.

Paired differences (the only comparison that means anything at this n):

- `stationary − adaptive` = **+4.2 pt, CI [0.1, 8.5]** → adaptive genuinely better.
- `selector − adaptive` = **−0.7 pt, CI [−3.1, 1.6]** → **not distinguishable from zero.**
- `selector − stationary` = **−4.9 pt, CI [−8.5, −1.5]** → genuinely better.
- `adaptive − constant-8` = **−6.3 pt, CI [−12.8, 0.1]** → **not distinguishable from a constant.**

### Where the old 22.0 / 22.5 / 29.7 numbers came from

They were produced by an uncommitted scratchpad `walkforward.js` with
`START_AT = 6` **and `vidyashree` excluded**. Reproduce it exactly:

```
node tools/backtest.js --warmup 6 --drop vidyashree
   stationary 28.6%   adaptive 20.8%   selector 21.6%
```

vidyashree is precisely the patient the selector exists to protect — adaptive
scores 37.3% on her, stationary 29.4%, selector 30.1%. **Dropping her flips the
sign of the headline conclusion.** With the full cohort the selector beats
always-adaptive by 0.7 pt; without her it loses by 0.8 pt. Neither is
significant. Do not quote a cohort-restricted number without saying so.

---

## 3. Why the stationary↔adaptive axis is exhausted

Sweep the weight placed on the most recent level (0 = patient mean, 1 = last
value) on the real cohort: 32.5% → 26.6% → **25.0% at 0.7** → 26.7% at 1.0.
**Flat from 0.4 to 1.0.** Even the oracle-tuned optimum is only −1.7 pt vs
always-adaptive, CI [−4.0, +0.6] — not distinguishable.

Mechanism: **lag-1 autocorrelation of log(C/D) deviation from the patient's own
mean is ρ ≈ +0.44.** More than half of within-patient C/D movement is transient
and reverts. No filter, forgetting factor, or selector can predict it. That is
*why* adaptive ≈ stationary ≈ selector, and it is a property of the biology plus
the assay, not of the code.

**Corollary: further selector / half-life / window tuning is not worth doing.**

---

## 4. The tuning leakage that makes every quoted number optimistic

These were all chosen by looking at the same 92-prediction score:

| knob | value | where |
|---|---|---|
| `RECENCY_REGIMES.adaptive` | 5 d / 0 | [pk-engine.js:109-120](pk-engine.js#L109-L120) |
| `REGIME_SELECT_TRIALS` | 5 | [pk-engine.js:131-135](pk-engine.js#L131-L135) |
| `REGIME_SWITCH_MARGIN` | 0 | [pk-engine.js:834-841](pk-engine.js#L834-L841) |
| `MATURATION_WINDOW_DAYS` | 90 | [pk-engine.js:96-127](pk-engine.js#L96-L127) |

Measured cost: tuning **one** knob on this cohort and reporting its winning
score is **+0.8 points optimistic** vs nested leave-one-patient-out. That single
knob's optimism already exceeds the 0.5-point gap that was being read as a
result. Four knobs were tuned this way.

Worse, the winner is unstable: bootstrapping the 10 patients, the "best"
shrinkage setting is 0.8 in 33% of resampled cohorts, 0.7 in 30%, 0.6 in 19%.
**Selecting a hyperparameter on this cohort is selecting noise.**

Every number the engine quotes is therefore a *resubstitution estimate of a
tuned pipeline*, not held-out performance. `POPULATION_BACKTEST_MAPE = 26` is
rendered to clinicians as a measured figure; it is an optimistic bound.

---

## 5. Leakage audit — verdict (done, keep it that way)

**No within-fit leakage.** Verified:

- `selectRecencyRegime` fits on `sorted.slice(0, i)` — the target never enters
  its own fit ([pk-engine.js:871](pk-engine.js#L871)).
- `predictAtTime` filters `d.time < t` — no future dose reaches a prediction.
- `cl_scalar` is 1.18 (the default) in the DB; no cohort-fitted override.

Guarded permanently by the section *"Walk-forward selection is leak-free"* in
[test/pk-engine.test.js](test/pk-engine.test.js) — a behavioural probe: triple a
held-out level and its own prediction must not move (with a live-probe
counter-assertion so it cannot pass trivially).

**Theoretical look-ahead, measured and found to be zero on this cohort:**
`fillHistoricalGaps` interpolates the dose regimen across gaps using the whole
log, so a dose entered *after* a target could in principle shape the inferred
regimen *before* it. `node tools/backtest.js --doses truncated` cuts the dose
history at each target before fitting; on the current cohort it reproduces the
`full` numbers exactly (24.6% selector MAPE either way) — this cohort's dose
logs have no gaps that straddle a held-out target. Re-run this check if the
cohort changes; it is not guaranteed to stay zero.

---

## 6. Traps — things that look like breakthroughs and are not

**Dose is endogenous.** `d log(dose_next) / d log(trough_prev) = −0.31`: the
clinician cuts the dose when the last trough was high. Consequence: naive
within-patient regression gives `d log(trough)/d log(dose) ≈ +0.4`, not +1.0, in
every patient. This looks like "PK is not dose-proportional!" **It is
regression-to-the-mean from a closed feedback loop.** Acting on the 0.4
elasticity would make the engine recommend ~2.5× larger dose changes than
needed — a safety bug, not an accuracy win. Never fit dose response from
clinician-adjusted observational data without breaking the loop.

**Time-post-transplant clearance curve — tested and rejected.** A pooled
`C/D = 0.60 × days^0.28` fits well per-patient (pooled R² = 0.56) and scored
**53.9% MAPE** in a forward-only backtest. The trend is confounded with the fact
that doses are tapered over the same period.

**Covariates carry no information in this cohort.** 9 of 10 patients have
identical albumin (3.5), bilirubin (1.0), genotype (unknown), inhibitor (none).
Hematocrit and albumin are collected and displayed but deliberately **not** in
the fit. Adding covariate sophistication cannot help until the data varies.

**Maturation submodel was already removed** after Sangha 2024 (n=103) showed it
made things worse: MAPE 40% → 53%. Do not reintroduce it.

---

## 7. Where the real headroom is

1. **Calibrated intervals instead of a point estimate.** Log-residual SD ≈ 0.38;
   empirical coverage of a nominal 90% band is already ~89%. An honest 90% band
   on a 9.0 ng/mL prediction spans **4.8–16.8**. Ship the band and per-dose
   in-range probabilities ("at 4 mg: 62% in range, 25% sub-therapeutic"), and
   score with CRPS / pinball / coverage. **Coverage is the one thing 92 points
   can actually validate.**
2. **Fit on log(concentration).** Residuals are multiplicative; raw-scale
   fitting is what produces the observed ~14% over-prediction bias.
3. **Score the decision, not the number.** Bucket accuracy is 55–59% and
   17–18 of 92 predictions miss a HIGH trough. A MAPE gain that changes no dose
   recommendation has delivered nothing. Missing HIGH (toxicity) and missing LOW
   (rejection) are not equal-cost.
4. **More patients.** n=92 from 10 patients is below the threshold at which any
   modelling decision is verifiable. The data pipeline outranks every algorithm
   change on this list.

---

## 8. How to work on this codebase

- **Never compare two independent MAPEs.** Same predictions, both configs,
  paired difference, bootstrap CI. `tools/backtest.js` only reports it that way
  on purpose.
- **Act only on a paired CI that excludes zero.** Otherwise say "not
  distinguishable" — that is a finding, not a failure.
- **State the cohort and warm-up with every number.** Dropping a patient is a
  claim, not a convenience; `--drop` echoes the exclusion in the header for that
  reason.
- **Check the humility baselines first.** "Always predict 8.0" scores 31.6% and
  carry-last-level scores 26.7%. A change that does not clear those is not doing
  what it claims.
- **Look at the per-patient table before believing a mean.** arun_chougle (n=5)
  and raj_bahadur_test (n=3) swing between 24% and 87% depending on config;
  8 predictions can move the aggregate a point.
- **Golden-value tests must constrain, not track.** Six hand-verified numbers
  were once rewritten to match new output. If a golden number changes, verify
  the new value independently or the test has stopped being a test.
- **Statistical code fails quietly — a wrong number looks like a right number.**
  That is the stated philosophy of [test/README.md](test/README.md); every
  assertion there encodes a mistake that was actually made. Keep it that way.
- **PHI**: `tacrolimus.db`, `patient_*.json`, `extracted_*`, `line_*_content.json`
  are real patient data and gitignored. Never commit them, never paste level
  data into a summary, never publish them.
- The app is clinical decision support, **not** a regulated device. External
  validation (Sangha 2024, n=103) reports ~40% MAPE. When the internal number
  and the external number disagree on screen, the external one is the honest
  guide.
