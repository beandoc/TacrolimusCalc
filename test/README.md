# Tests

```sh
npm install     # once
npm run setup   # once — enables the pre-commit hook
npm run verify  # syntax + wiring checks, then all tests
```

| Command | What it does |
|---|---|
| `npm run check` | Parses `pk-engine.js` and the inline `<script>` in `index.html`; verifies the engine stays DOM-free and nothing is defined twice. |
| `npm test` | Unit tests for the PK engine, plus an integration smoke test. |
| `npm run verify` | Both of the above. What the pre-commit hook runs. |

## Files

- **`check-html.js`** — the cheap, high-value one. Everything in `index.html`
  lives in a single `<script>` block, so one missing brace anywhere kills all 57
  functions and the page silently does nothing. That has happened. ~50 ms to
  rule out.
- **`pk-engine.test.js`** — 60 assertions over the numerical core. Includes an
  independent reimplementation of the Laplace posterior, so the engine is
  checked against maths rather than against its own past output.
- **`smoke.test.js`** — drives the real `runForecast()` from `index.html`
  behind DOM stubs. Deliberately thin: it only catches wiring failures the unit
  tests structurally cannot (signature drift, NaNs reaching the chart). It
  asserts on values, never on markup, so restyling the UI will not break it.

## Why these assertions exist

Statistical code fails quietly — a wrong number looks like a right number. Each
of these encodes a mistake that was actually made or actually proposed:

- **A posterior can never be wider than its prior.** `Σ_post = (Ω⁻¹ + JᵀJ/σ²)⁻¹`
  with `JᵀJ` positive semi-definite gives `Σ_post ⪯ Ω`, so `sd(ηV) ≤ √0.090 =
  0.300`. A proposed rewrite reported 0.52. It looks plausible in a diff; it is
  impossible.
- **CL and V are positively correlated in the trough-identified posterior.**
  `∂lnC/∂lnCL = −1.20` and `∂lnC/∂lnV = +0.20` have opposite signs, so raising V
  lowers `ke = CL/V` and raises the trough. An earlier revision hardcoded
  `ρ = −0.45`, the wrong sign.
- **Troughs barely identify V.** `sd(ηV)` stays near its prior no matter how many
  troughs accumulate. Any change that shrinks it is measuring something that
  isn't there.
- **3 SD on ηCL is 1.154**, not 1.26 — `3·√0.148`. Derive thresholds from ω;
  don't hardcode them.
- **`null < 0.5` is `true` in JavaScript.** `r2` is `null` at n=1, which once put
  a red "DO NOT TRUST FORECAST" banner on every single-trough patient.
- **A lagging dose log is not a drug washout.** Without bridge doses the
  optimizer modelled a patient who had stopped taking tacrolimus and recommended
  roughly double the correct dose (−58% predicted trough at a 3-day lag).
- **C0 is drawn 15 min pre-dose.** A sample timestamped 15 min *after* the dose
  reads ~30% high and drives the MAP fit toward slower clearance.

## Reference values

Exact Laplace posterior, 3 mg BID steady state, C0 at 06:45.
Prior: `sd(ηCL) = 0.385`, `sd(ηV) = 0.300`.

| observations | sd(ηCL) | sd(ηV) | ρ |
|---|---|---|---|
| 1 trough | 0.147 | 0.298 | +0.285 |
| 2 troughs | 0.112 | 0.298 | +0.399 |
| 3 troughs | 0.097 | 0.298 | +0.475 |
| 5 troughs | 0.081 | 0.298 | +0.576 |
| 4 mixed (C0 + C2) | 0.085 | 0.267 | +0.166 |

## Note on `pk-engine.js`

Loaded as a **classic** script, not `type="module"`. ES modules do not load over
`file://`, and this tool is designed to run straight off disk with no server
(see the `protocol === "file:"` branches in `index.html`). Keep it that way.
