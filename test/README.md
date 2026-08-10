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
- **C0 is drawn 15 min pre-dose, in the MORNING only.** A sample timestamped 15
  min *after* the dose reads ~30% high and drives the MAP fit toward slower
  clearance. Site protocol assays the 06:45 trough and never the evening one, so
  `nextC0Time` must never return 18:45 — offering it as a default invited
  timestamps for draws that do not happen.
- **Rebuild a params object by SPREADING, never by listing `{CL, V, KA}`.**
  `predictAtTime` applies the time-varying weight covariate only when
  `params.weight` is present. Four functions rebuild params — the MAP objective,
  `indParams`, the Laplace posterior and the Monte Carlo sampler — and every one
  that enumerated fields dropped the covariate silently. The fit then modelled a
  constant-weight patient while the plotted population curve modelled a varying
  one: MAP recovered `ηCL = −0.151` (CL 14% low) on data the model generated
  itself, and the individual curve jumped 30% the moment the first level was
  entered. Assert the round-trip through `mapBayesian`, not just `predictAtTime`
  — the original test only exercised the latter, which is why this survived.
- **The target range is in REPORTED units.** `suggestStartingDose` must apply the
  same `cmiaAdjust` the fit and the optimizer do. Without it a CMIA centre aiming
  at 10–11 was told 6 mg/day, which that lab reports as 11.98.
- **An observation the fit cannot use must not be scored.** The OFV's
  `predAdj <= 0.1` branch is a flat +1000, constant in the etas, so a level drawn
  before any logged dose cannot move the estimate — counting it as a 100% error
  put RMSE 5.0 and a "Severe Model Mismatch" banner on an otherwise exact fit.
- **0 is a value, not an absence.** `parseInt(x) || fallback` turned "No MPA"
  (`value="0"`) into MPA = yes on every save, applying a spurious 13% CL
  reduction; `r.dose || ''` opened a 0 mg held dose with a blank field and then
  deleted the record. Use `??` and `Number.isFinite`.

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

## Backend environment variables

The API holds patient-identifiable data. Three variables govern whether it is
safe to expose; the app starts without them and warns loudly, so that local
development keeps working — but a deployment missing them is not protected.

| Variable | Effect if unset |
|---|---|
| `TACRO_API_TOKEN` | **No authentication.** Every endpoint serves and overwrites patient data for anyone who can reach the host. Set it, and enter the same value under ⚙ Settings → API Access Token in the browser. |
| `TACRO_ALLOWED_ORIGINS` | CORS falls back to `*`. Set to a comma-separated list of the origins that serve the app. |
| `DATABASE_URL` | Falls back to local SQLite. On Vercel/Lambda the app now refuses to start rather than write to an ephemeral filesystem and lose records silently. |

`POST /api/patients/{mrn}/events` replaces the whole history, so it rejects any
write that would *reduce* the stored event count unless called with
`?allow_shrink=true`. Only a genuine deletion (or applying a starting protocol,
which clears doses by design) sends that flag. This is what stops a stale local
cache — loaded while the backend was briefly unreachable — from overwriting a
fuller server-side record on the next autosave.

## Note on `pk-engine.js`

Loaded as a **classic** script, not `type="module"`. ES modules do not load over
`file://`, and this tool is designed to run straight off disk with no server
(see the `protocol === "file:"` branches in `index.html`). Keep it that way.
