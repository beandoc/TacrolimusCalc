/**
 * pk-engine.js — Tacrolimus population-PK / MAP-Bayesian numerical core.
 *
 * Pure computation only: no DOM, no fetch, no app state. Everything here is a
 * function of its arguments (the single exception is PK_MODEL.INDIAN_CL_SCALAR,
 * which reads a centre-calibration override from localStorage when one exists).
 *
 * Extracted from index.html so the safety-critical maths can be unit-tested
 * directly — see test/. Before the split, testing it meant parsing the HTML and
 * running the whole script in a sandbox behind DOM stubs, which broke whenever
 * unrelated UI code changed.
 *
 * Loaded as a CLASSIC script (not type="module") on purpose: classic scripts
 * work over file://, ES modules do not, and this tool is designed to run
 * straight off disk with no server (see the protocol === "file:" branches in
 * index.html). In the browser it assigns its exports onto window, so calling
 * code needs no import. Under Node it is an ordinary CommonJS module.
 *
 * Requires dayjs (global in the browser, require()d under Node).
 */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) {
        module.exports = factory(require('dayjs'));      // Node / tests
    } else {
        Object.assign(root, factory(root.dayjs));        // browser globals
    }
}(typeof self !== 'undefined' ? self : this, function (dayjs) {
    'use strict';

    // ============================================================
    // VALIDATED POPULATION PK MODEL
    // Based on: Størset et al. (Transplantation 2014) — adult kidney
    // Adapted for Indian population (higher CYP3A5*1 allele frequency ~45%)
    //
    // All parameters are APPARENT: CL/F and V/F
    // Bioavailability F (~14%) is embedded — do NOT apply F separately
    //
    // Covariates applied:
    //   CYP3A5: effect on CL/F (validated)
    //   Weight: allometric scaling (0.75 on CL, 1.0 on V)
    //   MPA:    reduces CL/F by ~13% (drug interaction)
    //   HCT:    kept separate (whole-blood assay correction only,
    //           not a covariate on CL in whole-blood-based PK)
    // ============================================================
    const PK_MODEL = {
        TVCL: 19.2,         // L/hr  — apparent CL/F, ref: *3/*3, 60kg, no MPA (Størset 2014, Caucasian)
        TVV: 694,          // L     — apparent V/F
        TVKA: 4.53,         // hr⁻¹  — first-order absorption

        // Indian population calibration scalar (Sangha 2024, n=103 Indian renal Tx patients).
        // Caucasian Størset priors overpredict troughs by ~24% (pred > obs) in Indian patients.
        // ~12% attributable to lower oral bioavailability (higher apparent CL/F) in Indian patients,
        // ~12% to Abbott CMIA immunoassay inflation vs LC-MS/MS reference.
        // Net correction: increase TVCL by 18% so MAP Bayesian starts from a less-biased prior.
        // Override: saved by Center Calibration Report via localStorage key 'tacro_cl_scalar'.
        // Guarded so the engine also loads under Node (tests), where there is no
        // localStorage — tests then always see the uncalibrated 1.18 default.
        INDIAN_CL_SCALAR: (typeof localStorage !== 'undefined'
            && parseFloat(localStorage.getItem('tacro_cl_scalar'))) || 1.18,

        // Inter-individual variability (variance of ln-normal distribution)
        // ωCL = 40% CV  → var = ln(1 + 0.40²) = 0.148
        // ωV  = 31% CV  → var = ln(1 + 0.31²) = 0.090
        OMEGA_CL: 0.148,
        OMEGA_V: 0.090,

        // Proportional residual error (σ = 18% CV)
        SIGMA: 0.18,

        // Outlier alert threshold, in prior SD units. Deliberately DECOUPLED from
        // the MAP search bounds (ETA_*_BOUNDS): the search stays wide so genuine
        // outliers (severe hepatic impairment, strong CYP3A induction) converge
        // without artificial truncation, while the clinical banner fires at 3 SD.
        //   3 SD on ηCL = 3·√0.148 = 1.154
        //   3 SD on ηV  = 3·√0.090 = 0.900
        OUTLIER_SD: 3,

        // CYP3A5 covariate on CL/F
        // NOTE: Indian population has ~45% *1 allele frequency
        //       'unknown' defaults to intermediate-high (conservative for safety)
        CYP3A5: { '11': 1.79, '13': 1.35, '33': 1.00, 'unknown': 1.30 },

        // Weight allometric exponents
        WT_CL: 0.75,
        WT_V: 1.00,

        // MPA co-administration: inhibits CYP3A → reduces apparent CL by ~13%
        MPA_CL: 0.87,

        // v5.2 Advanced Clinical Weights
        BIL_REF: 1.0,
        INHIBITOR: {
            'none': 1.0,
            'moderate': 0.65, // ~35% reduction (e.g. Diltiazem)
            'strong': 0.25    // ~75% reduction (e.g. Voriconazole)
        },

        // HCT RBC binding (for quick correction tool only)
        HCT_BMAX: 418,
        HCT_KD: 3.8,

        // Reference values
        WT_REF: 60

        // NOTE: HCT_CL, ALB_REF and HCT_REF were removed here — they were dead
        // (referenced nowhere) and misleadingly implied HCT/albumin are
        // structural CL covariates in this whole-blood model. They are not:
        // HCT variability is already inside ωCL/ωV (see below), and albumin
        // has no fitted exponent in the source model.
    };

    // ============================================================
    // DOSING & SAMPLING CONVENTION
    //
    // Standard BD schedule: 07:00 and 19:00.
    // C0 (pre-dose trough) is drawn SAMPLING_LEAD_MIN minutes BEFORE the dose,
    // i.e. 06:45 and 18:45.
    //
    // Why this matters numerically: predictAtTime() includes a dose only when
    // dose.time < t (strict), so evaluating at exactly 07:00 happens to give a
    // pre-dose value — but it sits precisely on a discontinuity. At 3mg BID
    // steady state:
    //     06:45 (true C0)          →  9.15 ng/mL
    //     07:00 (on the boundary)  →  9.08 ng/mL   (−0.8%, harmless)
    //     07:01 (one minute late)  →  9.39 ng/mL   (+2.6%)
    //     07:15 (15 min late)      → 11.92 ng/mL   (+30.2%)
    // Defaulting every timestamp to 06:45/18:45 both matches the protocol and
    // keeps every evaluation a safe 15 minutes clear of that cliff.
    // ============================================================
    const DOSE_SLOTS = { AM: 7, PM: 19 };   // hour-of-day for BD dosing
    const SAMPLING_LEAD_MIN = 15;           // C0 drawn this many minutes pre-dose

    // Next scheduled C0 sample time strictly after `after` (a dayjs instance).
    function nextC0Time(after) {
        const base = dayjs(after);
        const candidates = [
            base.hour(DOSE_SLOTS.AM).minute(0),
            base.hour(DOSE_SLOTS.PM).minute(0),
            base.add(1, 'day').hour(DOSE_SLOTS.AM).minute(0)
        ].map(d => d.second(0).millisecond(0).subtract(SAMPLING_LEAD_MIN, 'minute'));
        return candidates.find(c => c.isAfter(base)) || candidates[candidates.length - 1];
    }

    // ============================================================
    // THERAPEUTIC RANGE (time-since-transplant based)
    // ============================================================
    function getTherapeuticRange(date, txDate) {
        const mo = dayjs(date).diff(dayjs(txDate), 'month', true);
        if (mo <= 1) return { low: 10, high: 11 };
        if (mo <= 3) return { low: 7, high: 9 };
        if (mo <= 6) return { low: 5, high: 7 };
        return { low: 4, high: 6 };
    }

    function getLevelStatus(level, date, txDate) {
        const r = getTherapeuticRange(date, txDate);
        if (isNaN(level)) return { text: 'N/A' };
        if (level < r.low) return { text: 'Low' };
        if (level > r.high) return { text: 'High' };
        return { text: 'Target' };
    }

    // ============================================================
    // PK ENGINE — 1-compartment oral Bateman equation
    // All parameters are APPARENT (CL/F, V/F)
    // ============================================================
    function getPopulationParameters(patientData) {
        const { weight, genotype, mpa, bilirubin, inhibitor } = patientData;
        const m = PK_MODEL;
        const cyp = m.CYP3A5[genotype] || m.CYP3A5['unknown'];
        const wtCL = Math.pow(weight / m.WT_REF, m.WT_CL);
        const wtV = Math.pow(weight / m.WT_REF, m.WT_V);
        const mpaF = parseInt(mpa) === 1 ? m.MPA_CL : 1.0;

        // ── HCT REMOVED FROM CL ──────────────────────────────────────────────
        // Størset 2014 is a whole-blood (WB) assay model. HCT variability is
        // already embedded in the population IIV (ωCL, ωV) by design.
        // Re-applying HCT as a structural CL covariate creates double-correction
        // and is absent from the validated parameter set.
        // Use the standalone HCT Correction Tool for assay-level adjustments.

        // ── ALBUMIN REMOVED FROM V ───────────────────────────────────────────
        // Tacrolimus partitions primarily into erythrocytes (~75–80% RBC-bound).
        // Albumin binding is a minor contributor and its effect on apparent V/F
        // in a WB model has no formal exponent in Størset 2014 or validations.
        // Removed to avoid unsubstantiated covariate inflation.

        // ── BILIRUBIN: continuous power-law hepatic impairment model ─────────
        // bilF = (BIL_REF / Bil)^0.30  when Bil > BIL_REF, else 1.0
        // Representative values: Bil=1→1.00 | Bil=2→0.81 | Bil=5→0.62 | Bil=10→0.50
        // Avoids the step-function discontinuity at Bil=2 and Bil=5.
        const bilVal = parseFloat(bilirubin) || m.BIL_REF;
        const bilF = bilVal > m.BIL_REF ? Math.pow(m.BIL_REF / bilVal, 0.30) : 1.0;

        // CYP3A4 Inhibitors. Falls back to 1.0 (no effect) for any value outside
        // none|moderate|strong instead of silently producing NaN — reachable via
        // CSV import or the /api/patients payload, even though the UI itself
        // only ever emits one of the three known keys.
        const inhibitorF = m.INHIBITOR[inhibitor || 'none'] ?? 1.0;

        return {
            CL: m.TVCL * m.INDIAN_CL_SCALAR * cyp * wtCL * mpaF * bilF * inhibitorF,
            V: m.TVV * wtV,
            // KA is population-fixed (4.53 hr⁻¹). Absorption variability (~30% CV)
            // is not individualised by MAP — negligible for trough-only TDM.
            KA: m.TVKA
        };
    }

    // ============================================================
    // TIME-VARYING MATURATION MODEL — REMOVED
    //
    // A heuristic early-post-transplant submodel (bioavailability
    // recovery / oedema / hepatic CL maturation factors) previously
    // lived here. It was disabled after the Sangha 2024 validation
    // (n=103 Indian patients) showed it made predictions worse:
    // MAPE 40% → 53%, F30 65% → 58%. The corrections overcorrected
    // early CL/V for Indian *3/*3 patients.
    //
    // The code has now been deleted rather than left dormant — it was
    // never called, and ~70 lines of detailed documentation for
    // behaviour the app does not have is actively misleading. Recover
    // it from git history if it is ever revisited.
    //
    // Consequence to keep in mind: the model has NO early-phase
    // adjustment, so predictions in the first ~21 days rest on
    // steady-state population PK. updateSummary() surfaces this as the
    // "EARLY POST-TRANSPLANT PHASE" caveat.
    // ============================================================

    function predictSingleDose(timeSinceDose, dose_mg, baseParams, doseTimeHours) {
        // Maturation model disabled: Sangha 2024 validation (n=103 Indian patients)
        // showed it increased MAPE from 40%→53% and reduced F30 from 65%→58%.
        // The heuristic corrections overcorrect early CL/V for Indian *3/*3 patients.
        const { CL, V, KA } = baseParams;
        if (timeSinceDose <= 0 || V <= 0 || dose_mg <= 0) return 0;
        const dose_mcg = dose_mg * 1000;  // mg → mcg; concentration in ng/mL = mcg/L
        const ke = CL / V;
        if (Math.abs(KA - ke) < 0.001) {
            // L'Hôpital / special case
            return (dose_mcg / V) * KA * timeSinceDose * Math.exp(-ke * timeSinceDose);
        }
        return (dose_mcg * KA / (V * (KA - ke))) *
            (Math.exp(-ke * timeSinceDose) - Math.exp(-KA * timeSinceDose));
    }

    function predictAtTime(t, allDoses, params) {
        return allDoses.reduce((sum, d) => {
            if (d.time < t && d.dose > 0) {
                return sum + predictSingleDose(t - d.time, d.dose, params, d.time);
            }
            return sum;
        }, 0);
    }

    // CMIA (Abbott immunoassay) correction applied to LC-MS/MS-scale predictions
    // so they are comparable with observed CMIA values:
    //     CMIA = 1.08 × LC-MS/MS + 0.55
    // The +0.55 intercept is only meaningful where drug is actually present.
    // Applying it at c = 0 (before the first dose, or pre-transplant) produced
    // a spurious 0.55 ng/mL floor on the chart and made the OFV's
    // `predAdj <= 0.1` degenerate-prediction guard unreachable under CMIA.
    function cmiaAdjust(c) {
        return c > 0.01 ? 1.08 * c + 0.55 : c;
    }

    function generateCurve(times, allDoses, params, bioassay) {
        return times.map(t => {
            const c = Math.max(0, predictAtTime(t, allDoses, params));
            return bioassay === 2 ? cmiaAdjust(c) : c;
        });
    }

    // ============================================================
    // NELDER-MEAD SIMPLEX OPTIMIZER (2-dimensional, unconstrained)
    //
    // Used as a gradient-free refinement step after the coarse grid
    // search. Guarantees convergence to the true MAP minimum to
    // machine-precision tolerance, eliminating the ~7% CL error
    // inherent in a 0.15-step grid.
    //
    // Parameters:  α=1.0 (reflect), γ=2.0 (expand),
    //              ρ=0.5 (contract), σ=0.5 (shrink)
    // ============================================================
    function nelderMead2D(fn, x0, y0, { maxIter = 600, tol = 1e-8, step = 0.08 } = {}) {
        const alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;

        // Build initial simplex around (x0, y0)
        let s = [
            [x0, y0, fn(x0, y0)],
            [x0 + step, y0, fn(x0 + step, y0)],
            [x0, y0 + step, fn(x0, y0 + step)]
        ];

        for (let iter = 0; iter < maxIter; iter++) {
            s.sort((a, b) => a[2] - b[2]);          // ascending OFV

            // Convergence: simplex diameter in both dimensions < tol
            if (Math.abs(s[2][0] - s[0][0]) < tol &&
                Math.abs(s[2][1] - s[0][1]) < tol) break;

            // Centroid of best 2 vertices
            const cx = (s[0][0] + s[1][0]) / 2;
            const cy = (s[0][1] + s[1][1]) / 2;

            // Reflection
            const rx = cx + alpha * (cx - s[2][0]);
            const ry = cy + alpha * (cy - s[2][1]);
            const rf = fn(rx, ry);

            if (rf < s[0][2]) {
                // Expansion
                const ex = cx + gamma * (rx - cx);
                const ey = cy + gamma * (ry - cy);
                const ef = fn(ex, ey);
                s[2] = ef < rf ? [ex, ey, ef] : [rx, ry, rf];
            } else if (rf < s[1][2]) {
                s[2] = [rx, ry, rf];
            } else {
                // Contraction toward better of reflected / worst
                const [bx, by, bf] = rf < s[2][2] ? [rx, ry, rf] : s[2];
                const cxc = cx + rho * (bx - cx);
                const cyc = cy + rho * (by - cy);
                const cf = fn(cxc, cyc);
                if (cf < bf) {
                    s[2] = [cxc, cyc, cf];
                } else {
                    // Shrink around best vertex
                    for (let i = 1; i < 3; i++) {
                        s[i][0] = s[0][0] + sigma * (s[i][0] - s[0][0]);
                        s[i][1] = s[0][1] + sigma * (s[i][1] - s[0][1]);
                        s[i][2] = fn(s[i][0], s[i][1]);
                    }
                }
            }
        }

        s.sort((a, b) => a[2] - b[2]);
        return { x: s[0][0], y: s[0][1], fval: s[0][2] };
    }

    // ============================================================
    // MAP BAYESIAN ESTIMATION — 2D (ηCL, ηV) with IIV penalty
    //
    // Three-phase optimizer:
    //   Phase 1 — Coarse global grid: basin identification
    //   Phase 2 — Nelder-Mead simplex: convergence to true minimum
    //   Phase 3 — Outlier / boundary flags: clinical reliability
    //
    // Proportional residual error model (NONMEM-standard WLS):
    //   w = predAdj × SIGMA  (predicted-value weighting, not observed)
    // ============================================================
    //
    // SEARCH BOUNDS vs OUTLIER ALERT — deliberately separate concepts.
    //
    // The search box is wide on purpose: ±2.5 on ηCL is ±6.5 prior SD
    // (ω = √0.148 = 0.385) and ±1.5 on ηV is ±5.0 SD (ω = 0.300). That
    // freedom lets genuine outliers — severe hepatic impairment, strong
    // CYP3A induction — converge without artificial truncation.
    //
    // Because the prior pulls so hard, the search bound is essentially
    // never reached (a patient with true CL at 15% of population, ηCL =
    // −1.90 ≈ 5 SD, converges to −1.71). An earlier revision fired the
    // clinical warning only on a bound hit and described that bound as
    // "±3 SD", so the warning could never appear. The alert now keys off
    // PK_MODEL.OUTLIER_SD (3 SD → ηCL 1.154, ηV 0.900), with a separate,
    // rarer flag for a genuinely truncated fit.
    const ETA_CL_BOUNDS = [-2.5, 2.5];
    const ETA_V_BOUNDS = [-1.5, 1.5];
    const ETA_BOUNDARY_TOL = 0.02;   // treat as truncated within this of a hard bound
    const clamp = (v, [lo, hi]) => Math.min(hi, Math.max(lo, v));

    function mapBayesian(popParams, measuredLevels, allDoses, bioassay) {
        const { OMEGA_CL, OMEGA_V, SIGMA } = PK_MODEL;

        if (measuredLevels.length === 0) {
            return {
                ...popParams, etaCL: 0, etaV: 0,
                rmse: null, mpe: null, mape: null, r2: null, pe: [],
                boundaryHit: false
            };
        }

        // ── MAP objective function ────────────────────────────────────────
        // OFV = ηCL²/ωCL + ηV²/ωV  +  Σ[(Cobs − Cpred)² / (Cpred·σ)²]
        const ofv = (eCL, eV) => {
            const p = {
                CL: popParams.CL * Math.exp(eCL),
                V: popParams.V * Math.exp(eV),
                KA: popParams.KA
            };
            const prior = (eCL * eCL) / OMEGA_CL + (eV * eV) / OMEGA_V;
            const lik = measuredLevels.reduce((s, obs) => {
                const pred = predictAtTime(obs.time, allDoses, p);
                const predAdj = bioassay === 2 ? cmiaAdjust(pred) : pred;
                if (predAdj <= 0.1) return s + 1000;
                const w = predAdj * SIGMA;   // WLS: predicted-value denominator
                return s + Math.pow((obs.level - predAdj) / w, 2);
            }, 0);
            return prior + lik;
        };

        // ── Phase 1: Coarse global grid — basin identification ─────────────
        // Step 0.25 across full eta space to locate the basin of attraction.
        // Kept coarse intentionally; Nelder-Mead refines from here.
        // Driven by an integer step count: accumulating `e += 0.25` in a float
        // loop drifts, so the final grid point could fall just past the bound
        // and be skipped, making the search grid asymmetric.
        const GRID_STEP = 0.25;
        const nCL = Math.round((ETA_CL_BOUNDS[1] - ETA_CL_BOUNDS[0]) / GRID_STEP);
        const nV = Math.round((ETA_V_BOUNDS[1] - ETA_V_BOUNDS[0]) / GRID_STEP);
        let bestCL = 0, bestV = 0, minOFV = Infinity;
        for (let i = 0; i <= nCL; i++) {
            const eCL = ETA_CL_BOUNDS[0] + i * GRID_STEP;
            for (let j = 0; j <= nV; j++) {
                const eV = ETA_V_BOUNDS[0] + j * GRID_STEP;
                const v = ofv(eCL, eV);
                if (v < minOFV) { minOFV = v; bestCL = eCL; bestV = eV; }
            }
        }

        // ── Phase 2: Nelder-Mead simplex refinement ───────────────────────
        // Starts at the grid winner; converges to sub-1e-8 precision,
        // eliminating the ~7% CL error from a 0.15-step grid.
        // nelderMead2D is UNCONSTRAINED, so clamp back into the search box —
        // otherwise indParams could carry an eta the grid would never allow.
        const nm = nelderMead2D(ofv, bestCL, bestV);
        bestCL = clamp(nm.x, ETA_CL_BOUNDS);
        bestV = clamp(nm.y, ETA_V_BOUNDS);

        // ── Phase 3a: Outlier alert (clinical) ───────────────────────────
        // Fires when the converged fit is beyond OUTLIER_SD prior SDs.
        // Thresholds derived from ω so they stay correct if the population
        // variances are ever recalibrated: 3·√0.148 = 1.154, 3·√0.090 = 0.900.
        const outlierCL = PK_MODEL.OUTLIER_SD * Math.sqrt(OMEGA_CL);
        const outlierV = PK_MODEL.OUTLIER_SD * Math.sqrt(OMEGA_V);
        const sdCL = bestCL / Math.sqrt(OMEGA_CL);
        const sdV = bestV / Math.sqrt(OMEGA_V);

        let outlierMsg = '';
        if (bestCL <= -outlierCL) outlierMsg += `ηCL ${sdCL.toFixed(1)} SD below population (unusually fast clearance). `;
        if (bestCL >= outlierCL) outlierMsg += `ηCL ${sdCL.toFixed(1)} SD above population (unusually slow clearance). `;
        if (bestV <= -outlierV) outlierMsg += `ηV ${sdV.toFixed(1)} SD below population (unusually small V). `;
        if (bestV >= outlierV) outlierMsg += `ηV ${sdV.toFixed(1)} SD above population (unusually large V). `;
        const outlierHit = outlierMsg.length > 0;

        // ── Phase 3b: Search-bound truncation (numerical) ────────────────
        // Distinct from the above: this means the optimiser was actually
        // stopped by the box, so the estimate is not a free optimum.
        const hitCLLo = bestCL <= ETA_CL_BOUNDS[0] + ETA_BOUNDARY_TOL;
        const hitCLHi = bestCL >= ETA_CL_BOUNDS[1] - ETA_BOUNDARY_TOL;
        const hitVLo = bestV <= ETA_V_BOUNDS[0] + ETA_BOUNDARY_TOL;
        const hitVHi = bestV >= ETA_V_BOUNDS[1] - ETA_BOUNDARY_TOL;
        const boundaryHit = hitCLLo || hitCLHi || hitVLo || hitVHi;

        let boundaryMsg = outlierMsg;
        if (boundaryHit) boundaryMsg += 'Estimate was truncated at the search boundary — treat the parameter values themselves as unreliable. ';

        const indParams = {
            CL: popParams.CL * Math.exp(bestCL),
            V: popParams.V * Math.exp(bestV),
            KA: popParams.KA,
            etaCL: bestCL,
            etaV: bestV,
            outlierHit,
            boundaryHit,
            boundaryMsg: boundaryMsg.trim()
        };
        const metrics = calculateAccuracyMetrics(measuredLevels, allDoses, indParams, bioassay);
        return { ...indParams, ...metrics };
    }

    // ============================================================
    // ACCURACY METRICS
    // RMSE, MPE (bias), MAPE (imprecision), R²
    // ============================================================
    function calculateAccuracyMetrics(measuredLevels, allDoses, params, bioassay) {
        if (measuredLevels.length === 0) return { rmse: null, mpe: null, mape: null, r2: null, pe: [] };
        const pe = measuredLevels.map((obs, i) => {
            const pred = (() => {
                const p = predictAtTime(obs.time, allDoses, params);
                return bioassay === 2 ? cmiaAdjust(p) : p;
            })();
            const error = obs.level - pred;
            const pct = obs.level > 0 ? (error / obs.level) * 100 : 0;
            return { n: i + 1, time: obs.time, observed: obs.level, predicted: pred, error, pct };
        });
        const n = pe.length;
        const rmse = Math.sqrt(pe.reduce((s, p) => s + p.error * p.error, 0) / n);
        const mpe = pe.reduce((s, p) => s + p.pct, 0) / n;
        const mape = pe.reduce((s, p) => s + Math.abs(p.pct), 0) / n;
        const meanObs = measuredLevels.reduce((s, o) => s + o.level, 0) / n;
        const ssTot = measuredLevels.reduce((s, o) => s + Math.pow(o.level - meanObs, 2), 0);
        const ssRes = pe.reduce((s, p) => s + p.error * p.error, 0);
        const r2 = ssTot > 1e-6 ? Math.max(0, 1 - ssRes / ssTot) : null;
        return { rmse, mpe, mape, r2, pe };
    }

    // ============================================================
    // LAPLACE POSTERIOR COVARIANCE
    //
    // Exact (to second order) posterior covariance of (ηCL, ηV) at the
    // MAP estimate:
    //
    //     Σ_post = (Ω⁻¹ + JᵀJ / σ²)⁻¹
    //
    // where J[i] = [∂ln C(tᵢ)/∂ηCL, ∂ln C(tᵢ)/∂ηV] is the sensitivity of
    // each PREDICTED observation to a log-scale change in the parameter.
    // Differentiating on the η (log) scale is essential — using absolute
    // CL/V produces a matrix in the wrong units and can yield a "posterior"
    // wider than the prior, which is impossible.
    //
    // This replaces an earlier closed-form heuristic
    //     var_post = Ω / (1 + Ω·nObs/σ²),  ρ = −0.45 (hardcoded)
    // which assumed every observation informs CL and V equally with unit
    // sensitivity. It does not. For a steady-state trough:
    //     ∂lnC/∂lnCL = −1.20      ∂lnC/∂lnV = +0.20
    // A trough carries almost no information about V, and because the two
    // sensitivities have OPPOSITE signs the posterior correlation is
    // POSITIVE, not negative. (Raising V lowers ke = CL/V and therefore
    // RAISES the trough, so CL and V move together along the ridge of
    // equally good fits.)
    //
    // Exact values for this model (3mg BID steady state, C0 at 06:45;
    // prior sd(ηCL) = 0.385, sd(ηV) = 0.300):
    //     n=1 trough   sd 0.147 / 0.298   ρ = +0.285
    //     n=3 troughs  sd 0.097 / 0.298   ρ = +0.475
    //     n=5 troughs  sd 0.081 / 0.298   ρ = +0.576
    //     n=4 mixed    sd 0.085 / 0.267   ρ = +0.166   (C0 + C2 peak)
    // sd(ηV) barely moves off the prior because troughs do not identify V —
    // corroborated by MAP round-trip tests (true ηV = +0.40 → MAP +0.03).
    // The old heuristic reported sd(ηV) = 0.154, understating it ~2×, so
    // the plotted 90% band was far narrower than the data supports.
    //
    // Sanity guard: JᵀJ/σ² is positive semi-definite, so Σ_post ⪯ Ω always.
    // A posterior sd exceeding its prior sd is a bug, never a result.
    // ============================================================
    function laplacePosterior(measuredLevels, allDoses, indParams, bioassay) {
        const { OMEGA_CL, OMEGA_V, SIGMA } = PK_MODEL;
        const sigma2 = SIGMA * SIGMA;
        const priorSdCL = Math.sqrt(OMEGA_CL);
        const priorSdV = Math.sqrt(OMEGA_V);
        const prior = { sdCL: priorSdCL, sdV: priorSdV, rho: 0 };

        if (!measuredLevels || measuredLevels.length === 0) return prior;

        // Predicted (bioassay-adjusted) concentration at time t for an
        // (ηCL, ηV) perturbation of the individual estimate.
        const predAt = (t, eCL, eV) => {
            const p = {
                CL: indParams.CL * Math.exp(eCL),
                V: indParams.V * Math.exp(eV),
                KA: indParams.KA
            };
            const c = predictAtTime(t, allDoses, p);
            return bioassay === 2 ? cmiaAdjust(c) : c;
        };

        // Accumulate JᵀJ by central differences on the log scale.
        const h = 1e-5;
        let Sxx = 0, Sxy = 0, Syy = 0;
        for (const obs of measuredLevels) {
            const c0 = predAt(obs.time, 0, 0);
            if (!(c0 > 0.1)) continue;   // uninformative / degenerate point
            const gx = (predAt(obs.time, h, 0) - predAt(obs.time, -h, 0)) / (2 * h) / c0;
            const gy = (predAt(obs.time, 0, h) - predAt(obs.time, 0, -h)) / (2 * h) / c0;
            if (!isFinite(gx) || !isFinite(gy)) continue;
            Sxx += gx * gx;
            Sxy += gx * gy;
            Syy += gy * gy;
        }

        // Σ_post = (Ω⁻¹ + JᵀJ/σ²)⁻¹  — closed-form 2×2 inverse
        const a = 1 / OMEGA_CL + Sxx / sigma2;
        const b = Sxy / sigma2;
        const c = 1 / OMEGA_V + Syy / sigma2;
        const det = a * c - b * b;
        if (!isFinite(det) || det <= 0) return prior;

        const varCL = c / det;
        const varV = a / det;
        const cov = -b / det;
        if (!(varCL > 0) || !(varV > 0)) return prior;

        // Clamp to the prior. Mathematically Σ_post ⪯ Ω, so this can only
        // trigger on numerical noise — but it makes the invariant explicit.
        const sdCL = Math.min(Math.sqrt(varCL), priorSdCL);
        const sdV = Math.min(Math.sqrt(varV), priorSdV);
        let rho = cov / Math.sqrt(varCL * varV);
        if (!isFinite(rho)) rho = 0;
        rho = Math.max(-0.99, Math.min(0.99, rho));

        return { sdCL, sdV, rho };
    }

    // ============================================================
    // MONTE CARLO CONFIDENCE INTERVALS (90%)
    //
    //   1. Box-Muller: independent (u1, u2) per variate pair, so z1 and z2
    //      are genuinely independent standard-normals before correlation.
    //
    //   2. Posterior draws use the Laplace covariance above, Cholesky-
    //      factored:  z_CL = z1 ;  z_V = ρ·z1 + sqrt(1−ρ²)·z2
    //
    //   3. CI floor at t < 0.5 hr: prevents the band from collapsing to
    //      zero at the transplant time origin, which is mathematically
    //      correct but visually misleading to clinical readers.
    //
    // NOTE: this is a CONFIDENCE interval on the model prediction
    // (parameter uncertainty). It does not include the σ=18% residual
    // error, so it is narrower than a prediction interval for a future
    // assay result.
    // ============================================================
    function monteCarloCI(timePoints, allDoses, indParams, measuredLevels, bioassay, nSamples = 300) {
        const { sdCL, sdV, rho } = laplacePosterior(measuredLevels, allDoses, indParams, bioassay);
        const rhoPerp = Math.sqrt(1 - rho * rho);   // sqrt(1 - ρ²)

        const curves = [];
        for (let s = 0; s < nSamples; s++) {
            // ── Fix 1: truly independent Box-Muller pairs ─────────────────
            // Each variate draws its OWN (u1, u2) so z_CL and z_V
            // are independent standard-normals before correlation.
            const ua1 = Math.random() + 1e-10;
            const ua2 = Math.random();
            const ub1 = Math.random() + 1e-10;
            const ub2 = Math.random();
            const z1 = Math.sqrt(-2 * Math.log(ua1)) * Math.cos(2 * Math.PI * ua2);
            const z2 = Math.sqrt(-2 * Math.log(ub1)) * Math.cos(2 * Math.PI * ub2);

            // ── Fix 2: Cholesky-correlated posterior samples ──────────────
            const etaCL = z1 * sdCL;
            const etaV = (rho * z1 + rhoPerp * z2) * sdV;

            const p = {
                CL: indParams.CL * Math.exp(etaCL),
                V: indParams.V * Math.exp(etaV),
                KA: indParams.KA
            };
            curves.push(generateCurve(timePoints, allDoses, p, bioassay));
        }

        // ── Fix 4: CI floor at t < 0.5 hr ────────────────────────────
        // At t ≈ 0 all curves are near-zero; the collapsed CI band
        // is mathematically correct but can mislead clinical staff
        // into thinking there is no uncertainty at that time.
        const FLOOR_TIME_HR = 0.5;   // apply floor for t < 0.5 hr

        return timePoints.map((t, i) => {
            const vals = curves.map(c => c[i]).sort((a, b) => a - b);
            const lo = Math.floor(0.05 * nSamples);
            const hi = Math.floor(0.95 * nSamples);
            let p5 = vals[lo];
            let p95 = vals[hi];
            if (t < FLOOR_TIME_HR && p95 - p5 < 0.1) {
                // Expand to a nominal ±0.05 ng/mL floor so the band
                // is always visible without distorting real CI data.
                p5 = Math.max(0, p5 - 0.05);
                p95 = p95 + 0.05;
            }
            return { p5, p95 };
        });
    }

    // ============================================================
    // TTR — ROSENDAAL METHOD (correct linear interpolation)
    // ============================================================
    function rosendaalFraction(v1, v2, low, high) {
        // Fraction of interval [0,1] where linearly interpolated value
        // v(s) = v1 + s*(v2-v1) lies in [low, high]
        if (Math.abs(v2 - v1) < 1e-10) return (v1 >= low && v1 <= high) ? 1.0 : 0.0;
        const s_low = (low - v1) / (v2 - v1);
        const s_high = (high - v1) / (v2 - v1);
        const sStart = Math.max(0, Math.min(s_low, s_high));
        const sEnd = Math.min(1, Math.max(s_low, s_high));
        return Math.max(0, sEnd - sStart);
    }

    // ============================================================
    // DOSE INTERPOLATION & EXTRAPOLATION
    // ============================================================
    function fillHistoricalGaps(historyLog, txDate) {
        const filled = [];
        // Anchor on dose >= 0, NOT dose > 0. A recorded 0 mg is a deliberate
        // "dose held" entry. Filtering it out dropped it from the anchor list,
        // and the 12-hourly gap-filler then spanned straight across it and
        // re-inserted a full dose in that slot — so recording a held dose had
        // the exact opposite of the intended effect.
        const logged = historyLog.filter(e => e.dose !== null && e.dose >= 0).sort((a, b) => a.time - b.time);
        if (logged.length === 0) return [];

        for (let i = 0; i < logged.length; i++) {
            filled.push(logged[i]);
            if (i === logged.length - 1) break;

            let nextTime = logged[i].recordDate.add(12, 'hour');
            // Fill 12-hour steps until we hit the next logged dose (with a 4hr buffer to avoid overlapping)
            while (nextTime.isBefore(logged[i + 1].recordDate.subtract(4, 'hour'))) {
                const isAm = nextTime.hour() < 12;
                // Carry the last same-shift POSITIVE dose. A single held dose
                // marks that one slot only; it must not zero the rest of the
                // regimen for every subsequent imputed slot.
                const priorSameShift = [...filled].reverse()
                    .find(d => (d.recordDate.hour() < 12) === isAm && d.dose > 0);
                const doseToGive = priorSameShift ? priorSameShift.dose : logged[i].dose;
                filled.push({ id: 'inter-' + i + '-' + nextTime.valueOf(), recordDate: nextTime, dose: doseToGive, level: null, time: nextTime.diff(txDate, 'hour', true) });
                nextTime = nextTime.add(12, 'hour');
            }
        }
        return filled;
    }

    function extrapolateDoses(filledHistory, untilDate, txDate) {
        if (filledHistory.length === 0) return [];
        const last = filledHistory[filledHistory.length - 1];
        const extra = [];
        let nextTime = last.recordDate.add(12, 'hour');
        while (nextTime.isBefore(untilDate)) {
            const isAm = nextTime.hour() < 12;
            // Same rule as fillHistoricalGaps: carry the last same-shift
            // POSITIVE dose, so a held dose at the end of the log does not
            // extrapolate a zero regimen indefinitely into the future.
            const priorSameShift = [...filledHistory].reverse()
                .find(d => (d.recordDate.hour() < 12) === isAm && d.dose > 0);
            const doseToGive = priorSameShift ? priorSameShift.dose : last.dose;
            extra.push({ id: 'ext-' + extra.length, recordDate: nextTime, dose: doseToGive, level: null, time: nextTime.diff(txDate, 'hour', true) });
            nextTime = nextTime.add(12, 'hour');
        }
        return extra;
    }

    // ============================================================
    // BRIDGE DOSES — cover a lagging dose log
    //
    // The dose optimizer can only apply a NEW regimen from the first dosing slot
    // the clinician can still act on. But the patient did not stop taking
    // tacrolimus while the log fell behind — they continued the existing
    // regimen. Everything in [fromDate, toDate) must therefore be filled in.
    //
    // Omitting this modelled a drug-free washout and made the optimizer
    // recommend roughly double the correct dose. Measured at 3mg BID steady
    // state, target = next trough: -22% predicted trough at a 1-day logging
    // lag, -58% at 3 days, -66% at 5 days.
    //
    // Same last-same-shift-POSITIVE-dose rule as fillHistoricalGaps, so a held
    // dose does not propagate a zero regimen across the whole bridge.
    // ============================================================
    function buildBridgeDoses(logged, fromDate, toDate, txDate) {
        const bridge = [];
        if (!logged || logged.length === 0) return bridge;
        const fallback = logged[logged.length - 1].dose;
        for (let t = fromDate; t.isBefore(toDate); t = t.add(12, 'hour')) {
            const isAm = t.hour() < 12;
            const priorSameShift = [...logged, ...bridge].reverse()
                .find(d => (d.recordDate.hour() < 12) === isAm && d.dose > 0);
            bridge.push({
                id: 'bridge-' + bridge.length,
                recordDate: t,
                dose: priorSameShift ? priorSameShift.dose : fallback,
                level: null,
                time: t.diff(txDate, 'hour', true)
            });
        }
        return bridge;
    }

    // ============================================================
    // Flags levels that look like they were drawn AFTER a dose rather than as a
    // pre-dose C0. Absorption is steep (tmax ≈ 1.2h), so a sample 15 min the
    // wrong side of a dose reads ~30% high; MAP then infers slower clearance and
    // under-doses the patient. Reported only — the level is still fitted exactly
    // as entered, because silently rewriting a clinician's sample time is not
    // acceptable in a record that drives dosing.
    function findPostDoseSamples(measuredLevels, allDoses) {
        const WINDOW_HR = 2;
        return (measuredLevels || []).map((obs, i) => {
            const priorDose = allDoses
                .filter(d => d.dose > 0 && d.time < obs.time && d.time >= obs.time - WINDOW_HR)
                .sort((a, b) => b.time - a.time)[0];
            return priorDose ? { n: i + 1, minutesAfter: Math.round((obs.time - priorDose.time) * 60) } : null;
        }).filter(Boolean);
    }


    return {
        // constants
        PK_MODEL, DOSE_SLOTS, SAMPLING_LEAD_MIN,
        ETA_CL_BOUNDS, ETA_V_BOUNDS, ETA_BOUNDARY_TOL,
        // time / protocol helpers
        nextC0Time, getTherapeuticRange, getLevelStatus,
        // PK model
        getPopulationParameters, predictSingleDose, predictAtTime,
        cmiaAdjust, generateCurve,
        // estimation
        nelderMead2D, mapBayesian, calculateAccuracyMetrics,
        laplacePosterior, monteCarloCI,
        // dose series
        fillHistoricalGaps, extrapolateDoses, buildBridgeDoses,
        // analytics / QC
        rosendaalFraction, findPostDoseSamples
    };
}));
