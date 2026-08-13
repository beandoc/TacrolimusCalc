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

        // ── Recency weighting (forgetting-factor MAP) ─────────────────────
        // Each observation's contribution to the fit decays exponentially
        // with its age relative to the MOST RECENT level in the fitted set,
        // not wall-clock "today" — so a forecast run long after the last
        // draw doesn't retroactively discount every observation equally.
        // This is the batch equivalent of a sequential/Kalman filter with a
        // forgetting factor: recomputing fresh each call (rather than
        // chaining posteriors visit-to-visit) is deliberate, because patient
        // history here is editable, and a chained filter would need a full
        // replay on every correction anyway — see mapBayesian.
        //
        // ONE global half-life cannot serve two genuinely different regimes,
        // confirmed on a real held-out-forecast test (troughs 4.5, 12.9,
        // 11.1, 7.4, 6.2 over ~3 weeks post-transplant; actual next trough
        // 6.8): a 45-day half-life left the forecast at 11.3 ng/mL (66%
        // error, barely moved from the unweighted 11.6) — the maturation
        // happened faster than that decay discounted it. A half-life short
        // enough to fix that (~5-10d) would make any single new level
        // dominate the fit for a stable, quarterly-monitored long-term
        // patient too — a real regression, not a hypothetical one.
        //
        // Split into named REGIMES instead. Which one a patient gets is
        // decided empirically per patient by selectRecencyRegime() — by
        // backtesting both against that patient's OWN recent troughs —
        // not by a global rule.
        //
        // A fixed post-transplant-day threshold was tried first and is kept
        // only as the cold-start fallback (< MIN_LEVELS_FOR_REGIME_SELECT
        // levels, nothing to backtest on yet). It was wrong often enough to
        // abandon as the primary: on a 92-prediction held-out backtest over
        // 10 real patients it improved the two genuinely-drifting patients a
        // lot (arun_chougle 107%->64% MAPE, raj_bahadur 71%->42%) but made
        // vidyashree materially WORSE (28%->37%) — her troughs oscillate
        // noisily (1.5, 2.8, 4.4, 6.3, 8.9, 11.4, 8.4 ...) without trending,
        // all inside the "early" window, so 'adaptive' chased noise as if it
        // were a new steady state. Post-transplant day says nothing about
        // whether THIS patient's clearance is actually moving; their own
        // prediction history does.
        //
        // 'adaptive' half-life/floor swept against a real held-out test
        // (3-week span, actual next trough 6.8 ng/mL):
        //   10d/0.05 -> 10.4 (53% err)   7d/0.05 -> 9.6 (41% err)
        //   5d/0.05  ->  8.5 (26% err)   7d/0    -> 9.1 (33% err)
        //   5d/0     ->  7.8 (15% err)   3d/0    -> 7.1  (5% err)
        // 3d/0 fits that one case almost exactly, which is itself reason to
        // be wary — that precision is as likely overfit as correct. 5d/0 was
        // chosen as a clear improvement without chasing one patient's best
        // number.
        RECENCY_REGIMES: {
            // Clearance may genuinely be moving now: discount hard by age.
            adaptive: { halfLifeDays: 5, floor: 0 },
            // Assume stability; covariate tracking (inhibitors, weight)
            // already handles the known reasons clearance would move. Long
            // memory so one new level can't override an established fit.
            stationary: { halfLifeDays: 180, floor: 0.4 }
        },
        // Cold-start fallback only — see selectRecencyRegime.
        MATURATION_WINDOW_DAYS: 90,
        // Below this many levels there is nothing to backtest a regime on
        // (each trial needs >=3 prior levels to fit plus a held-out target).
        MIN_LEVELS_FOR_REGIME_SELECT: 5,
        // How many of the most recent levels to hold out when scoring a
        // regime. Three was too jumpy on real backend data: it let a noisy
        // short window flip drifting patients back to stationary. Five gives
        // the selector more memory while still staying cheap enough for the UI.
        REGIME_SELECT_TRIALS: 5,

        // Outlier alert threshold, in prior SD units. Deliberately DECOUPLED from
        // the MAP search bounds (ETA_*_BOUNDS): the search stays wide so genuine
        // outliers (severe hepatic impairment, strong CYP3A induction) converge
        // without artificial truncation, while the clinical banner fires at 3 SD.
        //   3 SD on ηCL = 3·√0.148 = 1.154
        //   3 SD on ηV  = 3·√0.090 = 0.900
        OUTLIER_SD: 3,

        // CYP3A5 covariate on CL/F
        //
        // 'unknown' = 1.30 is the ALLELE-FREQUENCY-WEIGHTED POPULATION MEAN, not
        // a midpoint and not a safety margin. At p(*1) = 0.45 in the Indian
        // population, Hardy-Weinberg gives *1/*1 0.2025, *1/*3 0.4950, *3/*3
        // 0.3025, so the expected factor is
        //     0.2025(1.79) + 0.4950(1.35) + 0.3025(1.00) = 1.33.
        //
        // Note where that lands: 1.30 is BELOW the *1/*3 value of 1.35 — the
        // second-lowest CL on this list, not an "intermediate-high, conservative"
        // default as an earlier comment claimed. Lower assumed CL means a higher
        // predicted trough and therefore a LOWER recommended dose, so for an
        // ungenotyped patient in a population where ~70% carry a *1 allele this
        // errs toward under-exposure (rejection), not toxicity. Genotype the
        // patient rather than relying on this default where that matters.
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
    //
    // TDM sampling is MORNING ONLY. C0 is drawn SAMPLING_LEAD_MIN minutes before
    // the 07:00 dose — 06:45 — and never before the evening dose. This is site
    // protocol, not a modelling simplification: the evening trough is never
    // assayed, so offering 18:45 as a default sample time invited timestamps for
    // draws that do not happen. (Levels ARE still fitted at whatever time they
    // are entered — see findPostDoseSamples — this governs defaults only.)
    //
    // Why the 15 min matters numerically: predictAtTime() includes a dose only
    // when dose.time < t (strict), so evaluating at exactly 07:00 happens to
    // give a pre-dose value — but it sits precisely on a discontinuity. At 3mg
    // BID steady state:
    //     06:45 (true C0)          →  9.15 ng/mL
    //     07:00 (on the boundary)  →  9.08 ng/mL   (−0.8%, harmless)
    //     07:01 (one minute late)  →  9.39 ng/mL   (+2.6%)
    //     07:15 (15 min late)      → 11.92 ng/mL   (+30.2%)
    // Defaulting every timestamp to 06:45 both matches the protocol and keeps
    // every evaluation a safe 15 minutes clear of that cliff.
    // ============================================================
    const DOSE_SLOTS = { AM: 7, PM: 19 };   // hour-of-day for BD dosing
    const SAMPLING_LEAD_MIN = 15;           // C0 drawn this many minutes pre-dose

    // Next scheduled C0 sample time strictly after `after` (a dayjs instance).
    // Always a MORNING draw (06:45) — today's if it is still ahead, else
    // tomorrow's. See the protocol note above.
    function nextC0Time(after) {
        const base = dayjs(after);
        const morning = d => d.hour(DOSE_SLOTS.AM).minute(0).second(0).millisecond(0)
            .subtract(SAMPLING_LEAD_MIN, 'minute');
        const today = morning(base);
        return today.isAfter(base) ? today : morning(base.add(1, 'day'));
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
            KA: m.TVKA,
            weight: parseFloat(weight) || m.WT_REF,
            baselineInhibitor: inhibitor || 'none'
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
                let doseParams = params;
                const wtRatio = (d.weight != null && params.weight != null) ? (d.weight / params.weight) : 1.0;
                const baseInh = params.baselineInhibitor || 'none';
                const dInh = (d.inhibitor !== undefined && d.inhibitor !== null && d.inhibitor !== '') ? d.inhibitor : baseInh;

                if (wtRatio !== 1.0 || dInh !== baseInh) {
                    const wtCLScale = Math.pow(wtRatio, PK_MODEL.WT_CL);
                    const wtVScale = Math.pow(wtRatio, PK_MODEL.WT_V);

                    const baseInhF = PK_MODEL.INHIBITOR[baseInh] ?? 1.0;
                    const doseInhF = PK_MODEL.INHIBITOR[dInh] ?? 1.0;
                    const inhScale = doseInhF / baseInhF;

                    doseParams = {
                        ...params,
                        CL: params.CL * wtCLScale * inhScale,
                        V: params.V * wtVScale
                    };
                }
                return sum + predictSingleDose(t - d.time, d.dose, doseParams, d.time);
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
    // STARTING DOSE SUGGESTION — population-based, pre-Bayesian
    //
    // Before any level is drawn there is nothing for MAP to fit, so dosing has
    // to come from the population model alone. CPIC (Birdwell 2015) recommends
    // a flat 1.5-2x starting-dose increase for CYP3A5 expressers vs
    // non-expressers; this reproduces that guidance quantitatively from the
    // SAME population model used for the MAP Bayesian forecast (CYP3A5 factor,
    // weight, MPA, bilirubin, CYP3A4 inhibitor already applied via
    // getPopulationParameters), rather than a flat CPIC multiplier disconnected
    // from those covariates.
    //
    // Grid-searches BID total daily dose and simulates enough 12h-spaced doses
    // to reach steady state AT THIS PATIENT'S modelled half-life (not a fixed
    // dose count) — a strong CYP3A4 inhibitor or severe hepatic impairment can
    // quarter clearance and push steady state out by days, and a fixed count
    // tuned for the typical case would silently under-simulate those patients.
    //
    // `bioassay` matters and must be passed. The target range is expressed in
    // whatever the lab REPORTS, so the simulated trough has to be put on the same
    // scale before it is compared — exactly as mapBayesian and the dose optimizer
    // already do. Omitting it targeted the LC-MS/MS value: for a *3/*3 60 kg
    // patient aiming at the 10-11 band this returned 3 mg BID for a predicted
    // 10.58, which a CMIA lab reports as 1.08 x 10.58 + 0.55 = 11.98 — above the
    // band it was aiming for, at every CMIA centre.
    // ============================================================
    const STARTING_DOSE_MIN_TDD = 1.0;
    const STARTING_DOSE_MAX_TDD = 18.0;
    const STARTING_DOSE_STEP = 0.5;
    const STARTING_DOSE_HALF_LIVES = 6; // 2^-6 = 1.6% from true steady state
    const STARTING_DOSE_MAX_DOSES = 200; // 100-day safety cap on the simulation

    function suggestStartingDose(patientData, targetRange, bioassay) {
        const popParams = getPopulationParameters(patientData);
        const targetMid = (targetRange.low + targetRange.high) / 2;
        // Fall back to the profile's own setting so a caller that forgets the
        // argument still gets the patient's assay rather than silently LC-MS/MS.
        const assay = bioassay != null ? bioassay : patientData.bioassay;

        const halfLifeHr = Math.LN2 / (popParams.CL / popParams.V);
        // Dose index 0 is an AM dose, so even indices are AM and odd are PM.
        // The trough is evaluated immediately before dose `nDoses - 1`, which
        // must therefore be an AM dose for that trough to be the MORNING C0 the
        // patient is actually sampled at (site protocol: draws are 06:45 only,
        // never 18:45 — see the DOSING & SAMPLING CONVENTION note above).
        //
        // Rounding nDoses up to an odd count makes `nDoses - 1` even. Without
        // this the parity fell out of the patient's half-life, so an asymmetric
        // split was targeted at the pre-PM trough for some patients and the
        // pre-AM trough for others — 2.6% apart at 3.5 mg AM / 3.0 mg PM.
        let nDoses = Math.min(
            STARTING_DOSE_MAX_DOSES,
            Math.max(10, Math.ceil((STARTING_DOSE_HALF_LIVES * halfLifeHr) / 12) + 1)
        );
        if (nDoses % 2 === 0) nDoses += 1;
        const lastDoseTime = (nDoses - 1) * 12;

        let best = null;
        for (let tdd = STARTING_DOSE_MIN_TDD; tdd <= STARTING_DOSE_MAX_TDD + 1e-9; tdd += STARTING_DOSE_STEP) {
            // Standard tablet split: round each half to 0.5mg, remainder to AM.
            const amDose = Math.round((tdd / 2) * 2) / 2;
            const pmDose = Math.round((tdd - amDose) * 2) / 2;

            const doses = [];
            for (let i = 0; i < nDoses; i++) {
                doses.push({ time: i * 12, dose: i % 2 === 0 ? amDose : pmDose });
            }
            // Trough = concentration right before the final dose, which nDoses'
            // odd count guarantees is an AM dose — i.e. the 06:45 morning C0.
            // The strict `d.time < t` in predictAtTime excludes that dose itself.
            // Put it on the REPORTED assay scale before comparing with the target.
            const raw = predictAtTime(lastDoseTime, doses, popParams);
            const trough = assay === 2 ? cmiaAdjust(raw) : raw;
            const diff = Math.abs(trough - targetMid);
            if (!best || diff < best.diff) {
                best = { tdd, amDose, pmDose, trough, diff };
            }
        }

        return {
            tdd: best.tdd,
            amDose: best.amDose,
            pmDose: best.pmDose,
            predictedTrough: best.trough,
            targetRange,
            halfLifeHr,
            atBoundary: best.tdd <= STARTING_DOSE_MIN_TDD || best.tdd >= STARTING_DOSE_MAX_TDD,
            popParams
        };
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

    // ============================================================
    // RECENCY WEIGHTS — forgetting-factor equivalent, batch form.
    // Age is relative to the MOST RECENT level in the set passed in, not
    // wall-clock "today": see PK_MODEL.RECENCY_HALFLIFE_DAYS. Exported so
    // the UI can show a clinician which levels are actually driving the
    // current fit.
    // ============================================================
    // Rescales weights to sum to N (their count), preserving every RATIO
    // between individual weights while keeping their total magnitude
    // constant. Without this, downweighting old/disagreeing observations
    // shrinks the total likelihood term relative to the FIXED prior
    // penalty (eCL²/ωCL + eV²/ωV), which pulls even a clean, non-drifting
    // fit toward the population prior just because some of its weight sum
    // is below N — confirmed by a noise-free round-trip test recovering
    // ηCL=0.67 instead of the true 0.8 before this was added. Normalizing
    // keeps a recent point's weight relatively higher than an old point's
    // (the actual goal) without deflating how much total evidence the
    // patient's data carries against the prior.
    function normalizeWeights(weights) {
        if (!weights || weights.length === 0) return weights;
        const total = weights.reduce((s, w) => s + w, 0);
        if (!(total > 0)) return weights;
        const scale = weights.length / total;
        return weights.map(w => w * scale);
    }

    // Cold-start regime guess, used ONLY when there aren't yet enough levels
    // for selectRecencyRegime to measure which regime actually works for
    // this patient. Post-transplant day is a weak proxy — it says when
    // clearance COULD be moving, never whether it is — so this is a
    // starting assumption to be overridden by evidence, not a rule.
    function defaultRegimeFor(measuredLevels) {
        if (!measuredLevels || measuredLevels.length === 0) return 'stationary';
        const mostRecentTime = Math.max(...measuredLevels.map(o => o.time));
        return (mostRecentTime / 24) < PK_MODEL.MATURATION_WINDOW_DAYS ? 'adaptive' : 'stationary';
    }

    function recencyWeights(measuredLevels, regime) {
        if (!measuredLevels || measuredLevels.length === 0) return [];
        const name = regime || defaultRegimeFor(measuredLevels);
        const cfg = PK_MODEL.RECENCY_REGIMES[name] || PK_MODEL.RECENCY_REGIMES.stationary;
        const mostRecentTime = Math.max(...measuredLevels.map(o => o.time));
        // The regime is fixed for the whole fit, so weights within a single
        // fit stay a smooth function of age with no discontinuity.
        return normalizeWeights(measuredLevels.map(obs => {
            const ageDays = (mostRecentTime - obs.time) / 24;
            const decay = Math.exp(-Math.LN2 * ageDays / cfg.halfLifeDays);
            return cfg.floor + (1 - cfg.floor) * decay;
        }));
    }

    // ============================================================
    // ROBUST (RESIDUAL-BASED) WEIGHTS — Tukey biweight.
    //
    // recencyWeights alone (age only) turned out insufficient in practice:
    // on a real held-out-forecast test (patient with troughs 4.5, 12.9,
    // 11.1, 7.4, 6.2 over ~3 weeks, next trough actually 6.8), a 45-day
    // half-life barely moved the fit (11.3 vs 11.6 ng/mL unweighted) — the
    // drift happened faster than the age decay discounted it. Shortening
    // the half-life enough to fix that case (~5-7d, floor 0) would make
    // ANY single new level dominate the fit for every patient, including a
    // stable one seen quarterly — the exact over-reaction the floor existed
    // to prevent. A single global age constant cannot serve both.
    //
    // A level should be discounted because it DISAGREES with what the rest
    // of the data implies about current clearance, not because it is
    // chronologically old — a value from a stable long-term patient's last
    // quarterly visit is still fully informative; a value from three weeks
    // ago that no longer matches anything is not, regardless of its age.
    // Tukey's biweight (Mosteller & Tukey; c=4.685 for ~95% efficiency under
    // Gaussian residuals) downweights by residual size, smoothly to zero
    // past the cutoff, and is applied ON TOP of the (now secondary) age
    // weight so recency still breaks genuine ties.
    //
    // Floored at 0.05, not 0 — IRLS re-derives residuals from the CURRENT
    // fit each pass, so a point at 0 weight can never re-enter regardless
    // of what later iterations find; a small floor keeps it recoverable.
    // ============================================================
    const TUKEY_C = 4.685;
    function tukeyWeight(residualRatio) {
        const u = residualRatio / TUKEY_C;
        if (Math.abs(u) >= 1) return 0.05;
        return Math.pow(1 - u * u, 2);
    }

    function combinedWeights(measuredLevels, allDoses, params, bioassay, ageWeights) {
        return normalizeWeights(measuredLevels.map((obs, i) => {
            const pred = predictAtTime(obs.time, allDoses, params);
            const predAdj = bioassay === 2 ? cmiaAdjust(pred) : pred;
            if (!(predAdj > 0.1)) return ageWeights[i] * 0.05;
            const r = (obs.level - predAdj) / (predAdj * PK_MODEL.SIGMA);
            return ageWeights[i] * tukeyWeight(r);
        }));
    }

    // `regime` names which RECENCY_REGIMES entry to weight with. It is a
    // plain argument, never auto-selected here: selectRecencyRegime() calls
    // this function repeatedly to score each regime, so a fit that chose its
    // own regime would recurse. Callers that want the empirically-chosen
    // regime call selectRecencyRegime first and pass the winner through.
    function mapBayesian(popParams, measuredLevels, allDoses, bioassay, regime) {
        const { OMEGA_CL, OMEGA_V, SIGMA } = PK_MODEL;

        if (measuredLevels.length === 0) {
            return {
                ...popParams, etaCL: 0, etaV: 0,
                rmse: null, mpe: null, mape: null, r2: null, pe: [],
                boundaryHit: false
            };
        }

        const regimeUsed = regime || defaultRegimeFor(measuredLevels);
        const ageWeights = recencyWeights(measuredLevels, regimeUsed);

        // ── One MAP optimization pass at a fixed set of weights ────────────
        // OFV = ηCL²/ωCL + ηV²/ωV  +  Σ[wᵢ · (Cobs − Cpred)² / (Cpred·σ)²]
        // Factored out so IRLS (below) can re-run it against updated weights
        // without duplicating the grid+Nelder-Mead search.
        const fitAt = (weights) => {
            const ofv = (eCL, eV) => {
                // SPREAD popParams — do not enumerate {CL, V, KA}. predictAtTime
                // rescales each dose by (d.weight / params.weight) for time-varying
                // body weight, and that gate is silently skipped when `weight` is
                // absent. Enumerating the three fields dropped it here, so the fit
                // saw a constant-weight patient while the plotted population curve
                // saw a varying one: MAP then absorbed the missing correction as
                // slow clearance (ηCL = −0.151, CL 14% low, on data this very model
                // generated). Same reason applies at every other rebuild below.
                const p = {
                    ...popParams,
                    CL: popParams.CL * Math.exp(eCL),
                    V: popParams.V * Math.exp(eV)
                };
                const prior = (eCL * eCL) / OMEGA_CL + (eV * eV) / OMEGA_V;
                const lik = measuredLevels.reduce((s, obs, i) => {
                    const pred = predictAtTime(obs.time, allDoses, p);
                    const predAdj = bioassay === 2 ? cmiaAdjust(pred) : pred;
                    if (predAdj <= 0.1) return s + 1000 * weights[i];
                    const w = predAdj * SIGMA;   // WLS: predicted-value denominator
                    return s + weights[i] * Math.pow((obs.level - predAdj) / w, 2);
                }, 0);
                return prior + lik;
            };

            // ── Phase 1: Coarse global grid — basin identification ───────────
            // Step 0.25 across full eta space to locate the basin of attraction.
            // Kept coarse intentionally; Nelder-Mead refines from here.
            // Driven by an integer step count: accumulating `e += 0.25` in a float
            // loop drifts, so the final grid point could fall just past the bound
            // and be skipped, making the search grid asymmetric.
            const GRID_STEP = 0.25;
            const nCL = Math.round((ETA_CL_BOUNDS[1] - ETA_CL_BOUNDS[0]) / GRID_STEP);
            const nV = Math.round((ETA_V_BOUNDS[1] - ETA_V_BOUNDS[0]) / GRID_STEP);
            let gridCL = 0, gridV = 0, minOFV = Infinity;
            for (let i = 0; i <= nCL; i++) {
                const eCL = ETA_CL_BOUNDS[0] + i * GRID_STEP;
                for (let j = 0; j <= nV; j++) {
                    const eV = ETA_V_BOUNDS[0] + j * GRID_STEP;
                    const v = ofv(eCL, eV);
                    if (v < minOFV) { minOFV = v; gridCL = eCL; gridV = eV; }
                }
            }

            // ── Phase 2: Nelder-Mead simplex refinement ─────────────────────
            // Starts at the grid winner; converges to sub-1e-8 precision,
            // eliminating the ~7% CL error from a 0.15-step grid.
            // nelderMead2D is UNCONSTRAINED, so clamp back into the search box —
            // otherwise indParams could carry an eta the grid would never allow.
            const nm = nelderMead2D(ofv, gridCL, gridV);
            return { bestCL: clamp(nm.x, ETA_CL_BOUNDS), bestV: clamp(nm.y, ETA_V_BOUNDS) };
        };

        // ── IRLS: alternate fitting and robust reweighting ──────────────────
        // 4 passes: converges in practice within 2-3 for a problem this small
        // (few observations, 2 free parameters); the 4th confirms stability
        // rather than assuming it. Each pass re-derives Tukey weights from
        // the residuals of the PREVIOUS pass's fit, so a level's influence
        // reflects how well it agrees with where the data has converged, not
        // a single initial (possibly still-compromised) fit.
        const IRLS_ITERS = 4;
        let weights = ageWeights;
        let bestCL = 0, bestV = 0;
        for (let iter = 0; iter < IRLS_ITERS; iter++) {
            ({ bestCL, bestV } = fitAt(weights));
            const p = {
                ...popParams,
                CL: popParams.CL * Math.exp(bestCL),
                V: popParams.V * Math.exp(bestV)
            };
            weights = combinedWeights(measuredLevels, allDoses, p, bioassay, ageWeights);
        }

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

        // Spread popParams so covariates predictAtTime needs (currently
        // `weight`, for time-varying allometric scaling) reach every consumer of
        // indParams: the plotted individual curve, the dose optimizer, and the
        // accuracy metrics. Dropping it made the individual curve jump 30% the
        // moment the first level was entered — the n=0 branch above already
        // spreads, so the two disagreed.
        const indParams = {
            ...popParams,
            CL: popParams.CL * Math.exp(bestCL),
            V: popParams.V * Math.exp(bestV),
            etaCL: bestCL,
            etaV: bestV,
            outlierHit,
            boundaryHit,
            boundaryMsg: boundaryMsg.trim(),
            // Carried so every downstream consumer (laplacePosterior's CI,
            // the UI) reweights with the SAME regime this fit used, rather
            // than re-deriving and possibly disagreeing with it.
            regime: regimeUsed
        };
        const metrics = calculateAccuracyMetrics(measuredLevels, allDoses, indParams, bioassay);
        return { ...indParams, ...metrics };
    }

    // ============================================================
    // PER-PATIENT REGIME SELECTION
    //
    // Picks the recency regime by measuring which one has actually
    // predicted THIS patient's recent troughs better, instead of assuming
    // one from post-transplant day.
    //
    // Why: a global day-based rule was measurably wrong for some patients.
    // On a 92-prediction held-out backtest over 10 real patients, forcing
    // 'adaptive' inside the maturation window helped genuine drifters a lot
    // (arun_chougle 107%->64% MAPE, raj_bahadur 71%->42%) but hurt
    // vidyashree (28%->37%), whose troughs oscillate without trending. Post-
    // transplant day cannot distinguish those two patients; their own
    // prediction track record can.
    //
    // Method — walk-forward validation on the patient's own history: for
    // each of the last REGIME_SELECT_TRIALS levels, fit on only the levels
    // BEFORE it under each regime, predict it, and score by mean absolute
    // percentage error. Strictly out-of-sample: the target level never
    // enters the fit that predicts it.
    //
    // The selector used to make stationary win ties and required adaptive to
    // clear a 10% margin. That was too conservative on real centre data:
    // forcing adaptive beat the selector overall, because the selector fell
    // back to stationary on drifting patients exactly when adaptation mattered.
    // Now adaptive wins ties/small differences; stationary still wins when it
    // actually scores lower.
    // ============================================================
    const REGIME_SWITCH_MARGIN = 0;

    function selectRecencyRegime(popParams, measuredLevels, allDoses, bioassay) {
        const sorted = [...(measuredLevels || [])].sort((a, b) => a.time - b.time);
        const n = sorted.length;
        const fallback = {
            regime: defaultRegimeFor(sorted),
            basis: 'post-transplant day (too few levels to measure)',
            scores: null, trials: null, nTrials: 0
        };
        if (n < PK_MODEL.MIN_LEVELS_FOR_REGIME_SELECT) return fallback;

        // Each trial needs >=3 prior levels to fit against.
        const firstTarget = Math.max(3, n - PK_MODEL.REGIME_SELECT_TRIALS);
        const names = Object.keys(PK_MODEL.RECENCY_REGIMES);
        const scores = {};
        // Per-trial detail is kept, not just the mean: it IS the patient's
        // prediction track record ("last 3 forecasts for this patient were
        // off by 12%, 8%, 21%"), which the UI shows so a clinician can
        // calibrate trust from observed performance rather than from a
        // methodology note. Computing it here means one walk-forward pass
        // serves both regime selection and the displayed record.
        const trials = {};
        let nTrials = 0;

        for (const name of names) {
            const rows = [];
            for (let i = firstTarget; i < n; i++) {
                const target = sorted[i];
                if (!(target.level > 0)) continue;
                const fit = mapBayesian(popParams, sorted.slice(0, i), allDoses, bioassay, name);
                const pred = generateCurve([target.time], allDoses, fit, bioassay)[0];
                if (!(pred > 0) || !isFinite(pred)) continue;
                rows.push({
                    time: target.time,
                    observed: target.level,
                    predicted: pred,
                    pctError: (pred - target.level) / target.level * 100,
                    absPctError: Math.abs(pred - target.level) / target.level * 100
                });
            }
            trials[name] = rows;
            scores[name] = rows.length
                ? rows.reduce((s, r) => s + r.absPctError, 0) / rows.length / 100
                : null;
            nTrials = Math.max(nTrials, rows.length);
        }

        const stationary = scores.stationary;
        const adaptive = scores.adaptive;
        if (stationary == null || adaptive == null) return fallback;

        const adaptiveWins = adaptive <= stationary * (1 - REGIME_SWITCH_MARGIN);
        return {
            regime: adaptiveWins ? 'adaptive' : 'stationary',
            basis: 'walk-forward validation on this patient',
            scores, trials, nTrials
        };
    }

    // ============================================================
    // CAN THIS FORECAST SETTLE THE CLINICAL QUESTION?
    //
    // A trough number is only actionable if it can distinguish "in the
    // target band" from "outside it". This app's own bands are narrow —
    // 10-11, 7-9, 5-7, 4-6 ng/mL (getTherapeuticRange), i.e. 1-2 ng/mL
    // wide — while a 92-prediction held-out backtest over 10 real patients
    // put typical forecast error at ~26% MAPE. At a predicted 8 ng/mL that
    // is +/-2 ng/mL: wider than the entire 7-9 band. So for many patients
    // the honest answer is "this forecast cannot tell you whether they will
    // be in range — measure".
    //
    // Deliberately driven by the patient's OWN demonstrated error
    // (selectRecencyRegime's walk-forward record) when available, falling
    // back to the population backtest figure otherwise. Both are measured,
    // neither is asserted.
    //
    // This intentionally does NOT emit a "recheck in N days" interval:
    // monitoring frequency is center protocol and varies, and inventing
    // day numbers inside a clinical tool would dress a guess as guidance.
    // It answers only what the model can support — whether the number is
    // strong enough to act on, or whether a level is needed first.
    // ============================================================
    const POPULATION_BACKTEST_MAPE = 26;   // %, 92 held-out predictions / 10 patients

    function forecastResolution(predicted, range, expectedPctError) {
        if (!(predicted > 0) || !range || !(range.low > 0) || !(range.high > 0)) return null;
        const err = (typeof expectedPctError === 'number' && isFinite(expectedPctError) && expectedPctError > 0)
            ? expectedPctError
            : POPULATION_BACKTEST_MAPE;
        const lo = predicted * (1 - err / 100);
        const hi = predicted * (1 + err / 100);
        let status;
        if (lo >= range.low && hi <= range.high) status = 'in';
        else if (lo > range.high) status = 'above';
        else if (hi < range.low) status = 'below';
        else status = 'unresolved';
        return { status, lo, hi, err, resolves: status !== 'unresolved' };
    }

    // ============================================================
    // ACCURACY METRICS
    // RMSE, MPE (bias), MAPE (imprecision), R²
    //
    // UNFITTABLE OBSERVATIONS ARE EXCLUDED. A level drawn before any dose the
    // model knows about predicts ~0, and the MAP objective already ignores it:
    // its `predAdj <= 0.1` branch returns a FLAT +1000 penalty, constant in
    // (ηCL, ηV), so the point cannot move the estimate in any direction.
    //
    // Scoring it anyway was incoherent — the same observation was simultaneously
    // ignored by the fit and counted as a 100% error against it. One trough
    // entered before the dose log starts (routine: the log often begins at
    // admission, the level came from the referring unit) produced RMSE 5.0,
    // which tripped `rmse > 2.5` and put "🚨 CRITICAL: Severe Model Mismatch" on
    // an otherwise perfect fit, with nothing on screen explaining why.
    //
    // Such points are still RETURNED in `pe` with excluded: true so the UI can
    // list them and say what happened — they are just kept out of the summary
    // statistics. `nExcluded` gives the UI its count without re-filtering.
    // ============================================================
    const MIN_FITTABLE_PRED = 0.1;   // matches the OFV's degenerate-prediction guard

    function calculateAccuracyMetrics(measuredLevels, allDoses, params, bioassay) {
        if (measuredLevels.length === 0) return { rmse: null, mpe: null, mape: null, r2: null, pe: [], nExcluded: 0 };
        const all = measuredLevels.map((obs, i) => {
            const pred = (() => {
                const p = predictAtTime(obs.time, allDoses, params);
                return bioassay === 2 ? cmiaAdjust(p) : p;
            })();
            const error = obs.level - pred;
            const pct = obs.level > 0 ? (error / obs.level) * 100 : 0;
            // Excluded on the PREDICTION, not the observation: the question is
            // whether the model had any drug on board to predict with.
            const excluded = !(pred > MIN_FITTABLE_PRED);
            return { n: i + 1, time: obs.time, observed: obs.level, predicted: pred, error, pct, excluded };
        });
        const pe = all;
        const scored = all.filter(p => !p.excluded);
        const nExcluded = all.length - scored.length;
        // Every observation unfittable — there is no fit to describe.
        if (scored.length === 0) return { rmse: null, mpe: null, mape: null, r2: null, pe, nExcluded };
        const n = scored.length;
        const rmse = Math.sqrt(scored.reduce((s, p) => s + p.error * p.error, 0) / n);
        const mpe = scored.reduce((s, p) => s + p.pct, 0) / n;
        const mape = scored.reduce((s, p) => s + Math.abs(p.pct), 0) / n;
        const meanObs = scored.reduce((s, p) => s + p.observed, 0) / n;
        const ssTot = scored.reduce((s, p) => s + Math.pow(p.observed - meanObs, 2), 0);
        const ssRes = scored.reduce((s, p) => s + p.error * p.error, 0);
        const r2 = ssTot > 1e-6 ? Math.max(0, 1 - ssRes / ssTot) : null;
        return { rmse, mpe, mape, r2, pe, nExcluded };
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
            // Spread indParams — see the note in mapBayesian's ofv. Enumerating
            // {CL, V, KA} here dropped the time-varying weight covariate, so the
            // sensitivities J were computed on a different patient than the fit.
            const p = {
                ...indParams,
                CL: indParams.CL * Math.exp(eCL),
                V: indParams.V * Math.exp(eV)
            };
            const c = predictAtTime(t, allDoses, p);
            return bioassay === 2 ? cmiaAdjust(c) : c;
        };

        // Accumulate JᵀJ by central differences on the log scale, weighted
        // by the same combined (recency × robust) weights the fit itself
        // converged on (mapBayesian) — otherwise a discounted level would
        // still count in full toward narrowing the CI, understating
        // uncertainty on a patient whose effective sample size is smaller
        // than raw N. Recomputed from indParams (the converged fit) rather
        // than threaded through as a return value, since it's a pure
        // function of (measuredLevels, allDoses, indParams, bioassay) —
        // the same inputs this function already takes.
        const weights = combinedWeights(measuredLevels, allDoses, indParams, bioassay,
            recencyWeights(measuredLevels, indParams && indParams.regime));
        const h = 1e-5;
        let Sxx = 0, Sxy = 0, Syy = 0;
        measuredLevels.forEach((obs, i) => {
            const c0 = predAt(obs.time, 0, 0);
            if (!(c0 > 0.1)) return;   // uninformative / degenerate point
            const gx = (predAt(obs.time, h, 0) - predAt(obs.time, -h, 0)) / (2 * h) / c0;
            const gy = (predAt(obs.time, 0, h) - predAt(obs.time, 0, -h)) / (2 * h) / c0;
            if (!isFinite(gx) || !isFinite(gy)) return;
            Sxx += weights[i] * gx * gx;
            Sxy += weights[i] * gx * gy;
            Syy += weights[i] * gy * gy;
        });

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

            // Spread indParams — see the note in mapBayesian's ofv. Without it
            // the CI band was drawn for a constant-weight patient while the
            // individual curve it wraps used the time-varying one.
            const p = {
                ...indParams,
                CL: indParams.CL * Math.exp(etaCL),
                V: indParams.V * Math.exp(etaV)
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
    //
    // Site rule: THERE IS NEVER A DRUG-FREE DAY. The dose log is a CHANGE-POINT
    // log — a row is entered when the regimen changes, not for every tablet —
    // so every 12-hourly slot between two entries is reconstructed by carrying
    // the entered regimen forward until a revised dose is entered.
    // ============================================================

    // How far back a same-shift entry may sit and still count as part of the
    // CURRENT regimen. 24 h covers an AM and a PM entry made for the same
    // regimen change, even when they are entered a day apart.
    const REGIMEN_WINDOW_HR = 24;

    // Reconstructed slots must never be mistaken for clinician-entered ones.
    const isImputedDose = d => /^(inter|ext|bridge)-/.test(String(d.id));

    // ── Dose in force at `slotDate` ───────────────────────────────────────
    // The previous rule was "carry the last same-shift POSITIVE dose". That is
    // right when both shifts of a regimen are entered, but wrong when a change
    // is entered for ONE shift only: the other shift then kept a value that
    // could be arbitrarily old.
    //
    // Real case (raj bahadur): 1.5 mg pre-transplant desensitisation doses
    // entered at 04-Jul 07:00 and 19:00, then only the 07:00 slot entered at
    // 4.5 mg on 14-Jul. Every PM slot from 14-Jul to 23-Jul inherited the
    // 10-day-old pre-transplant 1.5 mg, so the engine modelled 4.5/1.5 =
    // 6 mg/day for a patient on 4.5 BID = 9 mg/day. MAP absorbed the missing
    // drug as slow clearance (CL 18.7 vs 25.9 L/h, 28% low), and the modelled
    // daily dose then jumped 33% the moment both shifts were entered again on
    // 24-Jul — pushing the 28-Jul forecast to 15.2 ng/mL against an observed
    // 7.4.
    //
    // New rule: a same-shift entry is carried only while it is still part of
    // the current regimen — i.e. not more than REGIMEN_WINDOW_HR older than the
    // most recent entry on ANY shift. Otherwise the most recent entry wins,
    // whichever shift it came from. A deliberate asymmetric split (5 mg AM /
    // 4.5 mg PM entered together) is preserved; a single revised entry replaces
    // a stale regimen outright.
    //
    // Only clinician-ENTERED doses may define the regimen. Letting imputed
    // slots define it re-seeds the stale value at every step, so the staleness
    // test could never fire.
    function regimenDoseAt(entered, slotDate, fallback) {
        // 0 mg is a deliberate "dose held" marker for that one slot. It never
        // defines the ongoing regimen, so it is excluded from both candidates.
        const prior = entered
            .filter(d => d.dose > 0 && d.recordDate.isBefore(slotDate))
            .sort((a, b) => a.recordDate.valueOf() - b.recordDate.valueOf());
        if (prior.length === 0) return fallback;

        const lastAny = prior[prior.length - 1];
        const isAm = slotDate.hour() < 12;
        const lastSame = [...prior].reverse()
            .find(d => (d.recordDate.hour() < 12) === isAm);

        const stillCurrent = lastSame && !lastSame.recordDate
            .isBefore(lastAny.recordDate.subtract(REGIMEN_WINDOW_HR, 'hour'));
        return stillCurrent ? lastSame.dose : lastAny.dose;
    }

    function fillHistoricalGaps(historyLog, txDate) {
        const filled = [];
        const sortedHistory = [...historyLog].sort((a, b) => a.time - b.time);

        // Time-varying covariates (body weight, haematocrit) are recorded on
        // whatever event the clinician happened to enter them on, and apply from
        // that moment until the next recorded value — the same carry-forward rule
        // the dose regimen uses. Returns the values in force at `t`.
        //
        // NOTE on the weight model: predictAtTime scales each dose by the weight
        // in force when that dose was GIVEN, and keeps it for that dose's whole
        // disposition. A dose given at 60 kg is therefore still eliminated at
        // 60 kg clearance ten days later at 45 kg. For slowly drifting weight
        // that is a fair approximation; for rapidly resolving post-operative
        // oedema — the main reason this feature exists — it is the weakest link.
        // A fully correct treatment needs time-varying ke integrated along the
        // elimination phase, not superposition of fixed-parameter doses.
        function getCovariatesAtTime(t) {
            let wt = null, hct = null, inh = null;
            for (const item of sortedHistory) {
                if (item.time <= t) {
                    if (item.weight != null) wt = item.weight;
                    if (item.hematocrit != null) hct = item.hematocrit;
                    if (item.inhibitor != null && item.inhibitor !== '') inh = item.inhibitor;
                } else {
                    break;
                }
            }
            return { weight: wt, hematocrit: hct, inhibitor: inh };
        }

        // Anchor on dose >= 0, NOT dose > 0. A recorded 0 mg is a deliberate
        // "dose held" entry. Filtering it out dropped it from the anchor list,
        // and the 12-hourly gap-filler then spanned straight across it and
        // re-inserted a full dose in that slot — so recording a held dose had
        // the exact opposite of the intended effect.
        const logged = historyLog.filter(e => e.dose !== null && e.dose >= 0).sort((a, b) => a.time - b.time);
        if (logged.length === 0) return [];

        for (let i = 0; i < logged.length; i++) {
            const cov = getCovariatesAtTime(logged[i].time);
            filled.push({
                ...logged[i],
                weight: logged[i].weight ?? cov.weight,
                hematocrit: logged[i].hematocrit ?? cov.hematocrit,
                inhibitor: logged[i].inhibitor ?? cov.inhibitor
            });
            if (i === logged.length - 1) break;

            let nextTime = logged[i].recordDate.add(12, 'hour');
            // Fill 12-hour steps until we hit the next logged dose (with a 4hr buffer to avoid overlapping)
            while (nextTime.isBefore(logged[i + 1].recordDate.subtract(4, 'hour'))) {
                // Carry the regimen in force at this slot. Reads `logged`, not
                // `filled`: only entered doses define the regimen.
                const doseToGive = regimenDoseAt(logged, nextTime, logged[i].dose);
                const slotTime = nextTime.diff(txDate, 'hour', true);
                const slotCov = getCovariatesAtTime(slotTime);
                filled.push({
                    id: 'inter-' + i + '-' + nextTime.valueOf(),
                    recordDate: nextTime,
                    dose: doseToGive,
                    level: null,
                    weight: slotCov.weight,
                    hematocrit: slotCov.hematocrit,
                    inhibitor: slotCov.inhibitor,
                    time: slotTime
                });
                nextTime = nextTime.add(12, 'hour');
            }
        }
        return filled;
    }

    function extrapolateDoses(filledHistory, untilDate, txDate) {
        if (filledHistory.length === 0) return [];
        const last = filledHistory[filledHistory.length - 1];
        // Future slots inherit the most recent known covariates — there is no
        // later measurement to carry forward from.
        const lastWt = last ? last.weight : null;
        const lastHct = last ? last.hematocrit : null;
        const lastInh = last ? last.inhibitor : null;
        // Same regimen rule as fillHistoricalGaps, driven off the entered doses
        // only — filledHistory already carries imputed slots by this point.
        const entered = filledHistory.filter(d => !isImputedDose(d));
        const source = entered.length ? entered : filledHistory;
        const extra = [];
        let nextTime = last.recordDate.add(12, 'hour');
        while (nextTime.isBefore(untilDate)) {
            const doseToGive = regimenDoseAt(source, nextTime, last.dose);
            extra.push({
                id: 'ext-' + extra.length,
                recordDate: nextTime,
                dose: doseToGive,
                level: null,
                weight: lastWt,
                hematocrit: lastHct,
                inhibitor: lastInh,
                time: nextTime.diff(txDate, 'hour', true)
            });
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
    // Same regimen rule as fillHistoricalGaps (see regimenDoseAt), so a held
    // dose does not propagate a zero regimen across the whole bridge.
    // ============================================================
    function buildBridgeDoses(logged, fromDate, toDate, txDate) {
        const bridge = [];
        if (!logged || logged.length === 0) return bridge;
        const fallback = logged[logged.length - 1].dose;
        const last = logged[logged.length - 1];
        const lastWt = last ? last.weight : null;
        const lastHct = last ? last.hematocrit : null;
        const lastInh = last ? last.inhibitor : null;
        // Callers pass the forecast's dose array, which mixes entered and
        // imputed slots; only the entered ones may define the regimen.
        const entered = logged.filter(d => !isImputedDose(d));
        const source = entered.length ? entered : logged;
        for (let t = fromDate; t.isBefore(toDate); t = t.add(12, 'hour')) {
            bridge.push({
                id: 'bridge-' + bridge.length,
                recordDate: t,
                dose: regimenDoseAt(source, t, fallback),
                level: null,
                weight: lastWt,
                hematocrit: lastHct,
                inhibitor: lastInh,
                time: t.diff(txDate, 'hour', true)
            });
        }
        return bridge;
    }

    // ============================================================
    // C/D RATIO SUPPORT
    //
    // The concentration/dose ratio needs the drug the patient ACTUALLY took in
    // the 24 h before the sample. For a change-point log that means the
    // RECONSTRUCTED regimen, not the rows that happen to be entered.
    //
    // Using raw entries made the C/D column and the IPV panel see a drug-free
    // patient — with only regimen-change rows entered, every 24 h window came
    // back empty, so the column printed '—' and IPV reported "insufficient dose
    // data" — while the MAP forecast, which already gap-fills, saw a continuous
    // BID regimen. Two halves of the app disagreeing about the dose history.
    // ============================================================

    // Logged doses + the 12-hourly doses implied between them, extended past
    // the last entered dose so a level drawn days later still has its
    // preceding 24 h covered. There is never a drug-free day.
    function buildEffectiveDoses(historyLog, txDate) {
        const filled = fillHistoricalGaps(historyLog, txDate);
        if (filled.length === 0) return [];

        const lastEvent = historyLog.reduce((m, e) => Math.max(m, e.time), -Infinity);
        const lastDose = filled[filled.length - 1];
        if (!(lastEvent > lastDose.time)) return filled;

        // extrapolateDoses stops strictly BEFORE untilDate, so a trough drawn
        // exactly on a dosing slot does not get that slot counted against it —
        // which is correct: a pre-dose sample precedes the dose.
        return filled.concat(
            extrapolateDoses(filled, txDate.add(lastEvent, 'hour'), txDate));
    }

    // Total mg given in the 24 h ending at `timeHr` (hours since transplant).
    // Inclusive at exactly −24 h, exclusive at timeHr itself, matching the
    // pre-dose sampling convention.
    function dailyDoseBefore(doses, timeHr) {
        return doses
            .filter(d => d.dose > 0 && d.time < timeHr && d.time >= timeHr - 24)
            .reduce((s, d) => s + d.dose, 0);
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

    // ============================================================
    // C/D RATIO VARIABILITY (IPV) — coefficient of variation of the
    // concentration/dose ratio across measured levels, using the
    // reconstructed (imputed) daily dose. Sapir-Pichhadze 2014: IPV >40% CV
    // is associated with increased rejection risk; the UI's IPV panel has
    // used this threshold since it was written.
    //
    // Pulled out as a pure function so the forecast panel can reuse the exact
    // same number the IPV panel shows, rather than a second inline
    // computation drifting from it. mapBayesian has no recency weighting —
    // every observation counts equally regardless of age — so a high IPV
    // (the C/D ratio genuinely swinging, not just trending) is precisely the
    // condition under which a single time-invariant MAP fit struggles, and
    // is worth flagging on the forecast itself, not only in the IPV tile.
    // ============================================================
    function calculateIPV(measuredLevels, allDoses) {
        const pairs = (measuredLevels || [])
            .map(lev => {
                const daily = dailyDoseBefore(allDoses, lev.time);
                return daily > 0 ? { time: lev.time, cd: lev.level / daily } : null;
            })
            .filter(Boolean);
        if (pairs.length < 2) return { ipv: null, mean: null, median: null, n: pairs.length, pairs };
        const mean = pairs.reduce((s, p) => s + p.cd, 0) / pairs.length;
        const sd = Math.sqrt(pairs.reduce((s, p) => s + Math.pow(p.cd - mean, 2), 0) / (pairs.length - 1));
        const ipv = mean > 0 ? (sd / mean) * 100 : null;
        // ── Two different central-tendency measures, deliberately ────────────
        // `mean` is arithmetic because IPV is defined as SD/mean × 100 — CV%
        // has no meaning on any other centre, so the CV keeps it.
        //
        // `median` exists for the metabolizer LABEL, which is a different job.
        // C/D is ratio data and its distribution is right-skewed; on a patient
        // whose C/D spanned 0.21–2.37 (11-fold) the arithmetic mean was dragged
        // to 0.96 by two high outliers, landing just above the fast-metabolizer
        // cutoff and printing "Normal Metabolizer" for a patient who could not
        // hold a therapeutic level on 9 mg/day and was switched off tacrolimus.
        // The median is unmoved by those two points and classifies correctly.
        const sortedCd = pairs.map(p => p.cd).sort((a, b) => a - b);
        const mid = Math.floor(sortedCd.length / 2);
        const median = sortedCd.length % 2 === 0
            ? (sortedCd[mid - 1] + sortedCd[mid]) / 2
            : sortedCd[mid];
        return { ipv, mean, median, n: pairs.length, pairs };
    }

    // ── Metabolizer classification from the C/D ratio ────────────────────────
    // Thresholds are Thölking 2014 (PLoS One 9:e111128), as replicated in
    // Thölking 2016 and Schütte-Nütgen 2019: fast <1.05, intermediate 1.05–2.0,
    // slow >2.0 ng/mL per mg/24h. This app previously used 0.9 / 1.5, which are
    // not the published values and shifted the fast/intermediate boundary far
    // enough that a genuine fast metabolizer at C/D 0.96 was reported "Normal".
    //
    // Takes the median C/D (see calculateIPV), not the arithmetic mean.
    function classifyMetabolizer(cdRatio) {
        if (typeof cdRatio !== 'number' || !isFinite(cdRatio) || cdRatio <= 0) return null;
        if (cdRatio < 1.05) return 'Fast Metabolizer';
        if (cdRatio <= 2.0) return 'Intermediate Metabolizer';
        return 'Slow Metabolizer';
    }

    return {
        // constants
        PK_MODEL, DOSE_SLOTS, SAMPLING_LEAD_MIN,
        ETA_CL_BOUNDS, ETA_V_BOUNDS, ETA_BOUNDARY_TOL,
        // time / protocol helpers
        nextC0Time, getTherapeuticRange, getLevelStatus,
        // PK model
        getPopulationParameters, predictSingleDose, predictAtTime,
        cmiaAdjust, generateCurve, suggestStartingDose,
        // estimation
        nelderMead2D, mapBayesian, calculateAccuracyMetrics,
        laplacePosterior, monteCarloCI, recencyWeights, combinedWeights,
        selectRecencyRegime, defaultRegimeFor,
        // dose series
        fillHistoricalGaps, extrapolateDoses, buildBridgeDoses,
        regimenDoseAt, buildEffectiveDoses, dailyDoseBefore, REGIMEN_WINDOW_HR,
        // analytics / QC
        rosendaalFraction, findPostDoseSamples, calculateIPV, classifyMetabolizer,
        forecastResolution, POPULATION_BACKTEST_MAPE
    };
}));
