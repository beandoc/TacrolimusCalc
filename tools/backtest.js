// ============================================================
// WALK-FORWARD BACKTEST HARNESS
//
// Why this exists: POPULATION_BACKTEST_MAPE = 26 (pk-engine.js) and every
// accuracy claim in the engine's comments ("arun_chougle 107%->64%",
// "vidyashree 28%->37%", "92 held-out predictions / 10 patients") were
// produced by a script that was never committed. They are unreproducible
// and unauditable, yet one of them is rendered to clinicians as a measured
// figure. This makes them reproducible.
//
// It also exists because the numbers being compared are far smaller than
// the noise on them. On this cohort the 95% CI on an absolute MAPE is
// roughly +/-5 points. A 0.5-point difference between two configurations
// is not a result. So this harness NEVER reports two independent MAPEs and
// invites you to subtract them: it reports the PAIRED difference and a
// bootstrap CI on that difference, because the same 92 predictions under
// two configurations are highly correlated and the paired CI is several
// times tighter than the difference of two independent ones.
//
// Usage:
//   node tools/backtest.js                 # all configs, default DB
//   node tools/backtest.js --db path.db
//   node tools/backtest.js --doses both    # also run the truncated-dose arm
//
// Requires the sqlite3 CLI (ships with macOS; `apt install sqlite3` on Linux)
// and a local tacrolimus.db. Both the DB and its patient data are gitignored
// PHI — this script reads them, and prints only aggregate numbers plus MRNs
// that already appear in pk-engine.js comments.
// ============================================================

const { execFileSync } = require('child_process');
const fs = require('fs');
const path = require('path');
const dayjs = require('dayjs');
const E = require('../pk-engine');

// ── args ────────────────────────────────────────────────────────────────────
const argv = process.argv.slice(2);
const argOf = (flag, dflt) => {
    const i = argv.indexOf(flag);
    return i >= 0 && argv[i + 1] ? argv[i + 1] : dflt;
};
const DB = path.resolve(argOf('--db', path.join(__dirname, '..', 'tacrolimus.db')));
const DOSE_MODE = argOf('--doses', 'full');       // full | truncated | both
const BOOT = parseInt(argOf('--boot', '4000'), 10);

// Warm-up must match selectRecencyRegime's `Math.max(3, ...)`: a fit needs
// >=3 prior levels before a held-out target means anything. Exposed as a flag
// ONLY so an older run can be reproduced and compared — raising it silently
// drops the hardest early predictions and flatters every arm, so the value in
// force is printed in the header rather than left implicit.
const WARMUP = parseInt(argOf('--warmup', '3'), 10);

// Test rows in the local DB. Excluded by name because they are fixtures, not
// patients — `raj_bahadur_test` is kept (6 real levels, it is the golden case
// in test/pk-engine.test.js) while the 2-event stub `raj bahadur` is dropped.
const EXCLUDE = new Set(['test_patient', 'test001', 'test_patient_mrn', 'raj bahadur', 'zz_desens_test_patient']);

// Dropping a patient from an accuracy cohort is a claim, not a convenience:
// whichever arm that patient was worst for gets quietly better. Exposed so the
// effect can be MEASURED, and every excluded MRN is echoed in the header.
const DROP = new Set((argOf('--drop', '') || '').split(',').filter(Boolean));
DROP.forEach(m => EXCLUDE.add(m));

// ── data ────────────────────────────────────────────────────────────────────
function q(sql) {
    const out = execFileSync('sqlite3', [DB, '-json', sql], { encoding: 'utf8' }).trim();
    return out ? JSON.parse(out) : [];
}

function loadCohort() {
    if (!fs.existsSync(DB)) {
        console.error(`\n  No database at ${DB}\n  Pass --db <path>, or run against a center export.\n`);
        process.exit(2);
    }
    const patients = q('SELECT * FROM patients');
    const cohort = [];
    for (const p of patients) {
        if (EXCLUDE.has(p.mrn)) continue;
        // Discontinued patients are excluded from calibration cohorts by the
        // same selection-bias argument documented in backend/models.py:28-42.
        if (p.tac_discontinued === 1) continue;

        const rows = q(`SELECT datetime, dose, level, weight, hematocrit, inhibitor
                        FROM clinical_events WHERE patient_mrn = '${p.mrn.replace(/'/g, "''")}'
                        ORDER BY datetime`);
        const txDate = dayjs(p.transplant_date);
        // Mirrors gatherData() in index.html:3082-3096 exactly.
        const historyLog = rows.map((r, i) => {
            const rd = dayjs(r.datetime);
            return {
                id: i, recordDate: rd, dose: r.dose, level: r.level,
                weight: r.weight != null ? r.weight : null,
                hematocrit: r.hematocrit != null ? r.hematocrit : null,
                inhibitor: r.inhibitor != null ? r.inhibitor : null,
                time: rd.diff(txDate, 'hour', true)
            };
        });
        const levels = historyLog.filter(e => e.level !== null && e.level > 0);
        if (levels.length < WARMUP + 1) continue;
        if (!historyLog.some(e => e.dose !== null && e.dose > 0)) continue;

        cohort.push({
            mrn: p.mrn, txDate, historyLog, levels,
            patientData: {
                weight: p.weight, hct: p.hematocrit, mpa: String(p.mpa),
                bioassay: p.bioassay, genotype: p.genotype, albumin: p.albumin,
                bilirubin: p.bilirubin, inhibitor: p.inhibitor, transplantDate: txDate
            }
        });
    }
    return cohort.sort((a, b) => b.levels.length - a.levels.length);
}

// ── dose series ─────────────────────────────────────────────────────────────
// `full` reproduces what the engine does today: the dose array is built from
// the patient's COMPLETE history, including doses recorded after the target.
// predictAtTime filters `d.time < t` so no future dose contributes to the
// predicted value directly — but fillHistoricalGaps interpolates across gaps,
// so a dose recorded AFTER the target can still shape the inferred dose level
// BEFORE it. `truncated` cuts the history at the target first, which is what a
// clinician standing at that moment actually had. Running both measures
// whether that subtlety is worth anything.
function buildDoses(patient, cutoffTime) {
    const log = cutoffTime == null
        ? patient.historyLog
        : patient.historyLog.filter(e => e.time < cutoffTime);
    const filled = E.fillHistoricalGaps(log, patient.txDate);
    const until = patient.txDate.add(
        Math.max(...patient.historyLog.map(e => e.time)) / 24 + 14, 'day');
    return [...filled, ...E.extrapolateDoses(filled, until, patient.txDate)];
}

// ── configurations under test ───────────────────────────────────────────────
// Each returns a predicted trough for `target` given only `priorLevels`.
const CONFIGS = {
    stationary: (pop, priorLevels, doses, bioassay, target) =>
        predictWith(pop, priorLevels, doses, bioassay, target, 'stationary'),

    adaptive: (pop, priorLevels, doses, bioassay, target) =>
        predictWith(pop, priorLevels, doses, bioassay, target, 'adaptive'),

    // The shipped behaviour: the selector runs its own nested walk-forward,
    // and it must only ever see the levels available at this point in time.
    // Passing the full level list here would be exactly the leak this harness
    // is meant to catch.
    selector: (pop, priorLevels, doses, bioassay, target) => {
        const choice = E.selectRecencyRegime(pop, priorLevels, doses, bioassay);
        return predictWith(pop, priorLevels, doses, bioassay, target, choice.regime);
    },

    // Reference floors. If a model cannot clear these comfortably, its
    // sophistication is not buying anything.
    'baseline:carry-last': (pop, priorLevels, doses, bioassay, target) => {
        const last = priorLevels[priorLevels.length - 1];
        const dNow = E.dailyDoseBefore(doses, target.time);
        const dThen = E.dailyDoseBefore(doses, last.time);
        return (dThen > 0 && dNow > 0) ? last.level * (dNow / dThen) : last.level;
    },
    'baseline:constant-8': () => 8.0
};

function predictWith(pop, priorLevels, doses, bioassay, target, regime) {
    const fit = E.mapBayesian(pop, priorLevels, doses, bioassay, regime);
    return E.generateCurve([target.time], doses, fit, bioassay)[0];
}

// ── the walk-forward run ────────────────────────────────────────────────────
function run(cohort, configNames, doseMode) {
    // preds[config] is an array aligned across configs: same patient, same
    // target, same index. That alignment is what makes the pairing valid.
    const preds = {}; configNames.forEach(c => { preds[c] = []; });
    const keys = [];
    let leakChecks = 0;

    for (const p of cohort) {
        const pop = E.getPopulationParameters(p.patientData);
        const bioassay = p.patientData.bioassay || 1;

        for (let i = WARMUP; i < p.levels.length; i++) {
            const target = p.levels[i];
            const priorLevels = p.levels.slice(0, i);

            // Forward-only guard. This is the assertion that would have caught
            // the class of bug this harness was written to rule out.
            for (const l of priorLevels) {
                if (!(l.time < target.time)) {
                    throw new Error(`LEAK: ${p.mrn} fit input at t=${l.time} >= target t=${target.time}`);
                }
                leakChecks++;
            }

            const doses = buildDoses(p, doseMode === 'truncated' ? target.time : null);
            keys.push({ mrn: p.mrn, time: target.time, observed: target.level });

            for (const c of configNames) {
                let v;
                try { v = CONFIGS[c](pop, priorLevels, doses, bioassay, target); }
                catch { v = NaN; }
                preds[c].push(v);
            }
        }
    }
    return { preds, keys, leakChecks };
}

// ── statistics ──────────────────────────────────────────────────────────────
const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
const median = a => { const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const ape = (p, o) => Math.abs(p - o) / o * 100;

// Deterministic RNG so a reported CI is reproducible run to run. A CI that
// moves when you rerun it is one more number nobody can check.
function mulberry32(a) {
    return function () {
        a |= 0; a = a + 0x6D2B79F5 | 0;
        let t = Math.imul(a ^ a >>> 15, 1 | a);
        t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
        return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
}

function bootstrapCI(values, n = BOOT, seed = 12345) {
    const rnd = mulberry32(seed);
    const means = [];
    for (let b = 0; b < n; b++) {
        let s = 0;
        for (let i = 0; i < values.length; i++) s += values[(rnd() * values.length) | 0];
        means.push(s / values.length);
    }
    means.sort((a, b) => a - b);
    return [means[Math.floor(n * 0.025)], means[Math.floor(n * 0.975)]];
}

// Clinical decision buckets. A model that improves MAPE without moving these
// has not changed a single dose decision.
const bucket = v => v < 6 ? 'LOW' : (v <= 10 ? 'IN' : 'HIGH');

function scoreConfig(preds, keys) {
    const errs = [], rows = [];
    for (let i = 0; i < keys.length; i++) {
        const p = preds[i], o = keys[i].observed;
        if (!(p > 0) || !isFinite(p)) continue;
        errs.push(ape(p, o));
        rows.push({ ...keys[i], predicted: p });
    }
    const missedHigh = rows.filter(r => bucket(r.predicted) !== 'HIGH' && bucket(r.observed) === 'HIGH').length;
    const missedLow = rows.filter(r => bucket(r.predicted) !== 'LOW' && bucket(r.observed) === 'LOW').length;
    const hit = rows.filter(r => bucket(r.predicted) === bucket(r.observed)).length;
    return {
        n: errs.length, mape: mean(errs), medianApe: median(errs), errs, rows,
        bucketAcc: hit / rows.length * 100, missedHigh, missedLow,
        ci: bootstrapCI(errs)
    };
}

// ── reporting ───────────────────────────────────────────────────────────────
const pad = (s, n) => String(s).padEnd(n);
const num = (v, w = 5, d = 1) => v.toFixed(d).padStart(w);

function report(cohort, doseMode) {
    const names = Object.keys(CONFIGS);
    const { preds, keys, leakChecks } = run(cohort, names, doseMode);
    const scored = {};
    names.forEach(c => { scored[c] = scoreConfig(preds[c], keys); });

    console.log(`\n${'='.repeat(78)}`);
    console.log(`WALK-FORWARD BACKTEST   dose series: ${doseMode}   warm-up: ${WARMUP} levels`);
    console.log(`${'='.repeat(78)}`);
    console.log(`cohort: ${cohort.length} patients, ${keys.length} held-out predictions`);
    console.log(`forward-only assertions passed: ${leakChecks}`);
    if (DROP.size) console.log(`DROPPED BY FLAG: ${[...DROP].join(', ')}  <- these patients are not in any number below`);
    console.log('');

    console.log(`  ${pad('config', 22)} ${pad('n', 4)} ${pad('MAPE', 7)} ${pad('95% CI', 16)} ${pad('median', 7)} ${pad('bucket', 8)} miss-HIGH miss-LOW`);
    console.log(`  ${'-'.repeat(88)}`);
    for (const c of names) {
        const s = scored[c];
        console.log(`  ${pad(c, 22)} ${pad(s.n, 4)} ${num(s.mape)}%  [${num(s.ci[0], 4)},${num(s.ci[1], 5)}]  ${num(s.medianApe)}%  ${num(s.bucketAcc)}%  ${String(s.missedHigh).padStart(7)}  ${String(s.missedLow).padStart(7)}`);
    }

    // ── the part that actually decides anything ─────────────────────────────
    console.log(`\n  PAIRED DIFFERENCES  (same predictions under both configs; this is the`);
    console.log(`  only comparison that means anything at this sample size)\n`);
    const pairs = [
        ['stationary', 'adaptive'],
        ['selector', 'adaptive'],
        ['selector', 'stationary'],
        ['adaptive', 'baseline:carry-last'],
        ['adaptive', 'baseline:constant-8']
    ];
    console.log(`  ${pad('comparison', 40)} ${pad('diff', 8)} ${pad('95% CI', 18)} verdict`);
    console.log(`  ${'-'.repeat(92)}`);
    for (const [a, b] of pairs) {
        const d = [];
        for (let i = 0; i < keys.length; i++) {
            const pa = preds[a][i], pb = preds[b][i], o = keys[i].observed;
            if (!(pa > 0) || !(pb > 0) || !isFinite(pa) || !isFinite(pb)) continue;
            d.push(ape(pa, o) - ape(pb, o));
        }
        const [lo, hi] = bootstrapCI(d);
        const sig = lo * hi > 0;
        const better = mean(d) < 0 ? a : b;
        console.log(`  ${pad(`${a} - ${b}`, 40)} ${num(mean(d), 6)}pt  [${num(lo, 5)},${num(hi, 6)}]  ${sig ? `${better} genuinely better` : 'NOT distinguishable from 0'}`);
    }

    // ── per-patient: one patient can carry the whole mean ───────────────────
    console.log(`\n  PER-PATIENT MAPE (a 10-patient mean is fragile; look here before believing it)\n`);
    console.log(`  ${pad('patient', 26)} ${pad('n', 4)} ${names.map(c => pad(c.replace('baseline:', 'b:'), 14)).join('')}`);
    console.log(`  ${'-'.repeat(26 + 4 + names.length * 14)}`);
    for (const p of cohort) {
        const idx = keys.map((k, i) => [k, i]).filter(([k]) => k.mrn === p.mrn).map(([, i]) => i);
        if (!idx.length) continue;
        const cells = names.map(c => {
            const e = idx.map(i => [preds[c][i], keys[i].observed])
                .filter(([v]) => v > 0 && isFinite(v)).map(([v, o]) => ape(v, o));
            return pad(e.length ? `${mean(e).toFixed(1)}%` : '-', 14);
        });
        console.log(`  ${pad(p.mrn, 26)} ${pad(idx.length, 4)} ${cells.join('')}`);
    }

    return scored;
}

// ── main ────────────────────────────────────────────────────────────────────
const cohort = loadCohort();
if (!cohort.length) { console.error('  No usable patients.'); process.exit(2); }

const modes = DOSE_MODE === 'both' ? ['full', 'truncated'] : [DOSE_MODE];
const results = {};
for (const m of modes) results[m] = report(cohort, m);

if (modes.length === 2) {
    console.log(`\n${'='.repeat(78)}`);
    console.log('DOSE-SERIES LEAKAGE CHECK');
    console.log(`${'='.repeat(78)}`);
    console.log('  `full` lets fillHistoricalGaps interpolate a pre-target dose using a dose');
    console.log('  recorded after the target. `truncated` cannot. A large gap here means the');
    console.log('  reported accuracy depends on information the clinician did not have.\n');
    for (const c of Object.keys(results.full)) {
        const f = results.full[c].mape, t = results.truncated[c].mape;
        console.log(`  ${pad(c, 22)} full ${num(f)}%   truncated ${num(t)}%   delta ${num(t - f, 6)}pt`);
    }
}

console.log(`\n${'='.repeat(78)}`);
console.log('HOW TO READ THIS');
console.log(`${'='.repeat(78)}`);
console.log('  The absolute MAPE column has a CI roughly +/-5 points wide at this cohort');
console.log('  size. Two configs differing by under ~2 points in that column are not');
console.log('  distinguishable — read the PAIRED table instead, and only act on a row');
console.log('  whose CI excludes zero.');
console.log('');
console.log('  These hyperparameters were all selected against this same cohort:');
console.log('    RECENCY_REGIMES.adaptive (5d/0)   pk-engine.js:109-120');
console.log('    REGIME_SELECT_TRIALS (5)          pk-engine.js:131-135');
console.log('    REGIME_SWITCH_MARGIN (0)          pk-engine.js:834-841');
console.log('    MATURATION_WINDOW_DAYS (90)       pk-engine.js:96-127');
console.log('  So every number above is a resubstitution estimate of a tuned pipeline,');
console.log('  not held-out performance. Treat it as an optimistic bound.\n');
