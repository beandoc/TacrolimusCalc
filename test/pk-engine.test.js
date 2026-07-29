// Unit tests for the PK / MAP-Bayesian core.
//
// These exist because this maths is easy to get plausibly, invisibly wrong.
// Two real examples this suite now guards against:
//   * a proposed "Laplace posterior" reporting sd(etaV) = 0.52 — impossible,
//     since a posterior can never be wider than its 0.300 prior;
//   * an outlier threshold of 1.26 described as 3 SD, when 3 SD is 1.154.
// Both look reasonable in a diff. Neither survives an assertion.

const { E, dayjs, ok, near, section, done, TX, POP, steadyStateDoses, c0Hour, referencePosterior } = require('./helpers');

const doses = steadyStateDoses(40);

section('PK core');
{
    const c = E.predictAtTime(c0Hour(20), doses, POP);
    ok('steady-state C0 on 3 mg BID is clinically plausible', c > 8 && c < 11, `${c.toFixed(2)} ng/mL`);

    const ke = POP.CL / POP.V;
    ok('elimination half-life is in the reported 12-36 h range',
        Math.log(2) / ke > 12 && Math.log(2) / ke < 36, `t1/2 = ${(Math.log(2) / ke).toFixed(1)} h`);

    // The KA ~= ke branch is a removable singularity; it must join smoothly.
    const atSingularity = { CL: 100, V: 100 / 4.5299, KA: 4.53 };
    const justOutside = { CL: 100, V: 100 / 4.52, KA: 4.53 };
    ok('KA=ke L\'Hopital branch is continuous',
        near(E.predictSingleDose(3, 5, atSingularity), E.predictSingleDose(3, 5, justOutside), 0.05));

    ok('no dose contributes before it is given', E.predictSingleDose(-1, 5, POP) === 0);
    ok('a zero (held) dose contributes nothing', E.predictSingleDose(6, 0, POP) === 0);

    // CMIA = 1.08 x LC-MS/MS + 0.55, but the intercept is meaningless at c=0.
    ok('cmiaAdjust adds no intercept at c=0', E.cmiaAdjust(0) === 0);
    ok('cmiaAdjust applies the affine correction at c>0', near(E.cmiaAdjust(10), 11.35, 1e-9));
}

section('MAP Bayesian recovery (noise-free)');
for (const [trueCL, trueV] of [[0, 0], [0.4, -0.3], [-0.5, 0.4], [0.8, 0.2]]) {
    const truth = { CL: POP.CL * Math.exp(trueCL), V: POP.V * Math.exp(trueV), KA: POP.KA };
    const obs = [10, 14, 18, 21].map(d => ({ time: c0Hour(d), level: E.predictAtTime(c0Hour(d), doses, truth) }));
    const fit = E.mapBayesian(POP, obs, doses, 1);
    ok(`recovers etaCL = ${trueCL}`, near(fit.etaCL, trueCL, 0.12),
        `MAP etaCL=${fit.etaCL.toFixed(3)} etaV=${fit.etaV.toFixed(3)}`);
}
ok('n=0 returns the population estimate untouched',
    (() => { const r = E.mapBayesian(POP, [], doses, 1); return r.etaCL === 0 && r.etaV === 0 && r.rmse === null; })());

section('Laplace posterior covariance');
{
    const priorSdCL = Math.sqrt(E.PK_MODEL.OMEGA_CL), priorSdV = Math.sqrt(E.PK_MODEL.OMEGA_V);
    console.log(`  prior sd(etaCL)=${priorSdCL.toFixed(3)}  sd(etaV)=${priorSdV.toFixed(3)}`);

    for (const n of [1, 2, 3, 5]) {
        const obs = Array.from({ length: n }, (_, k) =>
            ({ time: c0Hour(20 + k), level: E.predictAtTime(c0Hour(20 + k), doses, POP) }));
        const got = E.laplacePosterior(obs, doses, POP, 1);
        const exp = referencePosterior(obs.map(o => o.time), doses, POP);

        ok(`n=${n} matches the exact (Omega^-1 + J'J/sigma^2)^-1`,
            near(got.sdCL, exp.sdCL, 1e-4) && near(got.sdV, exp.sdV, 1e-4) && near(got.rho, exp.rho, 1e-3),
            `sd(etaCL)=${got.sdCL.toFixed(3)} sd(etaV)=${got.sdV.toFixed(3)} rho=+${got.rho.toFixed(3)}`);

        // The invariant that would have caught the claimed sd(etaV) = 0.52.
        ok(`n=${n} posterior is never wider than the prior`,
            got.sdCL <= priorSdCL + 1e-12 && got.sdV <= priorSdV + 1e-12,
            `sd(etaV)=${got.sdV.toFixed(3)} <= ${priorSdV.toFixed(3)}`);

        // dlnC/dlnCL and dlnC/dlnV have OPPOSITE signs at a trough, so CL and V
        // are positively correlated along the ridge. An earlier revision
        // hardcoded rho = -0.45, which has the wrong sign.
        ok(`n=${n} CL-V correlation is positive`, got.rho > 0, `rho=+${got.rho.toFixed(3)}`);

        // Troughs carry almost no information about V (dlnC/dlnV ~ +0.20).
        ok(`n=${n} trough data barely shrinks sd(etaV)`, got.sdV > 0.9 * priorSdV,
            `${got.sdV.toFixed(3)} vs prior ${priorSdV.toFixed(3)}`);
    }

    const rhos = [1, 3, 5].map(n => E.laplacePosterior(
        Array.from({ length: n }, (_, k) => ({ time: c0Hour(20 + k), level: E.predictAtTime(c0Hour(20 + k), doses, POP) })),
        doses, POP, 1).rho);
    ok('rho is derived from the data, not a constant',
        new Set(rhos.map(r => r.toFixed(3))).size === 3, `rho = ${rhos.map(r => '+' + r.toFixed(3)).join(', ')}`);

    const p0 = E.laplacePosterior([], doses, POP, 1);
    ok('n=0 falls back to exactly the prior',
        near(p0.sdCL, priorSdCL, 1e-12) && near(p0.sdV, priorSdV, 1e-12) && p0.rho === 0);
}

section('Outlier alert (3 SD), decoupled from search bounds');
{
    const thrCL = E.PK_MODEL.OUTLIER_SD * Math.sqrt(E.PK_MODEL.OMEGA_CL);
    const thrV = E.PK_MODEL.OUTLIER_SD * Math.sqrt(E.PK_MODEL.OMEGA_V);
    ok('3 SD on etaCL is 1.154, not 1.26', near(thrCL, 1.154, 5e-4), `= ${thrCL.toFixed(4)}`);
    ok('3 SD on etaV is 0.900', near(thrV, 0.900, 5e-4), `= ${thrV.toFixed(4)}`);

    const slow = { CL: POP.CL * 0.25, V: POP.V, KA: POP.KA };
    const obsSlow = [10, 14, 18, 21, 24].map(d => ({ time: c0Hour(d), level: E.predictAtTime(c0Hour(d), doses, slow) }));
    const fitSlow = E.mapBayesian(POP, obsSlow, doses, 1);
    ok('alert fires for a genuine 3+ SD patient', fitSlow.outlierHit === true,
        `etaCL=${fitSlow.etaCL.toFixed(3)} (${(fitSlow.etaCL / Math.sqrt(E.PK_MODEL.OMEGA_CL)).toFixed(1)} SD)`);

    const fitNorm = E.mapBayesian(POP, [{ time: c0Hour(20), level: E.predictAtTime(c0Hour(20), doses, POP) }], doses, 1);
    ok('alert stays silent for a typical patient', fitNorm.outlierHit === false, `etaCL=${fitNorm.etaCL.toFixed(3)}`);

    // Nelder-Mead is unconstrained; the result must be clamped back into the box.
    ok('MAP estimate stays inside the search bounds',
        fitSlow.etaCL >= E.ETA_CL_BOUNDS[0] && fitSlow.etaCL <= E.ETA_CL_BOUNDS[1] &&
        fitSlow.etaV >= E.ETA_V_BOUNDS[0] && fitSlow.etaV <= E.ETA_V_BOUNDS[1]);
}

section('Dose series: held doses and bridging');
{
    const mk = (d, h, dose) => {
        const rd = TX.add(d, 'day').hour(h).minute(0);
        return { id: `x${d}-${h}`, recordDate: rd, dose, level: null, time: rd.diff(TX, 'hour', true) };
    };
    // A 0 mg entry means "dose held" and must survive gap-filling verbatim.
    const filled = E.fillHistoricalGaps([mk(0, 7, 3), mk(0, 19, 3), mk(1, 7, 0), mk(1, 19, 3), mk(4, 7, 3)], TX);
    const held = filled.find(d => d.recordDate.date() === 2 && d.recordDate.hour() === 7);
    ok('an explicit 0 mg dose is preserved, not re-imputed', held && held.dose === 0,
        `day-1 07:00 dose = ${held ? held.dose : 'MISSING'}`);
    ok('imputed gaps still carry the last positive dose',
        filled.filter(d => String(d.id).startsWith('inter-')).every(d => d.dose === 3));

    // Bridging a lagging dose log.
    const logged = steadyStateDoses(20);
    const lastLogged = logged[logged.length - 1].recordDate;
    const from = lastLogged.add(12, 'hour');
    const to = from.add(3, 'day');
    const bridge = E.buildBridgeDoses(logged, from, to, TX);
    ok('bridge covers the whole gap at 12 h spacing', bridge.length === 6, `${bridge.length} doses`);
    ok('bridge carries the existing regimen', bridge.every(d => d.dose === 3));
    ok('bridge is labelled for UI transparency', bridge.every(d => String(d.id).startsWith('bridge-')));
    ok('bridge is empty when the log is current', E.buildBridgeDoses(logged, from, from, TX).length === 0);

    // The bug this prevents: omitting the bridge models a drug-free washout.
    const tgt = to.diff(TX, 'hour', true);
    const withBridge = E.predictAtTime(tgt, [...logged, ...bridge], POP);
    const withoutBridge = E.predictAtTime(tgt, logged, POP);
    ok('omitting the bridge would under-predict badly', withoutBridge < 0.6 * withBridge,
        `${withoutBridge.toFixed(2)} vs ${withBridge.toFixed(2)} ng/mL (${(100 * (withoutBridge / withBridge - 1)).toFixed(0)}%)`);
}

section('C0 sampling convention (15 min pre-dose)');
{
    const afterAm = E.nextC0Time(TX.add(20, 'day').hour(7).minute(0));
    ok('next C0 after an AM dose is 18:45', afterAm.hour() === 18 && afterAm.minute() === 45, afterAm.format('DD-MMM HH:mm'));
    const afterPm = E.nextC0Time(TX.add(20, 'day').hour(19).minute(0));
    ok('next C0 after a PM dose is 06:45 the next day',
        afterPm.hour() === 6 && afterPm.minute() === 45 && afterPm.date() === 22, afterPm.format('DD-MMM HH:mm'));
    ok('nextC0Time always moves forward', E.nextC0Time(TX.add(20, 'day').hour(6).minute(50)).isAfter(TX.add(20, 'day').hour(6).minute(50)));

    // Why the timing warning exists.
    const trueC0 = E.predictAtTime(c0Hour(20), doses, POP);
    const late = E.predictAtTime(20 * 24 + 7.25, doses, POP);
    ok('a sample 15 min post-dose reads ~30% high', (late / trueC0 - 1) > 0.25,
        `${trueC0.toFixed(2)} -> ${late.toFixed(2)} (+${(100 * (late / trueC0 - 1)).toFixed(0)}%)`);

    const flagged = E.findPostDoseSamples([{ time: 20 * 24 + 7.25 }, { time: c0Hour(21) }], doses);
    ok('post-dose sample flagged, pre-dose sample not',
        flagged.length === 1 && flagged[0].n === 1 && flagged[0].minutesAfter === 15, JSON.stringify(flagged));
}

section('Accuracy metrics');
{
    const one = E.calculateAccuracyMetrics([{ time: c0Hour(20), level: 9.0 }], doses, POP, 1);
    ok('r2 is null at n=1 (no observed variance)', one.r2 === null);
    // `null < 0.5` is true in JS — this is the coercion that produced a false
    // "DO NOT TRUST FORECAST" banner on every single-trough patient.
    ok('null r2 must not be treated as a number below 0.5',
        !(typeof one.r2 === 'number' && one.r2 < 0.5));
    ok('rmse is still computed at n=1', typeof one.rmse === 'number' && isFinite(one.rmse));

    // MPE sign convention: error = obs - pred, so positive = under-prediction.
    const under = E.calculateAccuracyMetrics([{ time: c0Hour(20), level: 20 }], doses, POP, 1);
    ok('positive MPE means the model under-predicts', under.mpe > 0, `MPE = ${under.mpe.toFixed(1)}%`);
}

section('Rosendaal TTR interpolation');
{
    ok('fully in range', E.rosendaalFraction(6, 7, 5, 8) === 1);
    ok('fully out of range', E.rosendaalFraction(1, 2, 5, 8) === 0);
    ok('half in range', near(E.rosendaalFraction(5, 9, 5, 7), 0.5, 1e-9));
    ok('flat series inside range', E.rosendaalFraction(6, 6, 5, 8) === 1);
    ok('flat series outside range', E.rosendaalFraction(2, 2, 5, 8) === 0);
}

section('Therapeutic range schedule');
{
    const r = (m) => E.getTherapeuticRange(TX.add(m * 30.44, 'day'), TX);
    ok('month 0-1 targets 10-11', r(0.5).low === 10 && r(0.5).high === 11);
    ok('month 1-3 targets 7-9', r(2).low === 7 && r(2).high === 9);
    ok('month 3-6 targets 5-7', r(5).low === 5 && r(5).high === 7);
    ok('beyond 6 months targets 4-6', r(12).low === 4 && r(12).high === 6);
}

section('Covariate edge cases');
{
    // An unrecognized inhibitor string used to make CL silently NaN with no
    // guard — reachable via CSV import or the /api/patients payload, even
    // though the UI itself only ever emits none|moderate|strong.
    const garbage = E.getPopulationParameters({ weight: 60, mpa: '0', genotype: '33', bilirubin: 1.0, inhibitor: 'not-a-real-value' });
    ok('unrecognized inhibitor falls back to no effect, not NaN',
        isFinite(garbage.CL) && near(garbage.CL, E.PK_MODEL.TVCL * E.PK_MODEL.INDIAN_CL_SCALAR, 1e-6));

    // Locks the corrected bilirubin doc-comment values (the code's exponent
    // is 0.30; a previous comment quoted numbers implying ~0.277).
    const bilF = (bil) => E.getPopulationParameters({ weight: 60, mpa: '0', genotype: '33', bilirubin: bil, inhibitor: 'none' }).CL
        / (E.PK_MODEL.TVCL * E.PK_MODEL.INDIAN_CL_SCALAR);
    ok('Bil=2 -> 0.8123', near(bilF(2), 0.8123, 1e-4));
    ok('Bil=5 -> 0.6170', near(bilF(5), 0.6170, 1e-4));
    ok('Bil=10 -> 0.5012', near(bilF(10), 0.5012, 1e-4));
}

section('Verified case regression lock (raj bahadur, code-review-verified numbers)');
{
    // This locks the exact scenario independently hand-verified during a code
    // review: tx date 09-Jul-26, 70.5 kg, CYP3A5 *1/*3, MPA yes, LC-MS/MS,
    // 6 logged BID doses across 22 days (04-Jul .. 25-Jul), 3 measured troughs.
    // Every one of these numbers reproduced the app's displayed output exactly
    // when this test was written — any future drift here means the verified
    // math was touched.
    const patTx = dayjs('2026-07-09T00:00');
    const pat = { weight: 70.5, mpa: '1', genotype: '13', bilirubin: 1.0, inhibitor: 'none' };
    const pop = E.getPopulationParameters(pat);
    ok('population CL/F', near(pop.CL, 30.031, 0.01));
    ok('population V/F', near(pop.V, 815.45, 0.01));

    // Reconstruct the historyLog exactly as gatherData() would build it:
    // 6 logged dose rows, integer ids (user-entered rows), dayjs recordDate.
    const loggedDoseHours = { '-113': 1.5, '-101': 1.5, '127': 4.5, '139': 4.5, '391': 4.0, '403': 4.0 };
    const historyLog = Object.entries(loggedDoseHours).map(([t, dose], i) => ({
        id: i, recordDate: patTx.add(Number(t), 'hour'), dose, level: null, time: Number(t)
    }));
    const measuredLevels = [
        { time: 126.75, level: 4.5 },
        { time: 270.75, level: 12.9 },
        { time: 342.75, level: 11.1 }
    ];

    const filled = E.fillHistoricalGaps(historyLog, patTx);
    const imputed = filled.filter(d => String(d.id).startsWith('inter-'));
    ok('imputes exactly 38 missing BID doses across the logging gaps', imputed.length === 38, `${imputed.length}`);

    const ind = E.mapBayesian(pop, measuredLevels, filled, 1);
    ok('individual CL/F', near(ind.CL, 25.223, 0.01), `${ind.CL.toFixed(3)}`);
    ok('individual V/F', near(ind.V, 826.00, 0.01), `${ind.V.toFixed(2)}`);
    ok('RMSE', near(ind.rmse, 0.8709, 1e-3), `${ind.rmse.toFixed(4)}`);
    ok('MPE%', near(ind.mpe, -0.2846, 1e-3), `${ind.mpe.toFixed(4)}`);
    ok('MAPE%', near(ind.mape, 7.9828, 1e-3), `${ind.mape.toFixed(4)}`);

    const trough = E.predictAtTime(patTx.add(414.75, 'hour').diff(patTx, 'hour', true), filled, ind);
    ok('forecast trough at 26-Jul-26 06:45', near(trough, 11.766, 0.01), `${trough.toFixed(3)}`);

    // TTR: duration-weighted Rosendaal fraction across the 3 troughs against
    // the month 0-1 band (10-11 ng/mL).
    let ttrHours = 0, totalHours = 0;
    for (let i = 0; i < measuredLevels.length - 1; i++) {
        const dt = measuredLevels[i + 1].time - measuredLevels[i].time;
        totalHours += dt;
        ttrHours += E.rosendaalFraction(measuredLevels[i].level, measuredLevels[i + 1].level, 10, 11) * dt;
    }
    ok('TTR', near((ttrHours / totalHours) * 100, 7.94, 0.01), `${((ttrHours / totalHours) * 100).toFixed(2)}%`);
}

section('Change-point dose log: regimen carry-forward (never a drug-free day)');
{
    const mk = (d, h, dose, level = null) => {
        const rd = TX.add(d, 'day').hour(h).minute(0);
        return { id: `e${d}-${h}`, recordDate: rd, dose, level, time: rd.diff(TX, 'hour', true) };
    };
    const at = (filled, d, h) => filled.find(x => x.recordDate.isSame(TX.add(d, 'day').hour(h).minute(0)));

    // The real raj bahadur log: 1.5 mg BID entered on day 0, then ONLY the
    // 07:00 slot re-entered at 4.5 mg on day 10. The old rule carried the
    // 10-day-stale 1.5 mg into every PM slot, modelling 6 mg/day for a patient
    // on 9 mg/day.
    const oneShift = E.fillHistoricalGaps(
        [mk(0, 7, 1.5), mk(0, 19, 1.5), mk(10, 7, 4.5), mk(16, 7, 4.5)], TX);
    ok('a single revised entry replaces the stale opposite shift',
        at(oneShift, 10, 19).dose === 4.5, `day-10 19:00 = ${at(oneShift, 10, 19).dose} mg`);
    ok('the replaced regimen holds for every later imputed slot',
        oneShift.filter(d => d.recordDate.isAfter(TX.add(10, 'day').hour(7)))
            .every(d => d.dose === 4.5));
    ok('slots BEFORE the revision keep the original regimen',
        oneShift.filter(d => d.recordDate.isBefore(TX.add(10, 'day').hour(7)))
            .every(d => d.dose === 1.5));

    // ...but a deliberate asymmetric split entered for both shifts must survive
    // untouched, which is what the same-shift rule was there for originally.
    const split = E.fillHistoricalGaps([mk(0, 7, 5), mk(0, 19, 4.5), mk(6, 7, 5)], TX);
    ok('a deliberate AM/PM split is preserved, not flattened',
        split.filter(d => d.recordDate.hour() === 7).every(d => d.dose === 5) &&
        split.filter(d => d.recordDate.hour() === 19).every(d => d.dose === 4.5));

    // A 0 mg entry marks ONE held slot and never becomes the ongoing regimen.
    const heldRun = E.fillHistoricalGaps([mk(0, 7, 3), mk(0, 19, 3), mk(2, 7, 0), mk(6, 7, 3)], TX);
    ok('a held dose marks its own slot only', at(heldRun, 2, 7).dose === 0);
    ok('a held dose never defines the regimen for later slots',
        heldRun.filter(d => d.recordDate.isAfter(TX.add(2, 'day').hour(7))).every(d => d.dose === 3));

    // Extrapolation past the end of the log obeys the same rule.
    const tail = E.fillHistoricalGaps([mk(0, 7, 1.5), mk(0, 19, 1.5), mk(4, 7, 4.5)], TX);
    const ext = E.extrapolateDoses(tail, TX.add(7, 'day'), TX);
    ok('extrapolation carries the revised regimen, not the stale shift',
        ext.length > 0 && ext.every(d => d.dose === 4.5), `${ext.map(d => d.dose).join(',')}`);

    // Bridging a lagging log, same rule again.
    const bridged = E.buildBridgeDoses(
        [mk(0, 7, 1.5), mk(0, 19, 1.5), mk(4, 7, 4.5)],
        TX.add(4, 'day').hour(19), TX.add(6, 'day').hour(7), TX);
    ok('bridge carries the revised regimen across both shifts',
        bridged.length === 3 && bridged.every(d => d.dose === 4.5),
        `${bridged.map(d => d.dose).join(',')}`);
}

section('C/D ratio uses the reconstructed regimen');
{
    const mk = (d, h, dose, level = null) => {
        const rd = TX.add(d, 'day').hour(h).minute(0);
        return { id: `c${d}-${h}`, recordDate: rd, dose, level, time: rd.diff(TX, 'hour', true) };
    };

    // Levels drawn long after the last ENTERED dose row. Reading raw rows gave
    // a 0 mg denominator and printed '—'; the patient was on 4 mg BID
    // throughout.
    const log = [mk(0, 7, 4), mk(0, 19, 4), mk(9, 7, null, 8.0), mk(14, 7, null, 8.8)];
    const eff = E.buildEffectiveDoses(log, TX);

    const d9 = E.dailyDoseBefore(eff, mk(9, 7, null).time);
    const d14 = E.dailyDoseBefore(eff, mk(14, 7, null).time);
    ok('a level after the last entered dose still gets a full 24 h denominator',
        d9 === 8 && d14 === 8, `day 9 = ${d9} mg, day 14 = ${d14} mg`);
    ok('C/D is computed, not suppressed', near(8.0 / d9, 1.0, 1e-9) && near(8.8 / d14, 1.1, 1e-9));

    // Window boundaries: inclusive at exactly -24 h, exclusive at the sample
    // itself (a trough precedes its dose).
    const t = mk(9, 7, null).time;
    ok('the dose at exactly -24 h counts',
        eff.some(d => near(d.time, t - 24, 1e-9) && d.dose > 0) && d9 === 8);
    ok('the dose at the sample time itself does not count',
        E.dailyDoseBefore([{ time: t, dose: 4 }], t) === 0);

    // A genuinely held dose must lower the denominator.
    const withHold = E.buildEffectiveDoses(
        [mk(0, 7, 4), mk(0, 19, 4), mk(8, 19, 0), mk(9, 7, null, 8.0)], TX);
    ok('a held dose reduces the 24 h denominator',
        E.dailyDoseBefore(withHold, t) === 4, `${E.dailyDoseBefore(withHold, t)} mg`);
}

section('Starting dose suggestion (pre-Bayesian, population-based)');
{
    const range = { low: 10, high: 11 };
    const base = { weight: 70, mpa: '1', bilirubin: 1.0, inhibitor: 'none' };

    const fast = E.suggestStartingDose({ ...base, genotype: '11' }, range);
    const inter = E.suggestStartingDose({ ...base, genotype: '13' }, range);
    const slow = E.suggestStartingDose({ ...base, genotype: '33' }, range);

    ok('CYP3A5 expresser (*1/*1) needs a higher TDD than intermediate (*1/*3)', fast.tdd > inter.tdd,
        `*1/*1=${fast.tdd} *1/*3=${inter.tdd}`);
    ok('intermediate (*1/*3) needs a higher TDD than non-expresser (*3/*3)', inter.tdd > slow.tdd,
        `*1/*3=${inter.tdd} *3/*3=${slow.tdd}`);

    // CPIC (Birdwell 2015): expressers need ~1.5-2x the non-expresser starting
    // dose. This isn't a flat multiplier here — it falls out of the population
    // model's own CYP3A5 factor — so check the resulting ratio lands in range
    // rather than asserting an exact number.
    const ratio = fast.tdd / slow.tdd;
    ok('*1/*1 vs *3/*3 starting-dose ratio matches CPIC 1.5-2x guidance',
        ratio > 1.5 && ratio < 2.5, `ratio=${ratio.toFixed(2)}`);

    [fast, inter, slow].forEach(r => {
        ok(`predicted trough (${r.tdd} mg/day) lands inside the requested range`,
            r.predictedTrough >= range.low - 0.5 && r.predictedTrough <= range.high + 0.5,
            `${r.predictedTrough.toFixed(2)} ng/mL`);
        ok(`${r.tdd} mg/day is not flagged at the search boundary`, !r.atBoundary);
    });

    // Extreme covariate stack (strong inhibitor + severe hepatic impairment +
    // slow metabolizer + low weight) drives clearance so low that even the
    // 1 mg/day floor overshoots target — this must surface as atBoundary, not
    // silently return a dose that undertreats/overtreats without a flag.
    const extreme = E.suggestStartingDose(
        { weight: 30, mpa: '1', genotype: '33', bilirubin: 15, inhibitor: 'strong' }, range);
    ok('extreme low-clearance stack hits the dose floor and is flagged',
        extreme.atBoundary && extreme.tdd === 1.0, `tdd=${extreme.tdd} atBoundary=${extreme.atBoundary}`);
    ok('extreme low-clearance stack overshoots target even at the floor dose',
        extreme.predictedTrough > range.high, `${extreme.predictedTrough.toFixed(2)} ng/mL`);
}

done('pk-engine');
