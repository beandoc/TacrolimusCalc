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

section('C0 sampling convention (morning draw, 15 min pre-dose)');
{
    // Site protocol: TDM levels are drawn in the MORNING only, 15 min before the
    // 07:00 dose. The evening trough is never assayed, so nextC0Time must never
    // return 18:45 — an earlier revision did, which offered clinicians a default
    // timestamp for a draw that does not happen.
    const afterAm = E.nextC0Time(TX.add(20, 'day').hour(7).minute(0));
    ok('next C0 after an AM dose is 06:45 the NEXT day, not 18:45',
        afterAm.hour() === 6 && afterAm.minute() === 45 && afterAm.date() === 22, afterAm.format('DD-MMM HH:mm'));
    const afterPm = E.nextC0Time(TX.add(20, 'day').hour(19).minute(0));
    ok('next C0 after a PM dose is 06:45 the next day',
        afterPm.hour() === 6 && afterPm.minute() === 45 && afterPm.date() === 22, afterPm.format('DD-MMM HH:mm'));
    ok('nextC0Time always moves forward', E.nextC0Time(TX.add(20, 'day').hour(6).minute(50)).isAfter(TX.add(20, 'day').hour(6).minute(50)));
    // Sweep the whole day: every answer must be a 06:45 morning slot.
    let allMorning = true;
    for (let h = 0; h < 24; h++) {
        const c0 = E.nextC0Time(TX.add(20, 'day').hour(h).minute(30));
        if (c0.hour() !== 6 || c0.minute() !== 45) allMorning = false;
    }
    ok('every hour of the day resolves to a 06:45 morning draw', allMorning);

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

    // Values below reflect mapBayesian AFTER recency+robust weighting was
    // added (see "Erratic C/D volatility" section) — the 4.5 ng/mL trough
    // (day 5.3, closest to "now") is trusted more than the later 12.9/11.1
    // spike, shifting CL up slightly vs the old unweighted fit (was 25.223).
    const ind = E.mapBayesian(pop, measuredLevels, filled, 1);
    ok('individual CL/F', near(ind.CL, 25.945, 0.01), `${ind.CL.toFixed(3)}`);
    ok('individual V/F', near(ind.V, 822.29, 0.01), `${ind.V.toFixed(2)}`);
    ok('RMSE', near(ind.rmse, 0.8148, 1e-3), `${ind.rmse.toFixed(4)}`);
    ok('MPE%', near(ind.mpe, 3.0441, 1e-3), `${ind.mpe.toFixed(4)}`);
    ok('MAPE%', near(ind.mape, 8.8032, 1e-3), `${ind.mape.toFixed(4)}`);

    const trough = E.predictAtTime(patTx.add(414.75, 'hour').diff(patTx, 'hour', true), filled, ind);
    ok('forecast trough at 26-Jul-26 06:45', near(trough, 11.355, 0.01), `${trough.toFixed(3)}`);

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

    // The target band is in REPORTED units. A CMIA lab reads 1.08c + 0.55 higher
    // than LC-MS/MS, so the same band needs a LOWER dose. Omitting the assay
    // targeted the LC-MS/MS value and overshot at every CMIA centre: 6 mg/day
    // predicted 10.58 on the LC-MS/MS scale but 11.98 as actually reported.
    const band = { low: 10, high: 11 };
    const pdA = { weight: 60, genotype: '33', mpa: '1', bilirubin: 1.0, inhibitor: 'none' };
    const lcms = E.suggestStartingDose(pdA, band, 1);
    const cmia = E.suggestStartingDose(pdA, band, 2);
    ok('CMIA needs a lower starting dose than LC-MS/MS for the same band',
        cmia.tdd < lcms.tdd, `LC-MS/MS=${lcms.tdd} mg/day vs CMIA=${cmia.tdd} mg/day`);
    ok('CMIA predicted trough lands inside the band ON THE REPORTED SCALE',
        cmia.predictedTrough >= band.low && cmia.predictedTrough <= band.high,
        `${cmia.predictedTrough.toFixed(2)} ng/mL as reported`);
    ok('the LC-MS/MS dose would have been reported above the band by a CMIA lab',
        E.cmiaAdjust(lcms.predictedTrough) > band.high,
        `${lcms.predictedTrough.toFixed(2)} -> ${E.cmiaAdjust(lcms.predictedTrough).toFixed(2)} reported`);
    ok('bioassay falls back to the patient profile when the argument is omitted',
        E.suggestStartingDose({ ...pdA, bioassay: 2 }, band).tdd === cmia.tdd);
}

section('Time-varying covariates (Weight & Hematocrit)');
{
    const historyLog = [
        { id: 0, recordDate: dayjs('2026-07-11 07:00'), dose: 3, level: null, weight: 60, hematocrit: 35, time: 0 },
        { id: 1, recordDate: dayjs('2026-07-11 19:00'), dose: 3, level: null, time: 12 },
        { id: 2, recordDate: dayjs('2026-07-13 07:00'), dose: 3, level: 8.5, weight: 52, hematocrit: 28, time: 48 },
        { id: 3, recordDate: dayjs('2026-07-14 07:00'), dose: 3, level: null, time: 72 }
    ];

    const filled = E.fillHistoricalGaps(historyLog, dayjs('2026-07-11'));
    ok('fillHistoricalGaps carries forward initial weight=60', filled.find(d => d.time === 12).weight === 60);
    ok('fillHistoricalGaps updates to weight=52 after time=48', filled.find(d => d.time === 72).weight === 52);
    ok('fillHistoricalGaps updates to hematocrit=28 after time=48', filled.find(d => d.time === 72).hematocrit === 28);

    // Verify predictAtTime with time-varying weight
    const pop60 = E.getPopulationParameters({ weight: 60, genotype: 'unknown', mpa: '1', bilirubin: 1.0, inhibitor: 'none', transplantDate: dayjs('2026-07-11') });
    const predBaseline = E.predictAtTime(71, filled, pop60);
    
    // Constant 60kg filled doses
    const filledConstant60 = filled.map(d => ({ ...d, weight: 60 }));
    const predConstant60 = E.predictAtTime(71, filledConstant60, pop60);

    ok('weight reduction (60kg -> 52kg) lowers clearance, resulting in higher predicted trough concentration',
        predBaseline > predConstant60, `time-varying pred=${predBaseline.toFixed(2)} vs constant pred=${predConstant60.toFixed(2)}`);
}

section('Time-varying inhibitor (Azoles / Diltiazem)');
{
    // Verify fillHistoricalGaps carries forward inhibitor
    const historyLog = [
        { id: 0, recordDate: dayjs('2026-07-11 07:00'), dose: 3, level: null, weight: 60, hematocrit: 35, inhibitor: 'none', time: 0 },
        { id: 1, recordDate: dayjs('2026-07-11 19:00'), dose: 3, level: null, time: 12 },
        { id: 2, recordDate: dayjs('2026-07-13 07:00'), dose: 3, level: 8.5, weight: 60, hematocrit: 35, inhibitor: 'moderate', time: 48 },
        { id: 3, recordDate: dayjs('2026-07-14 07:00'), dose: 3, level: null, time: 72 }
    ];
    const filled = E.fillHistoricalGaps(historyLog, dayjs('2026-07-11'));
    ok('fillHistoricalGaps carries forward initial inhibitor=none', filled.find(d => d.time === 12).inhibitor === 'none');
    ok('fillHistoricalGaps updates to inhibitor=moderate after time=48', filled.find(d => d.time === 72).inhibitor === 'moderate');

    // Verify predictAtTime with time-varying inhibitor
    const popBase = E.getPopulationParameters({ weight: 60, genotype: 'unknown', mpa: '1', bilirubin: 1.0, inhibitor: 'none', transplantDate: dayjs('2026-07-11') });
    const predInhibitor = E.predictAtTime(71, filled, popBase);

    // Constant none inhibitor filled doses
    const filledConstantNone = filled.map(d => ({ ...d, inhibitor: 'none' }));
    const predConstantNone = E.predictAtTime(71, filledConstantNone, popBase);

    ok('co-administration of moderate inhibitor reduces clearance, resulting in higher predicted trough concentration',
        predInhibitor > predConstantNone, `time-varying pred=${predInhibitor.toFixed(2)} vs constant pred=${predConstantNone.toFixed(2)}`);
}

section('Levels the fit cannot use are excluded, not scored');
{
    // A level drawn before any logged dose predicts ~0. The OFV already ignores
    // it (flat +1000, constant in the etas), so scoring it as a 100% error made
    // the same point simultaneously ignored by the fit and counted against it —
    // RMSE 5.0 on an otherwise exact fit, which tripped the "Severe Model
    // Mismatch" banner with nothing on screen explaining why.
    const TX3 = dayjs('2026-01-01 00:00');
    const hist = [{ id: 0, recordDate: TX3.add(6, 'hour'), dose: null, level: 5.0, time: 6 }];
    for (let i = 0; i < 20; i++) {
        const rd = TX3.add(1, 'day').add(i * 12, 'hour');
        hist.push({ id: i + 1, recordDate: rd, dose: 3, level: null, time: rd.diff(TX3, 'hour', true) });
    }
    const pp = E.getPopulationParameters({ weight: 60, genotype: '33', mpa: '1', bilirubin: 1.0, inhibitor: 'none' });
    const doses = E.fillHistoricalGaps(hist, TX3);
    const tGood = TX3.add(9, 'day').hour(7).minute(0).subtract(E.SAMPLING_LEAD_MIN, 'minute').diff(TX3, 'hour', true);
    const exact = E.predictAtTime(tGood, doses, pp);

    const fit = E.mapBayesian(pp, [hist[0], { time: tGood, level: exact }], doses, 1);
    ok('the pre-dose-log level is flagged excluded', fit.nExcluded === 1 && fit.pe[0].excluded === true);
    ok('the usable level is still scored', fit.pe[1].excluded === false);
    ok('RMSE reflects only the usable levels', near(fit.rmse, 0, 1e-6), `rmse=${fit.rmse.toFixed(4)} (was 5.00)`);
    ok('MAPE reflects only the usable levels', near(fit.mape, 0, 1e-6), `mape=${fit.mape.toFixed(2)}%`);
    ok('an exact fit is no longer condemned as a model mismatch', !(fit.rmse > 2.5));

    // Every level unfittable -> no fit to describe, rather than a fake one.
    const allBad = E.mapBayesian(pp, [hist[0]], doses, 1);
    ok('all-unfittable returns null metrics, not a fabricated RMSE',
        allBad.rmse === null && allBad.nExcluded === 1);
}

section('Time-varying weight survives every params rebuild');
{
    // predictAtTime applies the per-dose weight ratio only when params.weight is
    // present. Four functions rebuild a params object from CL/V/KA, and each one
    // that enumerates the fields instead of spreading silently drops the
    // covariate — so the fit, the posterior and the CI band would model a
    // constant-weight patient while the plotted population curve modelled a
    // varying one. These assertions fail if any of them regresses.
    const TX2 = dayjs('2026-07-11 00:00');
    const pd = { weight: 60, genotype: '33', mpa: '1', bilirubin: 1.0, inhibitor: 'none' };
    const pop = E.getPopulationParameters(pd);

    // 60 kg oedematous at transplant, down to 45 kg by day 10, 3 mg BID.
    const hist = [];
    for (let i = 0; i < 44; i++) {
        const rd = TX2.add(i * 12, 'hour');
        const e = { id: i, recordDate: rd, dose: 3, level: null, time: i * 12 };
        if (i === 0) e.weight = 60;
        if (i === 20) e.weight = 45;
        hist.push(e);
    }
    const doses = E.fillHistoricalGaps(hist, TX2);
    const tObs = 21 * 24 - 0.25;

    // The strongest available check: feed the model an observation it generated
    // itself. Any covariate the fit cannot see shows up as a non-zero eta.
    const truth = E.predictAtTime(tObs, doses, pop);
    const fit = E.mapBayesian(pop, [{ time: tObs, level: truth }], doses, 1);
    ok('MAP recovers etaCL = 0 on self-generated data (weight reaches the OFV)',
        near(fit.etaCL, 0, 1e-3), `etaCL=${fit.etaCL.toFixed(4)} (was -0.1508, CL 14% low)`);
    ok('indParams carries weight through to the individual curve',
        fit.weight === pop.weight, `indParams.weight=${fit.weight}`);

    // The n=0 branch spreads popParams, so n=0 and a perfectly-fitting n=1 must
    // put the individual curve in the same place.
    const tLate = 43 * 12 - 0.25;
    const curve0 = E.generateCurve([tLate], doses, E.mapBayesian(pop, [], doses, 1), 1)[0];
    const curve1 = E.generateCurve([tLate], doses, fit, 1)[0];
    ok('individual curve does not jump when the first level is added',
        near(curve0, curve1, 0.01), `n=0: ${curve0.toFixed(3)} vs n=1: ${curve1.toFixed(3)}`);

    // Laplace posterior and MC band must be built on the same patient as the fit.
    const post = E.laplacePosterior([{ time: tObs, level: truth }], doses, fit, 1);
    ok('posterior sd(etaCL) is tighter than the prior', post.sdCL < Math.sqrt(E.PK_MODEL.OMEGA_CL));
    const band = E.monteCarloCI([tLate], doses, fit, [{ time: tObs, level: truth }], 1, 400)[0];
    ok('MC band brackets the weight-aware individual curve',
        band.p5 <= curve1 && curve1 <= band.p95,
        `p5=${band.p5.toFixed(2)} curve=${curve1.toFixed(3)} p95=${band.p95.toFixed(2)}`);
}

section('calculateIPV (C/D ratio variability)');
{
    // Constant C/D ratio → zero variability, regardless of the absolute level.
    const flatDoses = [{ id: 0, time: -12, dose: 4 }, { id: 1, time: 0, dose: 4 }];
    const flatLevels = [{ time: 12, level: 8 }, { time: 24, level: 8 }, { time: 36, level: 8 }];
    const flat = E.calculateIPV(flatLevels, [...flatDoses, { id: 2, time: 12, dose: 4 }, { id: 3, time: 24, dose: 4 }]);
    ok('constant C/D ratio gives IPV ~0%', near(flat.ipv, 0, 1e-6), `ipv=${flat.ipv}`);

    // Fewer than 2 usable pairs → null, not NaN or a divide-by-zero.
    ok('n=0 levels returns null ipv, not NaN', E.calculateIPV([], []).ipv === null);
    ok('n=1 usable pair returns null ipv (need >=2 for a variance)',
        E.calculateIPV([{ time: 12, level: 8 }], flatDoses).ipv === null);

    // A level with no dose in the prior 24h contributes no C/D pair — matches
    // dailyDoseBefore's own `daily > 0` gate, not counted as a "0/0" pair.
    const withGap = E.calculateIPV(
        [{ time: -1000, level: 5 }, { time: 12, level: 8 }, { time: 24, level: 8 }],
        flatDoses.concat([{ id: 2, time: 12, dose: 4 }])
    );
    ok('an observation with no preceding dose is excluded, not counted as n', withGap.n === 2, `n=${withGap.n}`);

    // Regression case: arun_chougle (the real patient whose full-history MAP
    // fit motivated the volatility banner below). C/D ratios reconstructed
    // from the actual dose log via dailyDoseBefore: 2.63, 1.08, 1.94, 2.30,
    // 1.11, 1.03, 1.24, 0.56 — independently verified against the app's own
    // reported C/D column for this patient. Known IPV ≈ 48%, "High" per this
    // app's >40% threshold (Sapir-Pichhadze 2014).
    const cdRatios = [2.63, 1.08, 1.94, 2.30, 1.11, 1.03, 1.24, 0.56];
    const chougleLevels = cdRatios.map((cd, i) => ({ time: i * 100, level: cd * 8 }));
    const chougleDoses = cdRatios.map((_, i) => ({ id: i, time: i * 100 - 1, dose: 8 }));
    const chougle = E.calculateIPV(chougleLevels, chougleDoses);
    ok('arun_chougle case: IPV lands in the "High" band (>40%)', chougle.ipv > 40, `ipv=${chougle.ipv.toFixed(1)}%`);
    ok('arun_chougle case: IPV matches hand-computed CV% within rounding', near(chougle.ipv, 48.4, 1.0), `ipv=${chougle.ipv.toFixed(1)}%`);
}

section('classifyMetabolizer (Thölking 2014 C/D thresholds)');
{
    // Published cutoffs: fast <1.05, intermediate 1.05-2.0, slow >2.0.
    ok('C/D 0.96 is FAST, not "Normal"', E.classifyMetabolizer(0.96) === 'Fast Metabolizer');
    ok('C/D 1.04 is still fast (just under the cutoff)', E.classifyMetabolizer(1.04) === 'Fast Metabolizer');
    ok('C/D 1.05 is intermediate (boundary is inclusive above)', E.classifyMetabolizer(1.05) === 'Intermediate Metabolizer');
    ok('C/D 2.0 is intermediate (upper boundary inclusive)', E.classifyMetabolizer(2.0) === 'Intermediate Metabolizer');
    ok('C/D 2.01 is slow', E.classifyMetabolizer(2.01) === 'Slow Metabolizer');
    // Never label off a value that cannot be a ratio.
    ok('null C/D returns null, not a label', E.classifyMetabolizer(null) === null);
    ok('zero C/D returns null, not "Fast"', E.classifyMetabolizer(0) === null);
    ok('NaN C/D returns null', E.classifyMetabolizer(NaN) === null);

    // Regression case: CASE-B (de-identified). 15 troughs, C/D reconstructed from the real
    // dose-change log (each logged date carried forward until revised) and
    // verified against the app's own C/D column to 2dp. Clinical ground truth:
    // could not hold a therapeutic level even at 9 mg/day and was switched off
    // tacrolimus — a fast metabolizer by outcome.
    //
    // The arithmetic mean is 0.96 and the median 0.97, both below 1.05, so the
    // OLD 0.9 threshold reported "Normal Metabolizer" for this patient. That is
    // the bug this section locks down.
    const caseBCd = [1.02, 0.21, 0.40, 0.55, 0.79, 0.98, 1.27, 1.05, 0.89, 0.97, 0.69, 1.57, 1.09, 0.62, 2.37];
    const caseBLevels = caseBCd.map((cd, i) => ({ time: i * 100, level: cd * 8 }));
    const caseBDoses = caseBCd.map((_, i) => ({ id: i, time: i * 100 - 1, dose: 8 }));
    const caseB = E.calculateIPV(caseBLevels, caseBDoses);
    ok('CASE-B: median C/D is reported alongside the mean',
        typeof caseB.median === 'number', `median=${caseB.median}`);
    ok('CASE-B: median C/D 0.97 classifies as Fast Metabolizer',
        E.classifyMetabolizer(caseB.median) === 'Fast Metabolizer',
        `median=${caseB.median.toFixed(2)} -> ${E.classifyMetabolizer(caseB.median)}`);
    ok('CASE-B: the OLD 0.9 cutoff would have mislabelled this patient as Normal',
        caseB.mean > 0.9 && caseB.median > 0.9,
        `mean=${caseB.mean.toFixed(2)} median=${caseB.median.toFixed(2)}`);
    ok('CASE-B: IPV is "High" (>40%)', caseB.ipv > 40, `ipv=${caseB.ipv.toFixed(1)}%`);

    // The median must resist outliers that would drag an arithmetic mean over a
    // threshold — the whole reason the label switched statistic.
    const skew = E.calculateIPV(
        [0.8, 0.8, 0.8, 0.8, 9.0].map((cd, i) => ({ time: i * 100, level: cd * 8 })),
        [0.8, 0.8, 0.8, 0.8, 9.0].map((_, i) => ({ id: i, time: i * 100 - 1, dose: 8 }))
    );
    ok('one extreme outlier drags the mean over 1.05 but not the median',
        skew.mean > 1.05 && skew.median === 0.8, `mean=${skew.mean.toFixed(2)} median=${skew.median}`);
    ok('median-based label stays Fast despite the outlier',
        E.classifyMetabolizer(skew.median) === 'Fast Metabolizer');

    // Even n: median must average the two middle values.
    const evenN = E.calculateIPV(
        [1.0, 2.0, 3.0, 4.0].map((cd, i) => ({ time: i * 100, level: cd * 8 })),
        [1.0, 2.0, 3.0, 4.0].map((_, i) => ({ id: i, time: i * 100 - 1, dose: 8 }))
    );
    ok('even n averages the two middle C/D values', near(evenN.median, 2.5, 1e-9), `median=${evenN.median}`);
}

section('Erratic C/D volatility: automatic recency+robust weighting replaces the manual toggle');
{
    // Reproduces the arun_chougle case end-to-end through the real doses and
    // levels (not synthetic C/D ratios): 6 weeks of dosing changes, 8 troughs
    // oscillating instead of trending. Previously this required a manual
    // "recent-only fit" toggle (removed) to get a forecast connected to the
    // most recent, most clinically relevant trough; mapBayesian now does
    // this automatically via recencyWeights + combinedWeights (age-based
    // decay, split by maturation regime, refined by IRLS Tukey biweight).
    //
    // Historical reference, UNWEIGHTED (every observation counted equally,
    // the behavior before this feature): CL 16.77 L/hr, V 335 L, RMSE 6.89,
    // forecast 17.03 ng/mL — wildly disconnected from the last observed
    // trough of 4.5. That number is not asserted below (it belongs to a
    // fit mode this codebase no longer has); it is kept here only so a
    // future reader can see what changed.
    const epoch = dayjs('2026-06-23 00:00');
    const doseRows = [
        ['2026-06-23 19:00', 5.0], ['2026-06-24 07:00', 5.0], ['2026-06-24 19:00', 5.0],
        ['2026-06-26 07:00', 4.5], ['2026-06-26 19:00', 4.5],
        ['2026-07-06 07:00', 4.0], ['2026-07-06 19:00', 4.0],
        ['2026-07-10 07:00', 3.5], ['2026-07-10 19:00', 3.5],
        ['2026-07-14 07:00', 3.5], ['2026-07-14 19:00', 4.0],
        ['2026-07-20 07:00', 4.5], ['2026-07-20 19:00', 4.0],
        ['2026-07-27 07:00', 4.0], ['2026-07-27 19:00', 4.0],
        ['2026-08-04 07:00', 4.5], ['2026-08-04 19:00', 4.5],
    ];
    const levelRows = [
        ['2026-06-25 06:45', 26.3], ['2026-07-01 06:45', 9.7], ['2026-07-05 06:45', 17.5],
        ['2026-07-09 06:45', 18.4], ['2026-07-13 06:45', 7.8], ['2026-07-19 06:45', 7.7],
        ['2026-07-26 06:45', 10.5],
        ['2026-08-03 06:45', 4.5, 66, 40.1], // weight/Hct override entered on this lab row
    ];
    const rows = [
        ...doseRows.map(([dt, dose]) => ({ datetime: dt, dose, level: null, weight: null, hematocrit: null })),
        ...levelRows.map(([dt, level, weight, hct]) => ({ datetime: dt, dose: null, level, weight: weight ?? null, hematocrit: hct ?? null })),
    ].sort((a, b) => dayjs(a.datetime).diff(dayjs(b.datetime)));
    const historyLog = rows.map((r, i) => {
        const rd = dayjs(r.datetime);
        return { id: i, recordDate: rd, dose: r.dose, level: r.level, weight: r.weight, hematocrit: r.hematocrit, time: rd.diff(epoch, 'hour', true) };
    });

    const patientData = { weight: 68, genotype: '33', mpa: '1', bioassay: 1, bilirubin: 1.0, inhibitor: 'none' };
    const popParams = E.getPopulationParameters(patientData);
    const measuredLevels = historyLog.filter(e => e.level !== null);
    const filledHistory = E.fillHistoricalGaps(historyLog, epoch);
    const predDate = dayjs('2026-08-13 06:45');
    const allDoses = [...filledHistory, ...E.extrapolateDoses(filledHistory, predDate.add(14, 'day'), epoch)];

    const indParams = E.mapBayesian(popParams, measuredLevels, allDoses, 1);
    ok('weighted fit on the SAME 8-level history now converges to a faster clearance',
        near(indParams.CL, 33.46, 0.5), `CL=${indParams.CL.toFixed(2)} (unweighted reference was 16.77)`);

    const predTime = predDate.diff(epoch, 'hour', true);
    const fullPred = E.generateCurve([predTime], allDoses, indParams, 1)[0];
    ok('forecast is now close to the last observed trough (4.5), not the unweighted 17.03',
        near(fullPred, 8.84, 0.3), `forecast=${fullPred.toFixed(2)} last observed=4.5`);
    ok('this is a genuine improvement, not a coincidence of the new numbers',
        Math.abs(fullPred - 4.5) < Math.abs(17.03 - 4.5) / 2,
        `|forecast-4.5|=${Math.abs(fullPred - 4.5).toFixed(2)} vs unweighted |17.03-4.5|=${(17.03 - 4.5).toFixed(2)}`);

    const { ipv } = E.calculateIPV(measuredLevels, allDoses);
    ok('IPV on the real data is still "High" (>40%) — the input hasn\'t changed, only the fit', ipv > 40, `ipv=${ipv.toFixed(1)}%`);

    const weights = E.combinedWeights(measuredLevels, allDoses, indParams, 1, E.recencyWeights(measuredLevels));
    ok('the two most volatile mid-history levels (17.5, 18.4) end up least trusted',
        weights[measuredLevels.findIndex(l => l.level === 17.5)] < 0.1 &&
        weights[measuredLevels.findIndex(l => l.level === 18.4)] < 0.1,
        `weights=${weights.map(w => w.toFixed(3)).join(',')}`);
    ok('the most recent level (4.5) is among the most trusted',
        weights[measuredLevels.findIndex(l => l.level === 4.5)] > 0.3,
        `4.5's weight=${weights[measuredLevels.findIndex(l => l.level === 4.5)].toFixed(3)}`);
}

section('Per-patient regime selection (walk-forward, not a global rule)');
{
    // Regimes must stay distinguishable, or selection is meaningless.
    const early = [0, 5, 10, 15].map(d => ({ time: d * 24, level: 8 }));
    const wAdaptive = E.recencyWeights(early, 'adaptive');
    const wStationary = E.recencyWeights(early, 'stationary');
    ok('adaptive discounts the oldest level far harder than stationary does',
        wAdaptive[0] < wStationary[0] / 3,
        `adaptive=${wAdaptive[0].toFixed(3)} stationary=${wStationary[0].toFixed(3)}`);
    ok('both regimes still normalize to sum = N (prior influence unchanged)',
        near(wAdaptive.reduce((s, w) => s + w, 0), 4, 1e-9) &&
        near(wStationary.reduce((s, w) => s + w, 0), 4, 1e-9));

    // An unknown regime name must not silently produce NaN weights.
    const wBogus = E.recencyWeights(early, 'no-such-regime');
    ok('an unknown regime name falls back to stationary, not NaN',
        wBogus.every(w => isFinite(w)) && near(wBogus[0], wStationary[0], 1e-9));

    // Cold start: too few levels to measure anything, so fall back to the
    // day-based guess and SAY that is what happened.
    const tooFew = [{ time: 24, level: 8 }, { time: 48, level: 8 }];
    const cold = E.selectRecencyRegime(POP, tooFew, doses, 1);
    ok('below the minimum level count, selection falls back to the day heuristic',
        cold.nTrials === 0 && /post-transplant day/.test(cold.basis), cold.basis);
    ok('the cold-start fallback still names a usable regime',
        cold.regime in E.PK_MODEL.RECENCY_REGIMES, cold.regime);

    // A patient whose clearance genuinely SHIFTED partway through: the
    // stationary regime cannot represent it, so selection should pick
    // adaptive on the evidence.
    const shiftDoses = steadyStateDoses(60);
    const slow = { CL: POP.CL * Math.exp(-0.35), V: POP.V, KA: POP.KA };
    const fast = { CL: POP.CL * Math.exp(0.45), V: POP.V, KA: POP.KA };
    const shifted = [20, 24, 28, 32].map(d => ({ time: c0Hour(d), level: E.predictAtTime(c0Hour(d), shiftDoses, slow) }))
        .concat([40, 44, 48].map(d => ({ time: c0Hour(d), level: E.predictAtTime(c0Hour(d), shiftDoses, fast) })));
    const shiftSel = E.selectRecencyRegime(POP, shifted, shiftDoses, 1);
    ok('a genuine mid-history clearance shift selects the adaptive regime',
        shiftSel.regime === 'adaptive',
        `chose=${shiftSel.regime} adaptive=${shiftSel.scores.adaptive.toFixed(3)} stationary=${shiftSel.scores.stationary.toFixed(3)}`);
    ok('selection reports it measured, not guessed',
        /walk-forward/.test(shiftSel.basis) && shiftSel.nTrials > 0, `nTrials=${shiftSel.nTrials}`);

    // A true tie is no longer enough to keep the conservative regime. The
    // backend walk-forward audit showed the old stationary-wins-ties rule cost
    // accuracy by missing drifting patients, so adaptive now wins ties.
    const stable = [20, 24, 28, 32, 36, 40, 44].map(d =>
        ({ time: c0Hour(d), level: E.predictAtTime(c0Hour(d), shiftDoses, POP) }));
    const stableSel = E.selectRecencyRegime(POP, stable, shiftDoses, 1);
    ok('adaptive wins a walk-forward tie',
        stableSel.regime === 'adaptive',
        `chose=${stableSel.regime} adaptive=${stableSel.scores.adaptive.toFixed(3)} stationary=${stableSel.scores.stationary.toFixed(3)}`);

    // Stationary still wins when it is genuinely better, not merely tied.
    const noisyStable = stable.concat([
        { time: c0Hour(48), level: E.predictAtTime(c0Hour(48), shiftDoses, POP) * 0.65 },
        { time: c0Hour(52), level: E.predictAtTime(c0Hour(52), shiftDoses, POP) * 1.35 },
        { time: c0Hour(56), level: E.predictAtTime(c0Hour(56), shiftDoses, POP) }
    ]);
    const noisyStableSel = E.selectRecencyRegime(POP, noisyStable, shiftDoses, 1);
    ok('stationary still wins when its walk-forward score is lower',
        noisyStableSel.regime === 'stationary',
        `chose=${noisyStableSel.regime} adaptive=${noisyStableSel.scores.adaptive.toFixed(3)} stationary=${noisyStableSel.scores.stationary.toFixed(3)}`);

    // The walk-forward trials double as the displayed per-patient track
    // record, so they must carry enough detail to render it honestly.
    const t = shiftSel.trials[shiftSel.regime];
    ok('trials expose observed/predicted/error per prediction for the UI record',
        Array.isArray(t) && t.length > 0 &&
        t.every(r => r.observed > 0 && r.predicted > 0 && isFinite(r.pctError) && r.absPctError >= 0),
        `${t.length} trials`);
    ok('absPctError is the magnitude of pctError',
        t.every(r => near(r.absPctError, Math.abs(r.pctError), 1e-9)));
}

section('forecastResolution: can the number settle the clinical question?');
{
    const band = { low: 7, high: 9 };
    // Narrow uncertainty fully inside the band -> actionable.
    ok('a tight forecast inside the band resolves as "in"',
        E.forecastResolution(8, band, 5).status === 'in');
    // The same forecast with realistic (~26%) error spans 5.9-10.1 and
    // cannot distinguish in-range from out — the common real case.
    const real = E.forecastResolution(8, band, 26);
    ok('the same forecast at measured error straddles the band and is unresolved',
        real.status === 'unresolved' && !real.resolves,
        `${real.lo.toFixed(1)}-${real.hi.toFixed(1)} vs band ${band.low}-${band.high}`);
    ok('clearly high forecasts still resolve as above',
        E.forecastResolution(20, band, 26).status === 'above');
    ok('clearly low forecasts still resolve as below',
        E.forecastResolution(2, band, 26).status === 'below');

    // Missing/invalid patient-specific error must fall back to the measured
    // population figure, never to a silently perfect 0%.
    const fallback = E.forecastResolution(8, band, null);
    ok('a missing patient error falls back to the population backtest MAPE',
        near(fallback.err, E.POPULATION_BACKTEST_MAPE, 1e-9), `err=${fallback.err}`);
    ok('a zero/NaN error does not fake certainty',
        E.forecastResolution(8, band, 0).err === E.POPULATION_BACKTEST_MAPE &&
        E.forecastResolution(8, band, NaN).err === E.POPULATION_BACKTEST_MAPE);
    ok('degenerate inputs return null rather than a bogus verdict',
        E.forecastResolution(0, band, 10) === null && E.forecastResolution(8, null, 10) === null);
}

section('Walk-forward selection is leak-free (the accuracy numbers depend on it)');
{
    // Every accuracy figure this project quotes — POPULATION_BACKTEST_MAPE, the
    // per-patient track record shown to clinicians, forecastResolution's error
    // width — comes out of selectRecencyRegime's walk-forward loop. If a target
    // level could reach the fit that predicts it, all of those numbers would be
    // optimistic and nothing downstream would fail visibly.
    //
    // Tested behaviourally rather than by inspecting `sorted.slice(0, i)`:
    // change ONLY a held-out observation, and its own out-of-sample prediction
    // must not move. That is the definition of no leakage, and it keeps holding
    // if the implementation is rewritten.
    const lkDoses = steadyStateDoses(60);
    const truth = { CL: POP.CL * Math.exp(0.2), V: POP.V, KA: POP.KA };
    const mkLevels = () => [12, 16, 20, 24, 28, 32, 36]
        .map(d => ({ time: c0Hour(d), level: E.predictAtTime(c0Hour(d), lkDoses, truth) }));

    const baseLevels = mkLevels();
    const basePick = E.selectRecencyRegime(POP, baseLevels, lkDoses, 1);
    ok('the leak probe actually exercises the walk-forward path',
        basePick.nTrials >= 2 && /walk-forward/.test(basePick.basis),
        `nTrials=${basePick.nTrials}`);

    // Perturb the LAST level — a held-out target in every regime's trial set.
    const bumped = mkLevels();
    bumped[bumped.length - 1].level *= 3;
    const bumpedPick = E.selectRecencyRegime(POP, bumped, lkDoses, 1);

    const lastOf = sel => {
        const rows = sel.trials && sel.trials.stationary;
        return rows && rows.length ? rows[rows.length - 1] : null;
    };
    const a = lastOf(basePick), b = lastOf(bumpedPick);
    ok('the perturbed level is genuinely the final held-out target',
        a && b && near(a.time, b.time, 1e-9) && b.observed > a.observed * 2.5,
        a && b ? `observed ${a.observed.toFixed(2)} -> ${b.observed.toFixed(2)}` : 'no trial rows');
    ok('tripling a held-out level does NOT move its own prediction (no leakage)',
        a && b && near(a.predicted, b.predicted, 1e-9),
        a && b ? `pred ${a.predicted.toFixed(6)} vs ${b.predicted.toFixed(6)}` : '');

    // Guard against the assertion above passing trivially: an EARLIER level is
    // legitimate fit input for that same target, so it must change the number.
    const early = mkLevels();
    early[0].level *= 3;
    const earlyPick = E.selectRecencyRegime(POP, early, lkDoses, 1);
    const c = lastOf(earlyPick);
    ok('perturbing an earlier level DOES move the prediction (probe is live)',
        a && c && !near(a.predicted, c.predicted, 1e-6),
        a && c ? `pred ${a.predicted.toFixed(4)} vs ${c.predicted.toFixed(4)}` : '');

    // A future DOSE must not reach a prediction either: predictAtTime filters
    // `d.time < t`, and buildEffectiveDoses/fillHistoricalGaps must not smuggle
    // one in ahead of the sample.
    const t = c0Hour(20);
    const withFuture = [...lkDoses, { id: 'future', recordDate: TX.add(50, 'day'), dose: 99, level: null, time: 50 * 24 }];
    ok('a dose recorded after the sample time cannot change that prediction',
        near(E.predictAtTime(t, lkDoses, truth), E.predictAtTime(t, withFuture, truth), 1e-12));
}

done('pk-engine');
