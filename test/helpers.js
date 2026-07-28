// Shared test helpers: a tiny assertion runner and PK fixtures.
const dayjs = require('dayjs');
const E = require('../pk-engine');

let pass = 0, fail = 0;
const results = [];

function ok(name, cond, detail = '') {
    if (cond) { pass++; console.log(`  PASS  ${name}${detail ? '  ' + detail : ''}`); }
    else { fail++; console.log(`  FAIL  ${name}  ${detail}`); }
    results.push({ name, cond });
}
const near = (a, b, tol) => Math.abs(a - b) <= tol;
const section = title => console.log(`\n${title}`);

function done(suite) {
    console.log(`\n${fail === 0 ? 'ALL PASS' : 'FAILURES'}: ${pass} passed, ${fail} failed  (${suite})\n`);
    process.exit(fail === 0 ? 0 : 1);
}

// ── Fixtures ────────────────────────────────────────────────────────────────
const TX = dayjs('2025-01-01T00:00');
const POP = { CL: E.PK_MODEL.TVCL * E.PK_MODEL.INDIAN_CL_SCALAR, V: E.PK_MODEL.TVV, KA: E.PK_MODEL.TVKA };

// `days` of 3 mg BID at the standard 07:00 / 19:00 slots.
function steadyStateDoses(days = 30, mg = 3) {
    const out = [];
    for (let d = 0; d < days; d++) {
        for (const h of [E.PK_MODEL ? 7 : 7, 19]) {
            const rd = TX.add(d, 'day').hour(h).minute(0);
            out.push({ id: `d${d}-${h}`, recordDate: rd, dose: mg, level: null, time: rd.diff(TX, 'hour', true) });
        }
    }
    return out;
}

// Hour-offset of the C0 sample on a given post-transplant day (06:45).
const c0Hour = day => day * 24 + 7 - E.SAMPLING_LEAD_MIN / 60;

// Exact Laplace posterior, computed independently of pk-engine's implementation.
function referencePosterior(times, doses, ind) {
    const { OMEGA_CL, OMEGA_V, SIGMA } = E.PK_MODEL;
    const h = 1e-5, s2 = SIGMA * SIGMA;
    let Sxx = 0, Sxy = 0, Syy = 0;
    for (const T of times) {
        const f = (a, b) => E.predictAtTime(T, doses, { CL: ind.CL * Math.exp(a), V: ind.V * Math.exp(b), KA: ind.KA });
        const c = f(0, 0);
        const gx = (f(h, 0) - f(-h, 0)) / (2 * h) / c;
        const gy = (f(0, h) - f(0, -h)) / (2 * h) / c;
        Sxx += gx * gx; Sxy += gx * gy; Syy += gy * gy;
    }
    const a = 1 / OMEGA_CL + Sxx / s2, b = Sxy / s2, c = 1 / OMEGA_V + Syy / s2;
    const det = a * c - b * b, vCL = c / det, vV = a / det;
    return { sdCL: Math.sqrt(vCL), sdV: Math.sqrt(vV), rho: (-b / det) / Math.sqrt(vCL * vV) };
}

module.exports = { E, dayjs, ok, near, section, done, TX, POP, steadyStateDoses, c0Hour, referencePosterior };
