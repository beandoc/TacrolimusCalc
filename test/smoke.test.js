// Integration smoke test: drives the real runForecast() out of index.html.
//
// Deliberately thin. Unit tests in pk-engine.test.js cover the maths; this only
// catches wiring failures that unit tests structurally cannot — signature drift
// between runForecast and the engine, and NaNs reaching the chart. It asserts
// on values, never on markup, so restyling the UI does not break it.

const fs = require('fs');
const path = require('path');
const vm = require('vm');
const dayjs = require('dayjs');
const E = require('../pk-engine');
const { ok, done, section } = require('./helpers');

const root = path.join(__dirname, '..');
const lines = fs.readFileSync(path.join(root, 'index.html'), 'utf8').split('\n');
const start = lines.findIndex(l => l.trim() === '<script>');
const end = lines.length - 1 - [...lines].reverse().findIndex(l => l.includes('</script>'));

// Minimal DOM so the page script can evaluate and render.
const html = {};
const mkEl = id => ({
    id, style: {}, value: '', dataset: {}, classList: { add() { }, remove() { } }, addEventListener() { },
    get innerHTML() { return html[id] || ''; }, set innerHTML(v) { html[id] = v; },
    get textContent() { return html[id] || ''; }, set textContent(v) { html[id] = v; },
    getContext: () => ({ save() { }, restore() { }, fillRect() { }, strokeRect() { }, setLineDash() { } })
});
const els = new Proxy({}, { get: (t, k) => (t[k] = t[k] || mkEl(k)) });

let chart = null;
const ctx = Object.assign(Object.create(null), {
    dayjs, console, Math, JSON, String, Array, Number, Object, Set, Promise, Date,
    isNaN, isFinite, parseFloat, parseInt, setTimeout,
    localStorage: { getItem: () => null, setItem() { } },
    document: {
        getElementById: id => els[id], querySelectorAll: () => [],
        addEventListener() { }, createElement: () => mkEl('tmp')
    },
    window: { location: { hostname: 'localhost', protocol: 'http:' }, addEventListener() { } },
    fetch: () => Promise.reject(new Error('offline')),
    Chart: function (c, cfg) { chart = cfg; this.destroy = () => { }; },
    flatpickr: () => ({}), alert() { }, confirm: () => false,
    URL: { createObjectURL: () => '', revokeObjectURL() { } }, Blob: function () { }
});
ctx.globalThis = ctx;
Object.assign(ctx, E);                       // engine globals, as the browser sees them
vm.createContext(ctx);
vm.runInContext(lines.slice(start + 1, end).join('\n'), ctx);

const TX = dayjs().subtract(40, 'day').startOf('day');
const POP = { CL: E.PK_MODEL.TVCL * E.PK_MODEL.INDIAN_CL_SCALAR, V: E.PK_MODEL.TVV, KA: E.PK_MODEL.TVKA };
const PATIENT = {
    weight: 60, hct: 35, mpa: '1', bioassay: 1, genotype: '33',
    bilirubin: 1.0, inhibitor: 'none', transplantDate: TX
};

function history(nLevels) {
    const log = [];
    for (let t = TX.hour(7).minute(0); t.isBefore(dayjs()); t = t.add(12, 'hour')) {
        log.push({ id: log.length, recordDate: t, dose: 3, level: null, time: t.diff(TX, 'hour', true) });
    }
    for (let k = 0; k < nLevels; k++) {
        const st = TX.add(25 + k * 3, 'day').hour(6).minute(45);
        const time = st.diff(TX, 'hour', true);
        log.push({ id: 'L' + k, recordDate: st, dose: null, level: E.predictAtTime(time, log, POP) * 1.05, time });
    }
    return log.sort((a, b) => a.time - b.time);
}

async function forecast(nLevels) {
    Object.keys(html).forEach(k => delete html[k]);
    chart = null;
    await vm.runInContext('runForecast', ctx)(
        PATIENT, history(nLevels), E.nextC0Time(dayjs().add(1, 'day').startOf('day')));
    return html;
}

(async () => {
    section('runForecast wiring');
    const ciWidths = [];
    for (const n of [1, 3, 5]) {
        const out = await forecast(n);
        ok(`n=${n}: completes without throwing`,
            !/Forecast error/.test(out['alert-container'] || '') && (out['results-summary'] || '').length > 0);
        const m = (out['results-summary'] || '').match(/<strong>([\d.]+)<\/strong> – <strong>([\d.]+)<\/strong>/);
        if (m) ciWidths.push({ n, w: +m[2] - +m[1] });
    }
    ok('90% CI narrows monotonically as levels accumulate',
        ciWidths.length === 3 && ciWidths[0].w > ciWidths[1].w && ciWidths[1].w > ciWidths[2].w,
        ciWidths.map(c => `n=${c.n}: ${c.w.toFixed(2)}`).join(', '));

    section('Chart data integrity');
    await forecast(3);
    const ds = chart.data.datasets;
    const curve = ds.find(d => d.label === 'Individualized Forecast').data;
    const band = ds.find(d => d.label === '90% CI Lower').data;
    const span = curve[curve.length - 1].x - curve[0].x;
    const spacing = span / (curve.length - 1);

    ok('curve grid resolves the 12 h dosing cycle (<= 0.5 h)', spacing <= 0.5 + 1e-9,
        `${curve.length} pts, ${spacing.toFixed(3)} h spacing over ${span.toFixed(0)} h`);
    ok('CI band keeps its own coarse grid', band.length === 600);
    ok('CI band spans the same window as the curve', Math.abs(band[band.length - 1].x - span) < 1e-6);
    ok('no NaN or Infinity in the plotted curve', curve.every(p => isFinite(p.x) && isFinite(p.y)));
    ok('no NaN or Infinity in the CI band', band.every(p => isFinite(p.x) && isFinite(p.y)));
    ok('concentrations are non-negative', curve.every(p => p.y >= 0));

    section('Single-observation fit is not falsely condemned');
    const out1 = await forecast(1);
    ok('no spurious "DO NOT TRUST FORECAST"', !/DO NOT TRUST FORECAST/.test(out1['results-summary'] || ''));
    ok('no spurious "Severe Model Mismatch"', !/Severe Model Mismatch/.test(out1['model-params-output'] || ''));

    section('Mis-timed sample is surfaced');
    Object.keys(html).forEach(k => delete html[k]);
    const log = history(2);
    const aDose = log.filter(x => x.dose > 0).slice(-10)[0];
    log.push({ id: 'bad', recordDate: aDose.recordDate.add(15, 'minute'), dose: null, level: 12.5, time: aDose.time + 0.25 });
    await vm.runInContext('runForecast', ctx)(
        PATIENT, log.sort((a, b) => a.time - b.time), E.nextC0Time(dayjs().add(1, 'day').startOf('day')));
    ok('post-dose sample raises the C0 timing warning',
        /POSSIBLE POST-DOSE SAMPLE/.test(html['accuracy-metrics-grid'] || ''));

    done('smoke');
})();
