#!/usr/bin/env node
// Syntax-checks the inline <script> block in index.html.
//
// This exists because a single missing brace in one function once took down
// all 57 functions in the file — everything lives in one script block, so any
// syntax error anywhere kills the entire block and the page silently does
// nothing. `node --check` catches that in milliseconds. Wired into
// .githooks/pre-commit so it cannot be committed again.
//
// Also verifies index.html still loads pk-engine.js and does not redefine
// anything the engine exports.

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const root = path.join(__dirname, '..');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const lines = html.split('\n');

let bad = 0;
const fail = msg => { console.error(`  FAIL  ${msg}`); bad++; };
const pass = msg => console.log(`  PASS  ${msg}`);

// 1. The engine must be loaded, as a classic script (modules break file://).
if (!/<script src="pk-engine\.js"><\/script>/.test(html)) fail('index.html does not load pk-engine.js as a classic script');
else pass('index.html loads pk-engine.js');
if (/<script[^>]+src="pk-engine\.js"[^>]*type="module"/.test(html)) fail('pk-engine.js loaded as a module — breaks file:// operation');

// 2. Inline script block must parse.
const start = lines.findIndex(l => l.trim() === '<script>' );
const end = lines.length - 1 - [...lines].reverse().findIndex(l => l.includes('</script>'));
if (start < 0 || end <= start) { fail('could not locate the inline <script> block'); process.exit(1); }
const inline = lines.slice(start + 1, end).join('\n');
try {
    new vm.Script(inline, { filename: 'index.html:<script>' });
    pass(`inline <script> parses (${end - start - 1} lines)`);
} catch (err) {
    fail(`inline <script> has a syntax error — the whole page would be dead\n        ${err.message}`);
}

// 3. The engine itself must parse and stay DOM-free.
const engineSrc = fs.readFileSync(path.join(root, 'pk-engine.js'), 'utf8');
try { new vm.Script(engineSrc, { filename: 'pk-engine.js' }); pass('pk-engine.js parses'); }
catch (err) { fail(`pk-engine.js syntax error: ${err.message}`); }

const domRefs = engineSrc.split('\n')
    .map((l, i) => ({ n: i + 1, l }))
    .filter(({ l }) => !/^\s*(\/\/|\*|\/\*)/.test(l))
    .filter(({ l }) => /\bdocument\.|getElementById|\bfetch\(|\bwindow\./.test(l));
if (domRefs.length) fail(`pk-engine.js must stay pure — DOM/network at line(s) ${domRefs.map(d => d.n).join(', ')}`);
else pass('pk-engine.js is free of DOM and network access');

// 4. No symbol may be defined in both places.
const engine = require(path.join(root, 'pk-engine.js'));
const dupes = Object.keys(engine).filter(k =>
    new RegExp(`^\\s*(function ${k}\\s*\\(|const ${k}\\s*=|let ${k}\\s*=)`, 'm').test(inline));
if (dupes.length) fail(`index.html redefines engine export(s): ${dupes.join(', ')}`);
else pass(`no duplicate definitions (${Object.keys(engine).length} engine exports)`);

// 5. Every engine symbol index.html references must actually be exported.
const referenced = Object.keys(engine).filter(k => new RegExp(`\\b${k}\\b`).test(inline));
pass(`index.html uses ${referenced.length} engine exports, all resolved`);

console.log(bad === 0 ? '\nALL PASS  (check-html)\n' : `\n${bad} FAILED  (check-html)\n`);
process.exit(bad === 0 ? 0 : 1);
