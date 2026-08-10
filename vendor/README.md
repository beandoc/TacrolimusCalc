# Vendored frontend dependencies

These are served from disk rather than a CDN. Three reasons:

1. **file:// operation.** The app is designed to be opened straight off disk with
   no server (see the `protocol === "file:"` branches in `index.html` and the
   note at the top of `pk-engine.js`). CDN `<script>` tags silently defeated
   that — offline, the page lost Chart.js, flatpickr and dayjs, and dayjs is a
   hard requirement of the PK engine.
2. **It is a dosing tool.** A compromised or hijacked CDN could change what it
   recommends. Two of the old tags were also unpinned (`npm/flatpickr`,
   `dayjs@1`), so the shipped behaviour could change with no commit in this repo.
3. Hospital networks routinely block or MITM public CDNs.

## Contents

| File | Version | Source |
|---|---|---|
| `chart.umd.min.js` | 4.4.0 | `https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js` |
| `flatpickr.min.js` | 4.6.13 | `https://cdn.jsdelivr.net/npm/flatpickr@4.6.13/dist/flatpickr.min.js` |
| `flatpickr.min.css` | 4.6.13 | `https://cdn.jsdelivr.net/npm/flatpickr@4.6.13/dist/flatpickr.min.css` |
| `dayjs.min.js` | 1.11.13 | `https://cdn.jsdelivr.net/npm/dayjs@1.11.13/dayjs.min.js` |
| `dayjs-customParseFormat.js` | 1.11.13 | `https://cdn.jsdelivr.net/npm/dayjs@1.11.13/plugin/customParseFormat.js` |

Keep `dayjs` here in step with the `dayjs` devDependency in `package.json`: the
browser uses this copy, the tests use the npm one, and they must agree about
date parsing or the tests stop describing the app.

## Not vendored

- **Tailwind** — removed entirely. It was the dev-only in-browser JIT build
  (~400 KB; its own documentation says not to use it in production) and the whole
  page used exactly one Tailwind class, `text-white`. That is now a plain rule in
  the `<style>` block, alongside `.hidden`, which `switchTab()` relies on.
- **Google Fonts** — still fetched, loaded non-blockingly. Purely cosmetic: it
  cannot change a number on the page, and every `font-family` names concrete
  fallbacks, so an offline machine renders immediately in the fallback stack.

## Refreshing

```sh
curl -sSfL -o vendor/chart.umd.min.js  "https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"
curl -sSfL -o vendor/flatpickr.min.js  "https://cdn.jsdelivr.net/npm/flatpickr@4.6.13/dist/flatpickr.min.js"
curl -sSfL -o vendor/flatpickr.min.css "https://cdn.jsdelivr.net/npm/flatpickr@4.6.13/dist/flatpickr.min.css"
curl -sSfL -o vendor/dayjs.min.js      "https://cdn.jsdelivr.net/npm/dayjs@1.11.13/dayjs.min.js"
curl -sSfL -o vendor/dayjs-customParseFormat.js "https://cdn.jsdelivr.net/npm/dayjs@1.11.13/plugin/customParseFormat.js"
```

Then re-run `npm run verify` and load the page over `file://` to confirm the
chart and date pickers still render.
