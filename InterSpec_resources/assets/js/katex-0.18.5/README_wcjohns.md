# KaTeX 0.18.5 (vendored into InterSpec)

- **Upstream:** <https://katex.org/> / <https://github.com/KaTeX/KaTeX>
- **Release:** <https://github.com/KaTeX/KaTeX/releases/tag/v0.18.5>
- **Version:** 0.18.5
- **Downloaded:** 20260904, by wcjohns

The GitHub *source* archive contains only TypeScript sources (plus the fonts) - no browser build - so
the files here were taken from the npm tarball of the same release (`npm pack katex@0.18.5`), which
has the same content as the release's prebuilt `katex.zip` asset.  Verified that the npm tarball's
`LICENSE` and all twenty `.woff2` font files are byte-identical to those in the GitHub source archive.

## Files kept (unmodified)

| File | From |
|------|------|
| `katex.min.js`      | `dist/katex.min.js`      |
| `katex.min.css`     | `dist/katex.min.css`     |
| `fonts/*.woff2` (20 files) | `dist/fonts/*.woff2` |
| `LICENSE`           | `LICENSE`                |

`fonts/` must stay in a subdirectory next to `katex.min.css`, because the stylesheet's `@font-face`
rules reference the fonts as `url(fonts/KaTeX_*)`, relative to the stylesheet.

## Files deleted

- TypeScript sources (`src/`, `katex.ts`, `types/`), docs, tests, `package.json`, `cli.js`.
- Non-minified and ESM builds: `dist/katex.js`, `dist/katex.mjs`, `dist/katex.css`,
  `dist/katex-swap.css`, `dist/katex-swap.min.css`.
- All `contrib/` extensions and their builds (`auto-render`, `mhchem`, `copy-tex`,
  `mathtex-script-type`, `render-a11y-string`) - InterSpec drives KaTeX through
  `marked-katex-extension`, not through auto-render.
- The `.woff` and `.ttf` variants of every font.  Each `@font-face` rule lists `woff2` first with a
  `format()` hint, and every browser InterSpec supports (Safari 16.4+, Chromium-based Edge/Chrome,
  Android WebView) reads woff2, so the other two variants are never requested.  Dropping them saves
  about 2.4 MB.

**No modifications were made to any kept file.**

## Local-content verification (20260904)

- `katex.min.css`: every `url(...)` is a relative `url(fonts/KaTeX_*.{woff2,woff,ttf})`.  There are no
  `http`/`https` URLs and no `@import` of anything.
- `katex.min.js`: contains no `XMLHttpRequest`, `fetch` (the 28 occurrences of `fetch(` are all
  KaTeX's own `Parser.prototype.fetch()` token method, not `window.fetch`), `WebSocket`,
  `EventSource`, `sendBeacon`, `Image`, `Worker`, `importScripts`, `document.cookie`,
  `localStorage`/`sessionStorage`/`indexedDB`, `eval(`, `new Function` or `document.write`.  The only
  absolute URLs present are the XML namespace identifiers `http://www.w3.org/1998/Math/MathML` and
  `http://www.w3.org/2000/svg`, which are namespace names, never fetched.

InterSpec invokes KaTeX with its default `trust: false`, which disables the `\href`, `\url` and
`\includegraphics` commands - so LaTeX from the LLM cannot reference an external resource either.

## Copyright / license

Copyright (c) 2013-2020 Khan Academy and other contributors.

KaTeX is MIT licensed; see `LICENSE`.  The fonts are derived from the American Mathematical Society's
Computer Modern fonts and are covered by the same MIT license (upstream ships them under it).
