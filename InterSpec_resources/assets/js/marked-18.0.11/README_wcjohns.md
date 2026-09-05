# marked 18.0.11 (vendored into InterSpec)

- **Upstream:** <https://github.com/markedjs/marked>
- **Release:** <https://github.com/markedjs/marked/releases/tag/v18.0.11>
- **Version:** 18.0.11
- **Downloaded:** 20260904, by wcjohns

The GitHub release archive contains only TypeScript sources - no browser build - so the file here was
taken from the npm tarball of the same release (`npm pack marked@18.0.11`), whose `LICENSE` is
byte-identical to the one in the GitHub source archive.

## Files kept (unmodified)

| File | From |
|------|------|
| `marked.umd.js` | `lib/marked.umd.js` |
| `LICENSE`       | `LICENSE`          |

## Files deleted

Everything else: TypeScript sources, the ESM build (`lib/marked.esm.js`), type declarations
(`lib/marked.d.ts`), source maps (`*.js.map`), the CLI (`bin/`, `man/`), `package.json`, docs and
tests.

**No modifications were made to any kept file.**

Note that `marked.umd.js` ends with a `//# sourceMappingURL=marked.umd.js.map` comment and the map
file was not kept.  This is harmless: browsers only request a source map when developer tools are
open, and the request would simply 404 against the local docroot.

## Local-content verification (20260904)

`marked.umd.js` contains no `XMLHttpRequest`, `fetch`, `WebSocket`, `EventSource`, `sendBeacon`,
`Image`, `Worker`, `importScripts`, `document.cookie`, `localStorage`/`sessionStorage`/`indexedDB`,
`eval(`, `new Function` or `document.write`.  The only absolute URLs present are inert text: the
project's own `https://github.com/markedjs/marked` in the file banner and in an error message, and
the string `"http://"` used to prefix bare `www.` autolinks.  Nothing is ever fetched.

InterSpec additionally overrides marked's `link` and `image` renderers (see
`InterSpec_resources/LlmMarkdown.js`) so that no `<a>` or `<img>` element is ever emitted from LLM
output in the first place.

## Copyright / license

Copyright (c) 2018+, MarkedJS (<https://github.com/markedjs/>)
Copyright (c) 2011-2018, Christopher Jeffrey (<https://github.com/chjj/>)

marked itself is MIT licensed.  Note that the shipped `LICENSE` carries a second block as well: the
BSD-3-Clause-style license of John Gruber's original Markdown (Copyright (c) 2004, John Gruber),
which includes a no-endorsement clause on the name "Markdown".  Upstream's `package.json` declares
the package as MIT; the full `LICENSE` is shipped verbatim, so both grants travel with the code.
