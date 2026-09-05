# marked-katex-extension 5.1.12 (vendored into InterSpec)

- **Upstream:** <https://github.com/UziTech/marked-katex-extension>
- **Release:** <https://github.com/UziTech/marked-katex-extension/releases/tag/v5.1.12>
- **Version:** 5.1.12
- **Downloaded:** 20260904, by wcjohns

This is the glue that lets `marked` hand `$...$` / `$$...$$` spans to KaTeX.  The GitHub release
archive contains only `src/index.ts` - no browser build - so the file here was taken from the npm
tarball of the same release (`npm pack marked-katex-extension@5.1.12`), whose `LICENSE` and
`src/index.ts` are byte-identical to those in the GitHub source archive.

Requires `katex` to already be defined as a global when it loads, and `marked` to be loaded before it
is used; `InterSpec_resources/LlmMarkdown.js` and its callers load the three in that order.

## Files kept (unmodified)

| File | From |
|------|------|
| `index.umd.js` | `lib/index.umd.js` |
| `LICENSE`      | `LICENSE`          |

## Files deleted

Everything else: `src/index.ts`, the CJS and ESM builds (`lib/index.cjs`, `lib/index.esm.js`),
`package.json`, docs and tests.

**No modifications were made to any kept file.**

## Local-content verification (20260904)

`index.umd.js` is about 90 lines: two regular expressions, a `marked` tokenizer/renderer pair, and a
call to `katex.renderToString()`.  It contains no `XMLHttpRequest`, `fetch`, `WebSocket`,
`EventSource`, `sendBeacon`, `Image`, `Worker`, `importScripts`, `document.cookie`,
`localStorage`/`sessionStorage`/`indexedDB`, `eval(`, `new Function` or `document.write`, and no
absolute URLs of any kind.

## Copyright / license

Copyright (c) 2021 @markedjs (Tony Brix).

marked-katex-extension is MIT licensed; see `LICENSE`.

## Known upstream behaviour worth being aware of

The `$...$` inline rule's content class is `[^\\\n]`, which permits a `$` *inside* the span; only the
opening context (preceded by a space) and the closing context (followed by whitespace or punctuation)
constrain the match.  So a sentence that mixes currency with math, such as

    It cost $5 and $10, while $\chi^2/\mathrm{DOF}$ was 1.04.

opens the span at the currency `$` and closes it after `DOF}`, and KaTeX then shows a red parse error
for that stretch.  Plain currency on its own ("It cost $5 to make and $10 to ship.") is unaffected,
because neither `$` pair satisfies the closing lookahead.

This is upstream behaviour, shared with essentially every `$`-delimited Markdown math renderer, and is
left as-is: the failure is loudly visible rather than silent, the response's "Show Raw Text" menu item
shows the source, and currency figures are vanishingly rare in this application's replies.

Note that `InterSpec_resources/LlmMarkdown.js` registers one supplementary `$...$` rule of its own, for
the case upstream deliberately skips: a span hugged by a letter, as in `$^{235}$U`.  It only fires when
the span contains a LaTeX control character (`\ ^ _ { }`), so currency amounts never match it.
