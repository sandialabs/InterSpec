# Theme contrast check

`theme-contrast.js` drives a running InterSpec (default `http://127.0.0.1:8080/`) with
Playwright, switches between the light and dark colour themes through the app's own
`prefers-color-scheme` hook, and measures the WCAG contrast of representative text, controls and
borders in each.  It exits non-zero when anything is below AA (4.5:1 for text, 3:1 for borders
and large text).

The theme switch relies on the `AutoDarkFromOs` preference (on by default) and on the "Default"
colour theme being selected; the script checks the `--interspec-color-scheme` token before
measuring and explains what to reset if the expected theme is not active.

```
# one-time: share the phone harness' Playwright install
cd ../PlaywrightPhoneEmulation && npm install && npx playwright install chromium

# then, with InterSpec running
cd ../PlaywrightThemeContrast
NODE_PATH=../PlaywrightPhoneEmulation/node_modules node theme-contrast.js [--url URL] [--scheme light|dark|both] [--headed] [--screenshots DIR]
```

The token values live in `InterSpec_resources/themes/default/default.css` (light) and
`InterSpec_resources/themes/dark/dark.css`; tune a failing token there and re-run.
