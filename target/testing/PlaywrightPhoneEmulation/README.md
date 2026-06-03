# Playwright phone emulation for InterSpec

This directory contains a tiny Playwright script that launches a Chromium
window sized and styled as a chosen mobile device, then points it at a running
InterSpec server with the `?isphone=1` (or `?istablet=1`) URL parameter.

It is a more honest "phone-on-a-desktop" simulator than Chrome DevTools'
device-toolbar mode: the viewport actually matches the device size,
`window.innerWidth` reads correctly, and the User-Agent / device-pixel-ratio /
`navigator.maxTouchPoints` / `isMobile` are set together via Playwright's
built-in device descriptors.

## One-time setup

```bash
cd target/testing/PlaywrightPhoneEmulation
npm install                       # installs playwright into ./node_modules
npx playwright install chromium   # downloads the headless Chromium binary
```

(Or install Playwright globally with `npm install -g playwright` and skip the
local install.)

## Running InterSpec

In one terminal, start a local InterSpec server **from the `build_vscode`
directory** so the docroot has the `external_libs/SpecUtils/d3_resources/...`
symlink alongside `InterSpec_resources/`:

```bash
cd build_vscode
./InterSpec --docroot . --http-address 127.0.0.1 --http-port 8080 \
    -c ./data/config/wt_config_web.xml
```

(If you don't have the `external_libs` symlink, set up: `ln -s
../../external_libs/SpecUtils/d3_resources build_vscode/external_libs/SpecUtils/d3_resources`.)

## Launching the phone window

In another terminal:

```bash
cd target/testing/PlaywrightPhoneEmulation
node phone-test.js                              # iPhone 14 portrait
node phone-test.js --device "Pixel 7"
node phone-test.js --device "iPad Mini" --param istablet
node phone-test.js --url http://127.0.0.1:8080/?v=2  # extra query for cache-bust
node phone-test.js --headless                   # no window; for CI / scripted use
```

The window stays open until you close it. Inside the window you can:

- Hard-reload (Cmd-Shift-R / Ctrl-Shift-R) to pick up CSS / XML changes
  without restarting the InterSpec server. *Note:* InterSpec's `@import`
  chain doesn't cache-bust on its own; for stubborn CSS changes also try
  appending `&v=N` to the URL to force a different session token.
- Open DevTools (Cmd-Opt-I / Ctrl-Shift-I) for the usual inspector / console.

Closing the browser window terminates the script.

## Common device options

| Preset | viewport | DPR | UA family |
|---|---|---|---|
| `iPhone 14` (default) | 390×844 | 3 | iOS Safari |
| `iPhone SE` | 320×568 | 2 | iOS Safari |
| `iPhone 14 Pro Max` | 430×932 | 3 | iOS Safari |
| `iPad Mini` | 768×1024 | 2 | iPadOS Safari |
| `Pixel 7` | 412×915 | 2.625 | Android Chrome |
| `Galaxy S9+` | 320×658 | 4.5 | Android Chrome |

The complete list (~120 entries, including landscape variants):

```bash
node -e "console.log(Object.keys(require('playwright').devices).join('\n'))"
```

## What does `?isphone=1` do?

It forces `InterSpecApp::isPhone()` to return `true` on any build, not just
iOS. (`?istablet=1` is the equivalent for tablet mode.) `InterSpecApp.cpp`
checks the URL parameter at session startup; the rest of the code branches on
`isPhone()` / `isTablet()` / `isMobile()` to switch layout, menus, AuxWindow
behavior, and CSS classes.

Combined with Playwright's device descriptor (which also makes the page see
a mobile UA), this gives the closest possible match to a real phone session.

## Files in this directory

- `phone-test.js` — the launcher script. Read the header comment for full options.
- `package.json` — pins the Playwright dependency.
- `README.md` — this file.

The directory is intentionally outside the main CMake build — Playwright /
Node is a dev-tooling dependency that desktop builds don't need.
