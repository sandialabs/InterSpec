#!/bin/bash
#
# Sign and notarize a packaged InterSpec.app on macOS using @electron/osx-sign
# and @electron/notarize. Replaces the legacy hand-rolled codesign loop, which
# stopped tracking the changing list of frameworks Electron ships.
#
# Usage:
#   codesign-electron.sh path/to/InterSpec.app
#
# Required setup (one-time):
#   xcrun notarytool store-credentials "interspec-notary" \
#     --apple-id <your-apple-id> \
#     --team-id <YOUR_TEAM_ID> \
#     --password <app-specific-password>
#
#   This stores the notarization credential in the macOS login keychain
#   under the profile name "interspec-notary".  No Apple-ID password is
#   ever passed via env vars, command-line, or shell history after that.
#
# Required environment variables:
#   OSX_SIGN_IDENTITY              Developer ID Application identity
#                                  (e.g. "Developer ID Application: ...")
#
# Optional environment variables:
#   NOTARY_KEYCHAIN_PROFILE        Override the keychain profile name
#                                  (default: "interspec-notary")
#
# The `npm run package-mac` script invokes signing+notarization inline via
# @electron/packager flags using the same keychain profile; this script is
# provided for re-signing an existing .app bundle out of band (e.g., after
# a manual edit).

set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 path/to/InterSpec.app" >&2
  exit 1
fi

app="$1"

: "${OSX_SIGN_IDENTITY:?OSX_SIGN_IDENTITY env var is required}"
notary_profile="${NOTARY_KEYCHAIN_PROFILE:-interspec-notary}"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
entitlements="${script_dir}/entitlements.mac.plist"

echo "Signing ${app} with @electron/osx-sign..."
npx --yes @electron/osx-sign "${app}" \
  --identity="${OSX_SIGN_IDENTITY}" \
  --hardened-runtime \
  --entitlements="${entitlements}" \
  --entitlements-inherit="${entitlements}" \
  --gatekeeper-assess=false

echo "Verifying codesign..."
codesign --verify --deep --display --verbose=4 "${app}"

echo "Notarizing ${app} with @electron/notarize (keychain profile: ${notary_profile})..."
npx --yes @electron/notarize \
  --app-path="${app}" \
  --keychain-profile="${notary_profile}"

echo "Verifying Gatekeeper acceptance..."
spctl --ignore-cache --no-cache --assess --type execute --verbose=4 "${app}"

echo "Done."
