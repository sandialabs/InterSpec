/* electron-builder configuration for the InterSpec desktop app.
 *
 * This replaces the previous `@electron/packager` invocations.  The reason for the change is
 * Linux: electron-packager produces a bare directory, and shipping that as a .zip gave users an
 * app that could not be installed, had no desktop/menu entry or file associations, and usually
 * refused to start at all, because a .zip cannot carry the setuid bit that Chromium requires on
 * `chrome-sandbox`.  electron-builder gives us .deb/.rpm (which install to /opt, register the
 * desktop entry and icons, make chrome-sandbox setuid-root, and install an AppArmor profile),
 * an AppImage, and a plain directory we turn into a tarball.
 *
 * This is a .js rather than a .yml config because the `fpm` arguments that place files into
 * absolute system paths need absolute source paths, and resolving them from `__dirname` is the
 * only way to be independent of the working directory electron-builder happens to run fpm in.
 *
 * Windows and macOS: these paths are for development only.  Released Windows builds come from
 * target/wxWidgets and released macOS builds from target/osx (see .github/workflows/), so the
 * sections below are a convenience for building/trying the Electron variant, not a release path.
 *
 * The app directory itself is produced by CMake (`cmake-js ... --target install`, see
 * target/electron/CMakeLists.txt), and is passed in with
 *   electron-builder --linux -c.directories.app=/path/to/build/app
 */

const path = require('path');

/** Absolute path to a file in this directory, for fpm operands. */
const here = ( ...parts ) => path.join( __dirname, ...parts );

/* Files the package manager should place outside the app directory.  fpm takes these as
 * `source=destination` operands.
 *
 * Deliberately done this way rather than through electron-builder's `afterInstall` hook:
 * setting `afterInstall` *replaces* electron-builder's own post-install script, which is what
 * makes chrome-sandbox setuid-root and installs the AppArmor profile - losing that would undo
 * the main reason for packaging properly in the first place.  Letting the package manager own
 * these files is also better behaved: both dpkg and rpm have triggers on
 * /usr/share/mime/packages and /usr/share/applications that re-run update-mime-database and
 * update-desktop-database automatically, so no script is needed.
 */
const systemFiles = [
  here('linux', 'gov.sandia.InterSpec.xml') + '=/usr/share/mime/packages/gov.sandia.InterSpec.xml',
  here('linux', 'gov.sandia.InterSpec.metainfo.xml') + '=/usr/share/metainfo/gov.sandia.InterSpec.metainfo.xml',
];

/* Runtime libraries Electron needs that are not guaranteed to be installed.  libstdc++6 and
 * libgcc-s1 are on the list because InterSpec links the C++ runtime dynamically - see
 * cmake/CxxRuntimePolicy.cmake and https://github.com/sandialabs/InterSpec/issues/51.
 */
const debDepends = [
  'libstdc++6', 'libgcc-s1',
  'libgtk-3-0', 'libnotify4', 'libnss3', 'libxss1', 'libxtst6', 'xdg-utils',
  'libatspi2.0-0', 'libuuid1', 'libdrm2', 'libgbm1', 'libcups2', 'libxkbcommon0',
  'libasound2t64 | libasound2',
];

const rpmDepends = [
  'libstdc++', 'libgcc',
  'gtk3', 'libnotify', 'nss', 'libXScrnSaver', 'libXtst', 'xdg-utils',
  'at-spi2-atk', 'libuuid', 'libdrm', 'mesa-libgbm', 'cups-libs', 'libxkbcommon',
  'alsa-lib',
];

module.exports = {
  appId: 'gov.sandia.InterSpec',
  // Must not contain a space: electron-builder generates a syntactically invalid AppArmor
  // profile from a productName with whitespace (electron-userland/electron-builder#9230).
  productName: 'InterSpec',
  copyright: 'Copyright National Technology & Engineering Solutions of Sandia, LLC',

  /* The C++ code opens data/, InterSpec_resources/, example_spectra/ and Wt's resources/ as
   * ordinary filesystem paths, relative to the directory main.js lives in (see
   * InterSpecServer::start_server).  None of that can be read from inside an asar archive, so
   * asar has to stay off.  This also keeps the on-disk layout identical to what
   * electron-packager produced, which is one less thing to re-verify.
   */
  asar: false,

  // The native addon is built by cmake-js before electron-builder runs, and the app has no npm
  // dependencies, so there is nothing for electron-builder to rebuild.
  npmRebuild: false,
  nodeGypRebuild: false,

  /* Never let electron-builder publish.  When it detects CI it otherwise turns publishing on by
   * itself ("Implicit publishing triggered by CI detection"), infers a GitHub provider from the
   * repository URL, and then fails the build with
   *   ⨯ GitHub Personal Access Token is not set, neither programmatically, nor using env "GH_TOKEN"
   * - which is exactly what happened on the first green Linux build, *after* all the packages had
   * been produced.  It does not reproduce locally, because `CI` is not set there.
   *
   * Release assets are uploaded by the `release` job in .github/workflows/build_app.yml, not from
   * here.  This also stops the pointless latest-linux.yml auto-update metadata being generated.
   */
  publish: null,

  directories: {
    // `app` is supplied on the command line, pointing at the CMake-produced app directory.
    output: 'dist',
  },

  // Keep the InterSpec_app* prefix: the release workflow deletes previous assets by that
  // pattern before uploading new ones.
  artifactName: 'InterSpec_app_${os}_${arch}_${version}.${ext}',

  // Deep links.  Becomes x-scheme-handler/... entries in the .desktop MimeType list on Linux,
  // and CFBundleURLTypes on macOS.
  protocols: [
    {
      name: 'InterSpec',
      schemes: [ 'interspec', 'raddata' ],
    },
  ],

  /* Spectrum files.  The extension list matches the macOS declaration in macOS/Info.plist so
   * the two platforms agree.  Note that on Linux this only adds the MIME type to the .desktop
   * file - the type itself is defined by linux/gov.sandia.InterSpec.xml, installed above.
   */
  fileAssociations: [
    {
      ext: [ 'n42', 'pcf', 'cnf', 'spe', 'spc', 'mca', 'mps', 'lis', 'chn' ],
      name: 'Gamma energy spectrum file',
      description: 'Gamma energy spectrum file',
      mimeType: 'application/gamma-spectrum',
      role: 'Viewer',
    },
  ],

  linux: {
    /* Sit alongside the executable, so they are easy to find in an extracted tarball or in
     * /opt/InterSpec.  Scoped per-platform rather than set globally, because on macOS
     * `extraFiles` writes inside Contents/, which invalidates a code signature.
     * The CycloneDX SBOM is copied into the CMake app directory by CI, so it rides along with
     * the rest of the app payload and a local build without one still works.
     */
    extraFiles: [
      { from: '../../LICENSE.txt', to: 'LICENSE.txt' },
      { from: '../../NOTICE.html', to: 'NOTICE.html' },
    ],

    executableName: 'interspec',
    category: 'Science',
    synopsis: 'Spectral radiation analysis software',
    description: 'InterSpec assists in analyzing spectral nuclear radiation data, '
                 + 'using a peak-based methodology.',
    maintainer: 'InterSpec <InterSpec@sandia.gov>',
    vendor: 'Sandia National Laboratories',
    icon: 'linux/icons',
    // Note the `${arch}` a package format substitutes is its own (deb says amd64, rpm x86_64),
    // so the workflow globs on `InterSpec_app_*` rather than assuming a spelling.
    artifactName: 'InterSpec_app_Linux_${arch}_${version}.${ext}',
    // `dir` is here so CI can add the sandbox-probing launcher script and build the tarball
    // itself; see .github/workflows/build_app.yml.
    target: [ 'deb', 'rpm', 'AppImage', 'dir' ],
    // Keeps the installed .desktop filename, its StartupWMClass, and Electron's app_id all
    // equal to `desktopName` from the app's package.json.  Without that, desktop environments
    // cannot associate a running window with the launcher, and show a generic icon and name in
    // the dock and alt-tab.
    syncDesktopName: true,
    desktop: {
      entry: {
        Keywords: 'gamma;spectroscopy;radiation;nuclide;spectrum;isotope;',
      },
    },
  },

  /* `afterInstall` REPLACES electron-builder's own post-install template rather than appending to
   * it (app-builder-lib/out/targets/FpmTarget.js), so linux/fpm-after-install.sh is a fork of that
   * template.  It differs in exactly one place - how it decides whether chrome-sandbox needs to be
   * setuid-root - because upstream probes `unshare --user true` while running as root at install
   * time, and the Ubuntu 23.10+ userns restriction exempts root.  See the comment in that file.
   *
   * The "Check electron-builder postinst template drift" CI step fails the build when upstream's
   * template changes, so this fork gets re-synced deliberately instead of quietly rotting.
   */
  deb: {
    depends: debDepends,
    fpm: systemFiles,
    afterInstall: 'linux/fpm-after-install.sh',
  },

  rpm: {
    depends: rpmDepends,
    fpm: systemFiles,
    afterInstall: 'linux/fpm-after-install.sh',
  },

  appImage: {
    // AppImages run from a nosuid FUSE mount, so chrome-sandbox can never be setuid-root there;
    // electron-builder's AppRun passes --no-sandbox for this reason.
    license: '../../LICENSE.txt',
  },

  win: {
    // Development convenience only - released Windows builds come from target/wxWidgets.
    target: [ 'dir', 'zip' ],
    icon: 'windows/icon.ico',
    legalTrademarks: 'Sandia National Laboratories',
    extraFiles: [
      { from: '../../LICENSE.txt', to: 'LICENSE.txt' },
      { from: '../../NOTICE.html', to: 'NOTICE.html' },
    ],
  },

  mac: {
    /* Development convenience only - released macOS builds come from target/osx.
     *
     * Signing and notarization use electron-builder's standard environment variables
     * (CSC_NAME or CSC_LINK/CSC_KEY_PASSWORD, and APPLE_KEYCHAIN/APPLE_KEYCHAIN_PROFILE or
     * APPLE_ID/APPLE_APP_SPECIFIC_PASSWORD/APPLE_TEAM_ID).  macOS/codesign-electron.sh remains
     * available for signing a bundle out of band.
     */
    target: [ 'dir', 'dmg' ],
    category: 'public.app-category.developer-tools',
    icon: 'macOS/InterSpec.icns',
    hardenedRuntime: true,
    entitlements: 'macOS/entitlements.mac.plist',
    entitlementsInherit: 'macOS/entitlements.mac.plist',
    extendInfo: {
      // Kept from macOS/Info.plist; the document- and URL-type entries are generated by
      // electron-builder from `fileAssociations` and `protocols` above.
      NSHighResolutionCapable: true,
      NSSupportsAutomaticGraphicsSwitching: true,
    },
  },
};
