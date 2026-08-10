# Creating a Distributable Electron Package

In order to create the [Electron]([https://electronjs.org/](https://electronjs.org/) packaged version of `InterSpec`, we create a node add-on that contains all the C++ code, and then request to start the `InterSpec` server from [main.js](app/main.js).  To create the add-on we use [cmake-js](https://www.npmjs.com/package/cmake-js) to build things, and [node-addon-api](https://www.npmjs.com/package/node-addon-api) to actually interface the C++ to JS.


You can either manually build the InterSpec dependencies (boost, Wt, zlib, and for macOS libpng and libharu), or you can have CMake fetch and build these dependencies.  Having CMake fetch and build the dependencies takes maybe an hour to build things the first time (however, only a couple minutes of this requires your attention) and then subsequent builds are a couple minutes; if you plan to make a substantial number of changes and re-compilations, then the manually built dependencies will yield a slightly faster compile time after each change.

## Building with manually compiled dependencies
To compile the InterSpec code, and package the Electron app, with the manually compiled dependencies (see [patches/README.md](../patches/README.md) for instructions), the following commands are a good base to start with:

```bash
npm install -g cmake-js 

# For macOS only, you may want to define a deployment target
export MACOSX_DEPLOYMENT_TARGET=11.0

cd /path/to/InterSpec/target/electron

# Install dependency for compiling a node.js add-on
npm install --save-dev node-addon-api --arch=x64

# If boost and Wt are in standard locations, you can just run
cmake-js

# Or to have a little more control over things
cmake-js --generator "Visual Studio 17 2022" \
         --architecture x64 --arch=x64 \
         --CDCMAKE_PREFIX_PATH=C://Path/To/Wt_3.7.1_prefix \
         --CDBoost_USE_STATIC_RUNTIME=ON \
         --CDCMAKE_BUILD_TYPE="Release" \
         --CDLEAFLET_MAPS_KEY="..." \
         --out="build_win" \
         --target install

# If you make changes and want to recompile
CMAKE_BUILD_PARALLEL_LEVEL=12 cmake-js build --out="build_dir"
# Or you can use CMake to run the `make` command, which can be useful when
# the CMake generator isnt a command-line based system like (ex Xcode, MSVC)
cmake --build build_dir --config Release
# Or directly use the `make` command like:
ninja -C build_dir

# Make the "app" directory into the build directory
ninja -C build_dir install

# And also copy all the InterSpec resources to the 
# 'app' sub-directory of your build dir
cmake-js build --out=build_dir --target install
# Or
cmake --build build_dir --target install --config Release --clean-first
# Or
ninja -C build_dir install


# To run InterSpec without packaging everything, you
# need to install the Electron package
npm install electron --arch=x64

# Then to actually run things
npm start


# Packaging is done with electron-builder, configured in electron-builder.js.
# It is already in devDependencies, so `npm install` above has it.

# Then to actually create the distributable package, run one of the following
# (assuming your build directory is 'build_macos', 'build_win', or 'build_linux',
#  otherwise pass -c.directories.app=/your/build/app)
# Note: package-mac produces an arm64 (Apple Silicon) build; Intel macOS is no
# longer supported.
npm run package-mac
npm run package-win
npm run package-linux

# Results land in the 'dist' directory.
```

On Linux this produces a `.deb`, `.rpm`, `.AppImage`, and an unpacked `linux-unpacked/`
directory. `LICENSE.txt` and `NOTICE.html` are copied in automatically.

The runtime libraries Electron needs (`libnss3 libxss1 libasound2 libdrm2 libgbm1 libgtk-3-0
libxkbcommon0 libatk-bridge2.0-0 libcups2`, plus `libstdc++6`/`libgcc-s1`) are declared as
dependencies of the `.deb`/`.rpm`, so the package manager installs them. Users of the AppImage or
tarball on a minimal system may still need to install them by hand.

If you are cross-compiling, you can, for example, build the Linux package from macOS using a command like

```bash
npm --target_arch=x64 --target_platform=linux run package-linux
```

Note: cross-compiling to macOS arm64 from a Linux host is not supported here, as our native addon needs an arm64 toolchain in addition to whatever electron-builder would do. Building Linux `.deb`/`.rpm` packages also needs `dpkg-deb`/`rpmbuild` on the host, so those targets cannot be produced from macOS - use the Docker route below.


## Building with CMake fetched and compiled dependencies
The following commands will compile and package the InterSpec code, starting from the Fedora 35 Docker image; there is nothing special about Fedora, and any of the distributions compatible with npm, node.js, and Electron should work.

Using the CMake FetchContent will fail on macOS because it currently doesnt fetch/build libpng and libharu, but otherwise it seems to work well on Windows, and various flavors of Linux.

```bash
# From host OS terminal - grab InterSpec source code
git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec
cd InterSpec

# Start Fedora Docker container
docker run --rm -it -v `pwd`:/interspec fedora:35 sh

# Install the packages we need
dnf install -y make automake gcc gcc-c++ kernel-devel cmake git patch npm

# Change to directory we need
cd /interspec/target/electron/

# Install the NPM dependencies we need
npm install -g cmake-js
npm install --save-dev node-addon-api --arch=x64
npm install electron --arch=x64

# Configure and build the node-add-on shared library
#  This step may take half and hour or more to clone into the Boost and Wt github repos, and 
#  then another half hour or hour to build
CMAKE_BUILD_PARALLEL_LEVEL=8 cmake-js --architecture x64 --arch=x64 --CDInterSpec_FETCH_DEPENDENCIES=ON --out=build_linux --target install

# Package the code up
npm run package-linux

# If you are using Docker, you will need to install some additional libraries to run things, although the GUI still probably wont run
#   dnf install nss atk at-spi2-atk cups libdrm gtk3

# Run the executable
cd dist/linux-unpacked/
./interspec
```



## Using Docker to build Electron-based InterSpec package
These are the instructions for building the Electron-based InterSpec package for Linux, using the Python Many Linux docker image.

The work is split in three, which is also what CI does (see `.github/workflows/build_app.yml`):

1. **Build a dependency prefix inside the manylinux container.** boost, zlib, Wt, Eigen and Ceres
   are compiled once and installed into a prefix directory. It must be built in the same container
   the app is compiled in, because everything in it is linked into `InterSpecAddOn.node` and so has
   to come from the same gcc-toolset / glibc-2.28 world.
2. **Compile inside the manylinux container**, against that prefix. The image exists only because
   we need to build against an old glibc; it is not a good place to build Debian/RPM packages.
3. **Package on a normal Linux host**, where `dpkg`, `rpm`, `fakeroot` and a current Node are easy
   to get. electron-builder downloads the prebuilt Electron runtime itself, so the container does
   not need it.

Step 1 is optional: pass three arguments instead of four to step 2 and CMake's `FetchContent` will
download and build the dependencies into the build directory instead. That needs no prior setup but
adds roughly 25 minutes to *every* build, whereas a prefix is built once and reused. CI caches the
prefix, which is why it caches nothing else about the C++ build — see the comments in
`build_app.yml` for why caching the *build directory* turned out to be worthless.

### Building using script with manylinux container
```bash
cd /tmp
mkdir build_interspec
cd build_interspec

git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_code
mkdir build_electron
mkdir build_working_dir
mkdir deps_build
mkdir deps_prefix

# Step 1: build the dependency prefix.  Only needed once - re-use `deps_prefix` for later builds,
#  and delete `deps_build` afterwards (it is several GB of scratch).  Note the `bash -lc`: a login
#  shell is what puts gcc-toolset-14 on PATH in this image.
docker run --rm -it -v `pwd`/InterSpec_code:/interspec -v `pwd`/deps_build:/deps_build -v `pwd`/deps_prefix:/deps_prefix quay.io/pypa/manylinux_2_28_x86_64:latest bash -lc '/interspec/target/electron/build_linux_deps_from_docker.sh /interspec /deps_build /deps_prefix'

# Step 2: compile.  Leaves the app payload in build_electron/app, and fails the build if the
#  binaries would need a newer glibc/libstdc++ than we support (see check_elf_compat.py).
#  Drop the trailing `/deps_prefix` argument to use FetchContent instead of the prefix.
docker run --rm -it -v `pwd`/InterSpec_code:/interspec -v `pwd`/build_electron/:/build_app -v `pwd`/build_working_dir:/build_working_dir -v `pwd`/deps_prefix:/deps_prefix quay.io/pypa/manylinux_2_28_x86_64:latest /interspec/target/electron/build_linux_app_from_docker.sh /interspec /build_app /build_working_dir /deps_prefix

# Step 3: package.  Needs dpkg-deb, rpmbuild and fakeroot available.
cd InterSpec_code/target/electron
npm install --ignore-scripts
npx electron-builder --linux -c.directories.app=/tmp/build_interspec/build_electron/app
# Results are in target/electron/dist/
```

Note the `--ignore-scripts`: `package.json` has an `install` script of `cmake-js compile`, which
would otherwise try to build the C++ again on the host.

`build_linux_deps_from_docker.sh` is a thin wrapper: it checks the container has the toolchain it
expects (it deliberately installs nothing, so a missing tool is an error rather than a silently
different compiler), runs `target/patches/dep_build_linux.sh`, then verifies the resulting prefix
and writes a `.interspec_deps_complete` sentinel recording what built it. Step 2 refuses a prefix
without that sentinel. If you want a prefix outside a container — for ordinary development rather
than for building a redistributable package — run `target/patches/dep_build_linux.sh` directly; see
`target/patches/README.md`.

### Building manually using a manylinux container
```bash
# From your host OS terminal, run the following commands
git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_linux_electron_build
cd InterSpec_linux_electron_build/

# Grab the the oldest (currently) supported manylinux image to your machine
# (Electron 28+ requires glibc 2.28, so we use manylinux_2_28)
docker pull quay.io/pypa/manylinux_2_28_x86_64:latest

# Start a shell session within the image, mapping the InterSpec source 
#  directory to /interspec.  We'll also map port 8081 for testing.
docker run --rm -it -v `pwd`:/interspec quay.io/pypa/manylinux_2_28_x86_64:latest bash

# Get the dependancies we need to build InterSpec
yum update
yum install -y npm

# Make and cd into build directory - note this is in host filesystem incase we  
#  want to come back to things, but dont want to rebuild everything from scratch
cd /interspec/target/electron/

# Install the NPM dependancies we need
npm install uglify-js -g
npm install uglifycss -g
npm install cmake-js -g
npm install --save-dev node-addon-api --arch=x64
npm install node-api-headers


# This next command will take like 10 or 20 minutes to clone into the boost and Wt repositories
#  Note that libstdc++/libgcc are linked dynamically - see "Linux Considerations" below and
#  cmake/CxxRuntimePolicy.cmake; do not add `-static-libstdc++` here.
#  This command may take a long time to run.
CMAKE_BUILD_PARALLEL_LEVEL=6 cmake-js --architecture x64 --arch=x64 --CDCMAKE_BUILD_TYPE="Release" --CDInterSpec_FETCH_DEPENDENCIES=ON --CDBUILD_AS_LOCAL_SERVER=OFF --CDUSE_LEAFLET_MAP=ON --CDLEAFLET_MAPS_KEY="..." --CDUSE_REL_ACT_TOOL=ON --out=build_manylinux_electron --target install

# Then, on a normal Linux host (not in the manylinux container - it has no dpkg/rpmbuild):
npm install --ignore-scripts
npm run package-manylinux

# dist/ now holds the .deb, .rpm, .AppImage and linux-unpacked/.  To run without installing:
./dist/linux-unpacked/interspec
```


## Minimum supported OS
Electron 41 sets these floors for end-user systems:
- macOS 11+ (Big Sur), Apple Silicon (arm64) only — Intel macOS is not built.
- Windows 10 1809 (build 17763) or newer, x64 only — ia32 / 32-bit Windows is no longer built (Electron dropped it after Electron 21).
- Linux x64 with glibc 2.28+ — RHEL 8, Ubuntu 20.04+, Fedora 30+, Debian 11+. Older distros (CentOS 7, Ubuntu 18.04, etc.) are not supported.


## Code signing & notarization (macOS)

> Note: the **released** macOS app is built from `target/osx` (see
> `.github/workflows/build_macos.yml`), not from this Electron target. `npm run package-mac` is a
> development convenience, and its signing path is not exercised by CI.

`npm run package-mac` runs electron-builder, which signs and notarizes via `@electron/osx-sign`
and `@electron/notarize`. electron-builder reads its credentials from environment variables rather
than from CLI flags.

**Signing.** With nothing set, electron-builder searches the keychain for a valid Developer ID and
uses it; if it finds none it *skips signing entirely* - there is no ad-hoc fallback. To pin a
specific certificate, set `CSC_NAME` (e.g. `Developer ID Application: Your Org (TEAMID)`), or
`CSC_LINK` + `CSC_KEY_PASSWORD` for a `.p12`. `CSC_IDENTITY_AUTO_DISCOVERY=false` forces an
unsigned build.

**Notarization** only runs when one of these complete sets is present - a partial set means it is
silently skipped:

| Option | Variables |
| --- | --- |
| App Store Connect API key (recommended) | `APPLE_API_KEY`, `APPLE_API_KEY_ID`, `APPLE_API_ISSUER` |
| Apple ID | `APPLE_ID`, `APPLE_APP_SPECIFIC_PASSWORD`, `APPLE_TEAM_ID` |
| `notarytool` keychain profile | `APPLE_KEYCHAIN` **and** `APPLE_KEYCHAIN_PROFILE` |

Note the keychain option needs *both* variables: `APPLE_KEYCHAIN` is the keychain **path**
(`~/Library/Keychains/login.keychain-db` for the login keychain, which is where
`notarytool store-credentials` writes by default), and `APPLE_KEYCHAIN_PROFILE` is the profile
name. `APPLE_TEAM_ID` is not required for this option. Set `-c.mac.notarize=false` to skip
notarization even when credentials are present.

Ad-hoc signing (`-c.mac.identity=-`) works here because
`macOS/entitlements.mac.plist` already grants `com.apple.security.cs.disable-library-validation`;
without that entitlement, hardened runtime rejects Electron's pre-signed frameworks (different
Team ID) and the app will not launch.

### One-time setup: store the notarization credential in the macOS keychain
Apple's `notarytool` reads notarization credentials from a named profile in the login keychain, so the Apple-ID password never has to live in environment variables, shell history, or build scripts.

```bash
xcrun notarytool store-credentials "interspec-notary" \
  --apple-id <your-apple-id> \
  --team-id <YOUR_TEAM_ID> \
  --password <app-specific-password>
```

Generate the app-specific password at <https://appleid.apple.com/> → "Sign-In and Security" → "App-Specific Passwords"; do **not** use your iCloud password. The first time you run `store-credentials` macOS will prompt for your login-keychain password (i.e. your Mac account password) so it can encrypt the credential at rest. After that, no further prompts.

To use a different profile name (e.g. when sharing a build machine), point `APPLE_KEYCHAIN_PROFILE`
at it. `codesign-electron.sh` reads `NOTARY_KEYCHAIN_PROFILE` for the same purpose.

### Re-signing an existing bundle
For one-off re-signing of an existing bundle (out-of-band from the packaging step), `target/electron/macOS/codesign-electron.sh path/to/InterSpec.app` is a thin wrapper around the same tooling and uses the same keychain profile.


## Linux Considerations

### -fPIC
The `InterSpec` module is really a shared library that node.js loads, therefore you need the `-fPIC` C/C++ compiler flag enabled not just for the `InterSpec` code, but for all of the static libraries you link it against, including boost, Wt, and zlib - which isnt the default when compiling static libraries for any of them, so when building them you may want to add `-fPIC` to the flags.

### Do not statically link libstdc++
This one has bitten us: [issue #51](https://github.com/sandialabs/InterSpec/issues/51) was every
Electron 41 Linux build segfaulting before `main()`, and the cause was
`-static-libgcc -static-libstdc++` on `InterSpecAddOn.node`.

With libstdc++ linked statically, its locale-facet initialization becomes just another static
initializer inside our shared object, and the toolchain ordered it *after* Wt's. So the
namespace-scope `boost::basic_regex` in Wt's `CgiParser.C` ran first, called through
`std::collate<char>`, and landed in the wrong vtable slot
(`std::codecvt<char16_t,...>::do_unshift`). Linking libstdc++ dynamically makes it a real shared
object, and the loader is then required to run its `dl_init` before ours.

Linking dynamically does **not** raise the supported-OS floor, which is why the static link was
there in the first place. The release build uses `quay.io/pypa/manylinux_2_28_x86_64`, whose
compiler is Red Hat `gcc-toolset-14`. Developer Toolset links the *base* AlmaLinux 8
`libstdc++.so.6` and pulls only the newer ABI bits from `libstdc++_nonshared.a`, so the binaries
require at most `GLIBCXX_3.4.25` / `CXXABI_1.3.11` / `GLIBC_2.28` even though they are compiled by
GCC 14.

The policy lives in `cmake/CxxRuntimePolicy.cmake` (one place, applied to every consumer of
`InterSpecLib` regardless of link type) and is enforced by
`target/electron/linux/check_elf_compat.py`, which the container build script and CI both run. It
fails the build if a binary needs a newer runtime than we advertise, or if libstdc++ stops being a
dynamic dependency. To check a build by hand:

```bash
python3 target/electron/linux/check_elf_compat.py /path/to/build/app
```

### The Chromium sandbox
Chromium's sandbox needs `chrome-sandbox` to be owned by root with mode 4755, or working
unprivileged user namespaces. How each artifact copes:

| Artifact | Sandbox |
| --- | --- |
| `.deb` / `.rpm` | electron-builder's post-install script sets up `chrome-sandbox` and installs an AppArmor profile granting `userns` into `/etc/apparmor.d/` (needed on Ubuntu 23.10+, which sets `kernel.apparmor_restrict_unprivileged_userns=1`). The profile keys on the executable path, which is why the installed binary must not be wrapped in a launcher script. |
| `.AppImage` | Runs from a `nosuid` FUSE mount, so a setuid helper can never work; electron-builder's `AppRun` passes `--no-sandbox`. |
| `.tar.gz` | `linux/interspec-launcher.sh` ships as `interspec` and probes for a usable sandbox (setuid-root `chrome-sandbox` **and** a working `unshare -Ur`), exec'ing `interspec-bin` with or without `--no-sandbox` accordingly. An explicit `--no-sandbox`/`--enable-sandbox` from the user is passed through untouched. |

Note that `app.commandLine.appendSwitch('no-sandbox')` from `main.js` does **not** work on Linux:
Chromium initializes the sandbox and forks the zygote before the JS entry point is evaluated, so
the switch only counts if it is on the real argv. `main.js` therefore does not try. Instead it
detects a renderer death that is attributable to the sandbox and relaunches once via
`app.relaunch()` with `--no-sandbox` on the real argv.

"Attributable to the sandbox" is deliberately narrow, because the relaunch calls `app.exit()`,
which takes the whole application down without running `before-quit` or the window close handlers
— so a false positive discards unsaved work. It requires all of: the **first** window of the
launch (a later window's renderer dying means an established session is open), a death **within
60 s** of that window being created (a broken sandbox fails immediately), and a `reason` of
`launch-failed`, `crashed`, or `abnormal-exit` — explicitly **not** `oom` or `killed`, since
InterSpec can legitimately exhaust memory on a large file and that is not a sandbox problem.

Failures are counted in `sandbox_failures.json` in the user-data directory. It takes **two
consecutive** failed sandboxed launches before InterSpec stops attempting one, so a single
one-off crash cannot silently disable the sandbox; and any successful sandboxed session deletes
the file, so the state self-heals once the underlying cause is gone. Delete it by hand to retry
immediately.

As of Electron 41 the sandbox is broken on Ubuntu 24.04+ even *with* the AppArmor profile - the
renderer dies on `/dev/shm` access - which appears to be an upstream Chromium regression
(Electron 13 was fine). The relaunch fallback above is what keeps those systems usable.

### Desktop integration
- `linux/gov.sandia.InterSpec.xml` defines the `application/gamma-spectrum` MIME type and its
  globs. It is **not** redundant with the file electron-builder generates from
  `fileAssociations`: that one joins the whole extension list into a single malformed glob.
- `linux/gov.sandia.InterSpec.metainfo.xml` is the AppStream metadata, for GNOME Software /
  Discover. Validate with `appstreamcli validate`.
- `linux/icons/*.png` are extracted from `macOS/InterSpec.icns` (which carries
  512x512 and 1024x1024 PNG payloads) and installed into the hicolor theme.
- `desktopName` in the generated app `package.json` and `linux.syncDesktopName` keep the
  `.desktop` filename, `StartupWMClass`, and Electron's `app_id` equal, so window-to-launcher
  association works.


## Future Work
- In the future the build process may be improved to do the final packaging through CMake.
- Flatpak was considered and deferred: it needs `flatpak-builder` plus runtimes in CI, and
  Flathub distribution is a separate submission process. The AppStream metainfo above is a
  prerequisite for it, so most of the metadata work is already done.
- The AppStream metadata has no `<releases>` block, so GNOME Software cannot show "what's new".
  Adding one means keeping version entries up to date.
- Linux is x86_64 only. An arm64 build would need an arm64 manylinux container (QEMU, or an
  arm64 runner).
