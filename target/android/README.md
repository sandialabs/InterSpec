# Building InterSpec for Android

InterSpec uses CMake `FetchContent` to download and build Boost, Wt, Eigen, and Ceres
automatically. The first build takes 30–60 minutes; subsequent builds after a C++ change
take only a few minutes because compiled objects are cached.

Two build methods are available: **Docker** (reproducible, no Android Studio required)
and **local** (Android Studio or command-line Gradle).

---

## Requirements

| Component        | Version                 |
|------------------|-------------------------|
| NDK              | r28 (`28.0.13004108`)   |
| AGP              | 8.5.1                   |
| Gradle           | 8.9                     |
| compileSdk       | 35                      |
| minSdk           | 21 (Android 5.0)        |
| JDK              | 17                      |

### Product Flavors

| Flavor      | ABIs                      | Use case                          |
|-------------|---------------------------|-----------------------------------|
| `arm8`      | `arm64-v8a`               | Physical devices (fastest build)  |
| `universal` | `arm64-v8a`, `x86_64`     | Devices + emulator testing        |

---

## Docker Build (recommended for CI or first-time setup)

Docker provides a fully reproducible build environment. No Android SDK or NDK
installation is needed on the host — only Docker.

### 1. Build the SDK image (once)

```bash
cd target/android
./build_android_docker.sh image
```

This creates the `interspec-android-sdk` Docker image (~3–4 GB) containing JDK 17,
Android SDK, NDK r28, CMake, and Ninja.

### 2. Build the APK and AAB

```bash
./build_android_docker.sh build              # default: arm8Release
./build_android_docker.sh build arm8Debug
./build_android_docker.sh build universalRelease
```

On first run this will:
- Fetch Wt 3.7.1 source (for the `resources/` directory needed by InterSpec assets)
- Create the `interspec-assets.zip` that gets bundled into the APK
- Run Gradle, which triggers CMake `FetchContent` to download and compile all C++ dependencies

Build artifacts (`.apk` and `.aab`) are copied to `build_android_output/`.

Persistent caches are stored in `.docker-cache/` so subsequent builds are fast:
- `.docker-cache/gradle/` — Gradle daemon cache and downloaded dependencies
- `.docker-cache/app-cxx/` — CMake build tree (compiled `.o` files)
- `.docker-cache/app-build/` — Gradle/APK packaging output

### 3. Interactive shell (debugging)

```bash
./build_android_docker.sh shell
```

Drops you into the container with source and caches mounted, so you can run
`./gradlew` commands or inspect the build tree manually.

### 4. Install on a device

From the host (with `adb` available):

```bash
adb install build_android_output/app-arm8-release.apk
```

---

## Local Build (Android Studio or command line)

### Prerequisites

- Android Studio (or standalone Android SDK with NDK r28 and JDK 17)
- Wt resources directory available from a prior local build (e.g., `build_xcode/resources/`)

### Android Studio

Open `target/android/InterSpec` as a project. Android Studio will sync Gradle and
detect the NDK version. Select the desired build variant (e.g., `arm8Release`) and
build/run.

### Command Line

```bash
cd target/android/InterSpec

# Set JAVA_HOME if not already configured
export JAVA_HOME="/Applications/Android Studio.app/Contents/jbr/Contents/Home/"

# Build
./gradlew assembleArm8Release           # ARM64 release APK
./gradlew bundleArm8Release             # ARM64 release AAB (for Play Store)
./gradlew assembleUniversalDebug        # ARM64 + x86_64 debug APK

# Install to connected device
./gradlew installArm8Release
```

The `make_assets.sh` script runs automatically before the build. It looks for Wt
resources in these locations (in order):
1. `WT_RESOURCES_DIR` environment variable
2. `build_xcode/resources/` (from repo root)
3. `build_vscode/resources/`
4. `build/resources/`

---

## Signing

Signing is **optional** — unsigned/debug builds work by default.

To produce a signed release, provide keystore credentials via environment variables
or `local.properties`. Environment variables take precedence.

### Environment variables (recommended for CI/Docker)

```bash
export KEYSTORE_FILE=~/.android/release.keystore   # host path (Docker bind-mounts it)
export KEYSTORE_PASSWORD=changeit
export KEYSTORE_ALIAS_NAME=interspec
export KEYSTORE_ALIAS_PASSWORD=changeit
./build_android_docker.sh build
```

### local.properties (for local/Android Studio builds)

Add to `target/android/InterSpec/local.properties`:

```properties
keystore.location=/path/to/release.keystore
keystore.password=changeit
keystore.alias.name=interspec
keystore.alias.password=changeit
```

### Creating a keystore

If you don't have one yet:

```bash
keytool -genkeypair -v -keystore release.keystore -alias interspec \
    -keyalg RSA -keysize 2048 -validity 10000
```

### Verifying a signed APK

```bash
apksigner verify --verbose build_android_output/app-arm8-release.apk
```

---

## Leaflet Maps Key

The Leaflet maps API key defaults to the value in `app/build.gradle`. To override it
(e.g., for a different deployment), set the `LEAFLET_MAPS_KEY` environment variable
before building.

---

## Troubleshooting

**First build is very slow (30–60 min)**
This is expected — `FetchContent` downloads and compiles Boost, Wt, Eigen, and Ceres
from source. Subsequent builds reuse the cached objects.

**Android Studio is unresponsive during native builds**
The NDK compilation is CPU and memory intensive. Building from the command line
(`./gradlew`) is often faster and more predictable. The `gradle.properties` file
allocates 4 GB heap (`-Xmx4096m`).

**`make_assets.sh` fails with "Cannot find Wt resources directory"**
You need Wt resources available. Either:
- Build InterSpec locally first (e.g., in `build_xcode/`), or
- Set `WT_RESOURCES_DIR` to point to a Wt resources directory, or
- Use the Docker build (which fetches Wt source automatically)

**Clean rebuild**
```bash
# Docker: remove cached build trees
rm -rf target/android/.docker-cache/app-cxx target/android/.docker-cache/app-build

# Local:
cd target/android/InterSpec
./gradlew clean
```
