# target/sharedlib

**Status: reference / not actively maintained.**

Demonstrates integrating InterSpec into other applications by loading it as a shared
library and calling `extern "C"` functions (see `ElectronUtils.h`). A Python example is
included; it assumes the CMake option `BUILD_AS_ELECTRON_APP` is `ON` and
`INTERSPEC_LIBRARY_STATIC` is `OFF`.

Build options here are likely to change as the shared-library integration mode is
revisited; treat the contents of this directory as a starting point, not a stable API.
