It shouldnt be to hard to integrate InterSpec into other applications by loading the InterSpec shared library, and calling into some "extern C" defined functions (see ElectronUtils.h). 

Currently only a python example is provided, and assumes things were CMake option BUILD_AS_ELECTRON_APP as ON, and INTERSPEC_LIBRARY_STATIC is set to off.  Build options are likely to change in the future to better accomidate this mode.


