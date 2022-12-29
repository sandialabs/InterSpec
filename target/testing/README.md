This directory will contain the publicly avaiable tests, meant to be run as part of CI.  
These tests are planned to be added to over 2023, and only represent a portion of InterSpecs testing.  

InterSpecs testing consists of:
- The primary developer, is also probably the primary user, so any issue identified over the course of using the application, that could quietly cause an incorrect answer, is agressively addressed.
- There are a couple thousand `assert` statement throughout the code that help to ensure application logic is as expected.  These checks only get ran on 'debug' builds, which is what the primary developer/user normally uses.
- There is a compile-time `PERFORM_DEVELOPER_CHECKS` option that enables a few hundred sections of code that will do something like compute a quantity in an independent method, or perform a more expensive consistency-check, and log any issues found.  A number of these sections of code contain functions `static void equalEnough( const Type &lhs, const Type &rhs )` functions that test if class objects are (de)serialize correctly.
- The `SpecUtils` library that handles all the spectrum file related activity, and many of the string manipulation, filesystem, and time-related functions, and has its own [unit tests](https://github.com/sandialabs/SpecUtils/tree/master/unit_tests), [fuzz tests](https://github.com/sandialabs/SpecUtils/tree/master/fuzz_test), and a [regression test](https://github.com/sandialabs/SpecUtils/tree/master/regression_test).  The regression test ensures spectrum files parse identically to previously manually verified parsings, as well that every read file is able to be written to N42-2012 files, and read back in without losing any information; the test corpus of a couple hundred unique spectrum file format variants is not public.
- There are lots of consistency checks built into the code, that run on all builds, and issue errors to user if encountered.
- The [SpectrumViewerTester](../../InterSpec/SpectrumViewerTester.h) class implements a kind of end-to-end testing that read in N42 files saved from InterSpec with peak fits and activity/shielding fit results, and ensures those results can be reproduced (i.e., makes sure a bug isnt introduced that causes a wrong activity or shielding to then be fit for) through simulating a user session and doing things like double-clicking on the spectrum to fit peaks and comparing them to what the analyst used.
    - Note though that the session isnt rendered, its "driven" in a client-less mode from c++ - tests of the rendered session, using something like Puppeteer, has been experimented with, but isnt currently being used, but may be used in the future.
    - The few hundred spectra with analyst-verified solutions that are checked against, are not currently publicly avaiable, but some publicly avaiable ones are planned to be added over CY 2023.
- There are some more unit-like tests that are not public, as they use proprietary data in thier test cases.
- The tests in this directory.


To build:
```bash
cd InterSpec/target/testing
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/prefix ..

# Build the code 
make -j12
# Or
cmake --build . --config Debug -j12

# Run the tests
ctest
# Or you can run the tests more generically with
cmake --build . --target RUN_TESTS
```