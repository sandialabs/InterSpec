# target/peak_fit_improve

R&D sandbox for peak-fitting algorithm research. Contains genetic-algorithm experiments
for candidate-peak detection, initial/final peak-shape fits, detector-type classification,
and nuclide configuration tuning.

This directory has its own `CMakeLists.txt` but is **not** part of the main InterSpec
build — it is not added via `add_subdirectory()` from the top-level project. Build it
standalone if you want to iterate on the underlying fitting heuristics; the artifacts
produced here have historically been folded back into the main `PeakFit.cpp` codepath
by hand after tuning runs.

Despite living under `target/`, this is not a deployable build target; it is a research
workspace kept here for convenience.
