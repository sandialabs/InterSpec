# src/js_inline

These JavaScript files are compiled into the InterSpec executable. They are loaded by
the C++ via Wt's `LOAD_JAVASCRIPT(...)` macro, and the macro reads them at build time
through a direct `#include` of the .js file (see `src/InterSpec.cpp`, `src/SpectrumChart.cpp`,
`src/SpecFileQueryWidget.cpp`, `src/TerminalWidget.cpp`).

They are **not** runtime web assets. For JavaScript that is served to the browser as a
separate resource (D3, D3-based chart code, theme scripts, etc.), see `InterSpec_resources/`.

Because each file here is inlined into the C++ translation unit that includes it, editing
any of these files requires recompiling that unit.
