# BPFC 0.3.0 check notes

Local source-package checks used R 4.6.0 on macOS and NIMBLE 1.4.2.
The package includes tests comparing the conditional innovation density with
an independent dense Gaussian calculation, and tests for cluster selection
and posterior summaries.

Fitting compiles a NIMBLE model and requires a C++ toolchain. The distribution
registration exposes its NIMBLE density and simulation functions under the
names required by the engine.
