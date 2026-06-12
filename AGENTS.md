# Repository Guidelines

## Project Structure & Module Organization
Core simulation logic lives in `src/` (`pfqmc.cpp`, `skewMatUtils.cpp`), with shared headers and the PFAPACK Fortran/C bridge under `inc/` (build `inc/pfapack/{fortran,c_interface}` before linking). `main.cpp` drives CLI experiments, while executable and object artifacts land in `bin/` and `obj/`. GoogleTest fixtures are in `test/`, and reference derivations plus Hubbard study notes sit in `_research/` and the `HUBBARD_*.md` sheets—treat them as canonical physics context when adding operators or observables.

## Build, Test, and Development Commands
- `cd inc/pfapack/fortran && make mFC=ifx` then `cd ../c_interface && make` to refresh Pfaffian libraries whenever you change compilers.
- `./local-build.sh [--debug|--release|--clean-only]` sources Intel oneAPI, wipes `build/`, and runs CMake with the Eigen include path; defaults to Debug, so pass `--release` for benchmarks or `--` to forward extra CMake flags.
- `cmake -S . -B build -DEIGEN3_INCLUDE_DIR=/path/to/eigen3` is the manual alternative if you want a separate build tree; follow with `cmake --build build --target main` (or `main_test`).
- `ctest --test-dir build --output-on-failure` or `./build/main_test` executes the full GoogleTest suite after a configure/build.
- `make clean` only applies when using the legacy Makefile rather than CMake; the script already removes stale artifacts each run.

## Coding Style & Naming Conventions
Target C++17 with 4-space indentation, LLVM brace placement, and descriptive camel-case for types (`PfQMC`, `SpinlessTvUtils`) plus lower camelCase for methods (`leftSweep`). Keep Eigen aliases (`MatType`, `DataType`) in headers, prefer `const` references for matrices, and document stabilization logic with brief comments rather than block prose. Run `clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4}" -i <files>` before committing and ensure headers under `inc/` stay self-contained. Clangd users should keep the repo-root `.clangd` to pick up Intel oneAPI include paths so diagnostics mirror the `mpiicpx`/MKL toolchain.

## Testing Guidelines
Unit tests rely on GoogleTest (see `test/main_test.cpp`, `test/squareLatticeTest.cpp`). Mirror production namespaces in test fixture names (`SquareLatticeTest.CheckPropagation`) and cover both numerical accuracy and sign-handling paths. Any change to propagator math must add or extend a deterministic test; stochastic pieces should expose a seeded path. Execute `ctest -R <pattern>` for focused runs and avoid skipping tests in CI—coverage is thin, so every new feature should introduce at least one assertion that protects it.

## Commit & Pull Request Guidelines
History shows concise, imperative commits ("fix sign in Number", "remove two point term when mu=0"); follow that style, referencing issues as `Fix: ... (#42)` when relevant. Each PR should summarize the physics change, list required environment variables (e.g., lattice size, `beta`), and include `ctest` output plus any generated `.out` files that prove convergence. When touching numerics, attach small tables or plots comparing prior runs to justify deltas before requesting review.

## Security & Configuration Tips
Keep MKL and PFAPACK builds in sync with the compiler used by CMake/Make to avoid ABI mismatches. Do not commit Intel oneAPI environment dumps or raw lattice datasets larger than ~5 MB; store them under a reproducible script in `_research/` instead. Validate new environment variables in `README.md` and scrub secrets from shell snippets.
