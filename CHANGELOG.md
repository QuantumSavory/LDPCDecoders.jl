# News

## v0.3.4 - dev

- Migrate test infrastructure from SafeTestsets to TestItemRunner with tag-based filtering
- Make JET a conditional test dependency (only installed when `JET_TEST=true`)
- bump julia compat to 1.12
- fix deprecated `parity_checks_x` API call in BP-OTS tests (now `parity_matrix_x`)
- Add generic `batchdecode!` fallback and widen `BPOTSDecoder` `decode!` to accept `AbstractVector`
- Use `decode!` instead of deprecated `syndrome_decode!` in BP `batchdecode!`
- Set up Documenter.jl docs page with API reference and decoder guides
- Add docstrings to all public types and methods
- Fix pre-existing docstring bugs (broken code blocks, undefined variables)
- Add benchmark suite (`benchmark/benchmarks.jl`) for AirspeedVelocity CI
- Add BP-OTS `decode!` benchmark for tracking message-passing performance
- Fix BP-OTS decoder and benchmark script to correctly handle integer syndromes
- Implement OSD-0 fast shortcut (Algorithm 2 from Panteleev & Kalachev) for BP-OSD decoder
- Fix `rowswap!` to use zero-allocation element-wise swap instead of slice-based copy
- Clean up outdated developer comments and TODOs across the codebase
- Remove 10 deprecated/legacy decoder files (superseded by modern API) and old examples
- Finalize legacy cleanup by completely removing unused syndrome decoders and generator
- Add a read-only organization audit for Documenter deployment setup and enable pull-request documentation previews

## v0.3.3 - 2025-04-15

- Added an implementation for BP-OTS decoder.

## v0.3.2 - 2024-11-15

- Add a (still unoptimized) implementation of a BP OSD decoder.

## Older - before 2021-10-28 unrecorded
