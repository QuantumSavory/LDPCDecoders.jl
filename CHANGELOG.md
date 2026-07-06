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
- Clean up outdated developer comments and TODOs across the codebase

## v0.3.3 - 2025-04-15

- Added an implementation for BP-OTS decoder.

## v0.3.2 - 2024-11-15

- Add a (still unoptimized) implementation of a BP OSD decoder.

## Older - before 2021-10-28 unrecorded
