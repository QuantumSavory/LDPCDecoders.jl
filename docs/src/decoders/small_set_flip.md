# Small-Set Flip (SSF) Decoder

The [`SmallSetFlipDecoder`](@ref) implements the Small-Set Flip algorithm for quantum LDPC codes (Leverrier et al., 2015). It's a greedy iterative decoder: instead of flipping one bit at a time like classical bit-flip decoding, SSF looks at small subsets of qubits in each check's support and flips the subset that gives the biggest drop in total syndrome weight.

For CSS quantum codes, [`CSSSmallSetFlipDecoder`](@ref) wraps two SSF decoders to handle X-errors and Z-errors independently.

## How it Works

1. Compute syndrome $s = H \cdot e \pmod 2$.
2. For each check (row in $H$), test all qubit subsets of size $\le t$ within its support.
3. Find the candidate subset that flips the maximum number of unsatisfied checks back to satisfied.
4. Apply the flip and repeat until the syndrome is zero or no flip improves the syndrome weight.

Setting `t = 1` makes this equivalent to standard single-bit-flip decoding. Setting `t = 0` (or leaving it default) searches all subsets up to the check weight.

## Usage

### Standard Decoder

```julia
using LDPCDecoders

H = LDPCDecoders.parity_check_matrix(100, 5, 4)
decoder = SmallSetFlipDecoder(H, 50; t=2)

error = rand(100) .< 0.02
syndrome = Bool.((H * error) .% 2)

guess, converged = decode!(decoder, syndrome)
```

### CSS Decoder

Decodes X and Z errors separately using $H_z$ for X-errors and $H_x$ for Z-errors:

```julia
using LDPCDecoders

decoder = CSSSmallSetFlipDecoder(Hx, Hz, 50; t=2)

# Pass separate X and Z syndromes
(guess_x, guess_z), converged = decode!(decoder, syn_x, syn_z)

# Or pass a combined syndrome vector [syn_x; syn_z]
guess_combined, converged = decode!(decoder, vcat(syn_x, syn_z))
```

### Batch Decoding

```julia
samples = 100
errors = rand(100, samples) .< 0.02
syndromes = Bool.((H * errors) .% 2)

guesses, successes = batchdecode!(decoder, syndromes, zeros(Int, size(H, 2), samples))
```

## Parameters

- `H` — Parity check matrix (`SparseMatrixCSC{Bool,Int}` or `BitMatrix`).
- `max_iters::Int` — Maximum greedy flip iterations. A value around `50` is usually enough.
- `t::Int` — Maximum subset size evaluated per check (default `0` = full check weight).
- `rng` — Optional random number generator for randomized tie-breaking when multiple flips yield equal syndrome weight reduction.

| `t` value | Subsets tested per check (degree-5) | Trade-off |
|---|---|---|
| `t=1` | 5 | Fastest, single-bit flip |
| `t=2` | 15 | Good balance between speed and correction |
| `t=3` | 25 | Better error correction, slightly slower |
| `t=0` (full) | 31 | Maximum correction power |

## References

- Leverrier, Tillich & Zémor, *Quantum Expander Codes* ([arXiv:1504.00822](https://arxiv.org/abs/1504.00822), 2015)
- Fawzi, Grospellier & Leverrier, *Efficient decoding of random errors for quantum expander codes* ([arXiv:1711.08351](https://arxiv.org/abs/1711.08351), 2018)
- Grospellier & Krishna, *Numerical study of hypergraph product codes* ([arXiv:1810.03681](https://arxiv.org/abs/1810.03681), 2019)

## API

See the [API Reference](../api.md) for full docstrings.
