# BP-OTS Decoder

The [`BPOTSDecoder`](@ref) implements Belief Propagation with Ordered Trapping Set (OTS) biasing. This algorithm extends standard BP by detecting oscillating variable nodes (trapping sets) that prevent convergence, and periodically applying a bias to break them out of decoding failures.

The algorithm follows a biasing schedule controlled by the period `T` and bias constant `C`:
1. Run standard BP message passing.
2. Every `T` iterations, identify the variable node with the highest oscillation count.
3. Apply a negative bias (`-C`) to that node's prior belief to push it toward being flipped.
4. Continue BP with the modified beliefs.

## Usage

```julia
using LDPCDecoders

H = BitMatrix(LDPCDecoders.parity_check_matrix(1000, 10, 9))
decoder = BPOTSDecoder(H, 0.01, 100; T=9, C=2.0)

error = rand(1000) .< 0.01
syndrome = Bool.((H * error) .% 2)

guess, converged = decode!(decoder, syndrome)
```

## Parameters

- `H` — Parity check matrix (`BitMatrix` or `SparseMatrixCSC{Bool,Int}`).
- `per::Float64` — Physical error rate.
- `max_iters::Int` — Maximum number of BP iterations.
- `T::Int` — Biasing period (default `9`). How often the OTS bias is applied.
- `C::Float64` — Bias constant (default `2.0`). Strength of the bias applied to oscillating nodes.

## API

