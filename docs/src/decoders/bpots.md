# BP-OTS Decoder

The [`BPOTSDecoder`](@ref) runs belief propagation with Ordered Trapping Set (OTS) biasing. It detects oscillating variable nodes (trapping sets) that stall convergence and applies a corrective bias periodically to get unstuck.

Biasing schedule:
1. Run standard BP message passing.
2. Every `T` iterations, find the variable node with the highest oscillation count.
3. Apply a negative bias (`-C`) to that node's prior belief to push it toward flipping.
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
- `max_iters::Int` — Maximum BP iterations.
- `T::Int` — Bias period (default `9`).
- `C::Float64` — Bias strength (default `2.0`).

## API

See the [API Reference](../api.md) for full docstrings.
