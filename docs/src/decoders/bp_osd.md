# Belief Propagation + OSD Decoder

The [`BeliefPropagationOSDDecoder`](@ref) runs BP to get soft decisions (log-likelihood ratios), then sorts parity check matrix columns by reliability and applies Gaussian elimination on the unreliable columns to find a minimum-weight correction.

OSD post-processing lets the decoder recover a valid correction even when BP doesn't converge.

## Usage

```julia
using LDPCDecoders

H = BitMatrix(LDPCDecoders.parity_check_matrix(1000, 10, 9))
decoder = BeliefPropagationOSDDecoder(H, 0.01, 100; osd_order=0)

error = rand(1000) .< 0.01
syndrome = Bool.((H * error) .% 2)

guess, converged = decode!(decoder, syndrome)
```

## Parameters

- `H::BitMatrix` — Parity check matrix (must be a `BitMatrix`).
- `per::Float64` — Physical error rate.
- `max_iters::Int` — Maximum BP iterations.
- `osd_order::Int` — OSD post-processing order (default `0`). Higher orders check more
  candidate corrections but scale exponentially.

## API

See the [API Reference](../api.md) for full docstrings.