# Belief Propagation + OSD Decoder

The [`BeliefPropagationOSDDecoder`](@ref) combines standard belief propagation with Ordered Statistics Decoding (OSD) as a post-processing step. After BP produces soft decisions (log-likelihood ratios), columns of the parity check matrix are sorted by reliability and Gaussian elimination is performed on the least reliable subset to find a minimum-weight correction.

This decoder is particularly effective when standard BP fails to converge, as OSD can still recover a valid correction from the soft information.

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

- `H::BitMatrix` — The parity check matrix (must be a `BitMatrix`).
- `per::Float64` — Physical error rate.
- `max_iters::Int` — Maximum BP iterations.
- `osd_order::Int` — Order of OSD post-processing (default `0`). Higher orders search more candidate corrections but scale exponentially.

## API

