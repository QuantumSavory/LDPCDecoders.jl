# Belief Propagation Decoder

The [`BeliefPropagationDecoder`](@ref) implements the standard sum-product belief propagation algorithm for decoding LDPC codes. It operates on the Tanner graph representation of the parity check matrix, passing log-likelihood ratio (LLR) messages between variable and check nodes until convergence or a maximum iteration count is reached.

## Usage

```julia
using LDPCDecoders

H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
decoder = BeliefPropagationDecoder(H, 0.01, 100)

error = rand(1000) .< 0.01
syndrome = (H * error) .% 2

guess, converged = decode!(decoder, syndrome)
```

## Batch Decoding

For running many decoding trials efficiently, use [`batchdecode!`](@ref) which re-uses the same pre-allocated scratch space:

```julia
samples = 100
errors = rand(1000, samples) .< 0.01
syndromes = H * errors .% 2

guesses, successes = batchdecode!(decoder, syndromes, zero(errors))
```

## API

