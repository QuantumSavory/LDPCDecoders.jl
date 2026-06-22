# Belief Propagation Decoder

The [`BeliefPropagationDecoder`](@ref) runs the sum-product algorithm on the Tanner graph of the parity check matrix, passing LLR messages between variable and check nodes until the syndrome is zero or the iteration limit is reached.

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

For many decoding trials, [`batchdecode!`](@ref) reuses the same scratch space across
calls:

```julia
samples = 100
errors = rand(1000, samples) .< 0.01
syndromes = H * errors .% 2

guesses, successes = batchdecode!(decoder, syndromes, zero(errors))
```

## API

See the [API Reference](../api.md) for full docstrings.