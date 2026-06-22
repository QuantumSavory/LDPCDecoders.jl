# Iterative BitFlip Decoder

The [`BitFlipDecoder`](@ref) runs the classical hard-decision bit-flip algorithm. Each iteration computes a syndrome, counts unsatisfied check votes per bit, and flips the top candidate — repeating until the syndrome is zero or the vote counts stop changing.

It's the simplest iterative decoder and a reasonable baseline for comparing against heavier algorithms.

## Usage

```julia
using LDPCDecoders

H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
decoder = BitFlipDecoder(H, 0.01, 100)

error = rand(1000) .< 0.01
syndrome = (H * error) .% 2

guess, converged = decode!(decoder, syndrome)
```

## Batch Decoding

```julia
samples = 100
errors = rand(1000, samples) .< 0.01
syndromes = H * errors .% 2

guesses, successes = batchdecode!(decoder, syndromes, zero(errors))
```

## API

See the [API Reference](../api.md) for full docstrings.