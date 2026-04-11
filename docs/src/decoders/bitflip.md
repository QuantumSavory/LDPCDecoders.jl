# Iterative BitFlip Decoder

The [`BitFlipDecoder`](@ref) implements the classical hard-decision bit-flip decoding algorithm. At each iteration, it computes a syndrome from the current error estimate, counts "votes" from unsatisfied checks for each bit, and flips the bit with the most votes. This process repeats until the syndrome matches or no further improvement is possible.

BitFlip decoding is the simplest iterative decoder and serves as a useful baseline for benchmarking more sophisticated algorithms.

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

