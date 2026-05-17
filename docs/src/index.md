# LDPCDecoders.jl

Iterative decoders for Low-Density Parity-Check (LDPC) codes in Julia, aimed at
quantum error correction applications.

## Overview

LDPCDecoders.jl provides several classical LDPC decoding algorithms, usable standalone or plugged into [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl) via the `QuantumCliffordLDPCDecodersExt` extension for quantum syndrome decoding and benchmarking.

### Available Decoders

| Decoder | Type | Description |
|:--------|:-----|:------------|
| [`BeliefPropagationDecoder`](@ref) | Soft-decision | Standard sum-product message passing on the Tanner graph |
| [`BeliefPropagationOSDDecoder`](@ref) | Soft-decision | Belief propagation with Ordered Statistics Decoding post-processing |
| [`BitFlipDecoder`](@ref) | Hard-decision | Classical greedy bit-flip algorithm |
| [`BPOTSDecoder`](@ref) | Soft-decision | Belief propagation with Ordered Trapping Set biasing |

### Architecture

All decoders follow the same structure:

1. **`AbstractDecoder`** — Every decoder subtypes this abstract type.
2. **`ScratchSpace`** — Each decoder preallocates mutable working arrays to avoid
   allocations during the decoding loop.
3. **`decode!(decoder, syndrome)`** — Decodes a syndrome vector in-place using the
   preallocated scratch space.
4. **`batchdecode!(decoder, syndromes, errors)`** — Batch variant that reuses the same
   scratch space across calls.

## Quick Start

```julia
using LDPCDecoders

# Generate a random LDPC parity check matrix
H = LDPCDecoders.parity_check_matrix(1000, 10, 9)

# Create a decoder
decoder = BeliefPropagationDecoder(H, 0.01, 100)

# Simulate an error and compute its syndrome
error = rand(1000) .< 0.01
syndrome = (H * error) .% 2

# Decode
guess, converged = decode!(decoder, syndrome)
```

## Installation

```julia
using Pkg
Pkg.add("LDPCDecoders")
```