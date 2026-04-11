# LDPCDecoders.jl

A Julia package providing iterative decoders for Low-Density Parity-Check (LDPC) codes, with a focus on quantum error correction applications.

## Overview

LDPCDecoders.jl implements several classical LDPC decoding algorithms that can be used standalone or bridged into [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl) via the `QuantumCliffordLDPCDecodersExt` package extension for quantum syndrome decoding and benchmarking.

### Available Decoders

| Decoder | Type | Description |
|:--------|:-----|:------------|
| [`BeliefPropagationDecoder`](@ref) | Soft-decision | Standard sum-product message passing on the Tanner graph |
| [`BeliefPropagationOSDDecoder`](@ref) | Soft-decision | Belief propagation with Ordered Statistics Decoding post-processing |
| [`BitFlipDecoder`](@ref) | Hard-decision | Classical greedy bit-flip algorithm |
| [`BPOTSDecoder`](@ref) | Soft-decision | Belief propagation with Ordered Trapping Set biasing |

### Architecture

All decoders follow a common pattern:

1. **`AbstractDecoder`** — Every decoder subtypes this abstract type.
2. **`ScratchSpace`** — Each decoder pre-allocates mutable working arrays to prevent memory allocations during the decoding loop.
3. **`decode!(decoder, syndrome)`** — The mutating entry point that decodes a syndrome vector in-place using the pre-allocated scratch space.
4. **`batchdecode!(decoder, syndromes, errors)`** — Batch variant that re-uses the same scratch space across multiple decoding trials.

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
