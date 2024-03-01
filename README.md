# LDPCDecoders.jl

[![Build Status](https://github.com/krishna-praneet/LDPC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/krishna-praneet/LDPC.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package of LDPC decoders for decoding certain LDPC Quantum Error Correcting Codes using Julia. It currently has a simple iterative decoder and belief propagation (BP) decoder, which also has variation with post processing using Ordered Statistics Decoding (BP+OSD). 

# Setting up

To use decoders and structs from LDPCDecoders.jl in your project, add it to your `Project.toml`, and automatically to `Manifest.toml`

Inside your project folder (which contains the `Project.toml` and `Manifest.toml`), run julia in package mode (enabled by pressing `]` in Julia REPL)

```bash
pkg> activate .
(YourProject) pkg> add LDPCDecoders.jl
```

When prompted, press `Y/y` for yes.

# Usage

Using each of the decodersis discussed below

## Belief Propagation
The code and user interface for belied propagation decoder lies in `src/decoders/belief_propagation.jl`. First step is to set up the decoder, with your parity check matrix `H`, physical error rate `per` and max number of decoding iterations for belied propagation algorithm. 

```julia
julia> using LDPCDecoders
julia> decoder = BeliefPropagationDeocder(H, per, max_iters)
```

There are two available methods for decoding - `decode!` which takes a syndrome and `batchdecode!`, which takes a batch of syndromes at once. See code for more docs.

```julia
julia> decode!(decoder, syndrome, error)
```


```julia
julia> batchdecode!(decoder, syndromes, errors)
```
