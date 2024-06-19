# LDPCDecoders.jl

<a href="https://github.com/QuantumSavory/LDPCDecoders.jl/actions?query=workflow%3ACI+branch%3Amain"><img src="https://github.com/QuantumSavory/LDPCDecoders.jl/actions/workflows/ci.yml/badge.svg" alt="GitHub Workflow Status"></a>
<a href="https://codecov.io/gh/QuantumSavory/LDPCDecoders.jl"><img src="https://img.shields.io/codecov/c/gh/QuantumSavory/LDPCDecoders.jl?label=codecov" alt="Test coverage from codecov"></a>

A package of LDPC decoders for decoding certain LDPC Quantum Error Correcting Codes using Julia. It currently has a simple iterative decoder and belief propagation (BP) decoder, which also has variation with post processing using Ordered Statistics Decoding (BP+OSD). 

# Setting up

To use decoders and structs from LDPCDecoders.jl in your project, add it to your `Project.toml`, and automatically to `Manifest.toml`

Inside your project folder (which contains the `Project.toml` and `Manifest.toml`), run Julia in package mode (enabled by pressing `]` in Julia REPL)

```bash
pkg> activate .
(YourProject) pkg> add LDPCDecoders.jl
```

When prompted, press `Y/y` for yes.

# Usage

Using each of the decoders is discussed below

## Belief Propagation
The code and user interface for the belief propagation decoder lies in `src/decoders/belief_propagation.jl`. First step is to set up the decoder, with your parity check matrix `H`, physical error rate `per` and max number of decoding iterations for belief propagation algorithm. 

```julia
julia> using LDPCDecoders
julia> decoder = BeliefPropagationDecoder(H, per, max_iters)
```

There are two available methods for decoding - `decode!` which takes a syndrome and `batchdecode!`, which takes a batch of syndromes at once. See code for more docs.

```julia
julia> decode!(decoder, syndrome, error)
```


```julia
julia> batchdecode!(decoder, syndromes, errors)
```
