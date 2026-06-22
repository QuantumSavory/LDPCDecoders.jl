# API Reference

This page lists the complete public API of LDPCDecoders.jl.

## Types

```@docs
AbstractDecoder
BeliefPropagationDecoder
BeliefPropagationOSDDecoder
BitFlipDecoder
BPOTSDecoder
```

## Functions

```@docs
decode!
batchdecode!
```

## Internal Utilities

These are not exported but may be useful for advanced usage.

```@docs
LDPCDecoders.reset!(::BeliefPropagationDecoder)
LDPCDecoders.reset!(::BitFlipDecoder)
LDPCDecoders.parity_check_matrix
```
