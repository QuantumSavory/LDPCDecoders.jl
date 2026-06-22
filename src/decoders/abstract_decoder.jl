"""
    AbstractDecoder

Abstract supertype for all LDPC decoders in this package.

Concrete decoders (e.g., [`BeliefPropagationDecoder`](@ref), [`BitFlipDecoder`](@ref))
subtype `AbstractDecoder` and must implement [`decode!`](@ref).
"""
abstract type AbstractDecoder end
