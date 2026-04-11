"""
    AbstractDecoder

Abstract supertype for all LDPC decoders in this package.

Every concrete decoder (e.g., [`BeliefPropagationDecoder`](@ref), [`BitFlipDecoder`](@ref))
subtypes `AbstractDecoder` and implements the [`decode!`](@ref) interface.
"""
abstract type AbstractDecoder end
