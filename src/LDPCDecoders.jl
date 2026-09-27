module LDPCDecoders

using Combinatorics: combinations
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using Statistics

using RowEchelon


export
    decode!, batchdecode!,
    AbstractDecoder,
    BeliefPropagationDecoder,
    BeliefPropagationOSDDecoder,
    BitFlipDecoder,
    BPOTSDecoder,
    SmallSetFlipDecoder,
    CSSSmallSetFlipDecoder

include("parity_generator.jl")

include("decoders/abstract_decoder.jl")
include("decoders/belief_propagation.jl")
include("decoders/belief_propagation_osd.jl")
include("decoders/iterative_bitflip.jl")
include("decoders/bpots_decoder.jl")
include("decoders/small_set_flip.jl")
include("decoders/small_set_flip_css.jl")


end
