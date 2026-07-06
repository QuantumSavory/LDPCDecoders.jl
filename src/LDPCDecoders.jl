module LDPCDecoders

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
    BPOTSDecoder

include("parity_generator.jl")

include("decoders/abstract_decoder.jl")
include("decoders/belief_propagation.jl")
include("decoders/belief_propagation_osd.jl")
include("decoders/iterative_bitflip.jl")
include("decoders/bpots_decoder.jl")


end
