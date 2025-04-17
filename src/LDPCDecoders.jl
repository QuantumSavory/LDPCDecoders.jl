module LDPCDecoders

using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using Statistics

using RowEchelon


export
    decode!, batchdecode!,
    BeliefPropagationDecoder,
    BeliefPropagationOSDDecoder,
    BitFlipDecoder

include("generator.jl")
include("bp_decoder.jl")
include("bp_simulator.jl")
include("it_decoder.jl")
include("parity_generator.jl")

include("decoders/abstract_decoder.jl")
include("decoders/belief_propagation.jl")
include("decoders/belief_propagation_osd.jl")
include("decoders/iterative_bitflip.jl")
include("decoders/bpots_decoder.jl")

include("syndrome_bp_decoder.jl")
include("syndrome_simulator.jl")
include("syndrome_it_decoder.jl")
include("syndrome_it_simulate.jl")


end
