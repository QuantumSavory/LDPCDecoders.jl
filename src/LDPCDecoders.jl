module LDPCDecoders

using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using Statistics

using RowEchelon


export
    decode!, batchdecode!,
    BeliefPropagationDecoder

include("generator.jl")
include("parity.jl")
include("bp_decoder.jl")
include("bp_simulator.jl")
include("it_decoder.jl")
include("parity_generator.jl")

include("decoders/abstract_decoder.jl")
include("decoders/belief_propagation.jl")
include("syndrome_bp_decoder.jl")
include("syndrome_simulator.jl")
include("syndrome_it_decoder.jl")
include("syndrome_it_simulate.jl")


end
