module LDPCDecoders

using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using Statistics

using RowEchelon


export
    parity_check_matrix, decode!,
    BeliefPropagationDecoder, batchdecode!

include("generator.jl")
include("util.jl")
include("parity.jl")
include("bp_decoder.jl")
include("bp_simulator.jl")
include("it_decoder.jl")
include("parity_generator.jl")
include("syndrome_decoder.jl")
include("syndrome_simulator.jl")
include("syndrome_it_decoder.jl")
include("syndrome_it_simulate.jl")
end
