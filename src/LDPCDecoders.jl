module LDPCDecoders

using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using Statistics

using RowEchelon


# Write your package code here.
export
    parity_to_generator, hamming_to_parity, repetition_to_parity, 
    bp_decode, bp_simulate, it_decode, parity_check_matrix, save_pcm, 
    load_pcm, syndrome_decode, syndrome_simulate, syndrome_it_decode, 
    syndrome_it_simulate, syndrome_decode!, decode!, BeliefPropagationSetup, 
    BeliefPropagationDecoder, BeliefPropagationDecoderSparse, batchdecode!

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
include("decoders/belief_propogation.jl")
end
