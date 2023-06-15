module LDPC

# Write your package code here.
export parity_to_generator
include("generator.jl")
export hamming_to_parity
include("util.jl")
export repetition_to_parity
include("parity.jl")
end
