module LDPCDecoders

# Write your package code here.
export parity_to_generator
include("generator.jl")
export hamming_to_parity
include("util.jl")
export repetition_to_parity
include("parity.jl")
export bp_decode
include("bp_decoder.jl")
export bp_simulate
include("bp_simulator.jl")
export it_decode
include("it_decoder.jl")
export parity_check_matrix, save_pcm, load_pcm
include("parity_generator.jl")
export syndrome_decode
include("syndrome_decoder.jl")
export syndrome_simulate
include("syndrome_simulator.jl")
export syndrome_it_decode
include("syndrome_it_decoder.jl")
export simulate_bp, plot_per_vs_ler
include("experiments.jl")
export syndrome_it_simulate
include("syndrome_it_simulate.jl")
export simulate_it
include("experiments_it.jl")
end
