#BPOTS Tests


#Plotting, can try 10000 or more iterations too
#Can plot error bars (1/\sqrt(sample))

using Test
using SparseArrays
using Random
using QuantumClifford.ECC
import LDPCDecoders: BPOTSDecoder, BeliefPropagationDecoder, decode!, reset!
 

#= 
Tests syndrome matching
create_cycle_matrix(n)- creates a cycle graph with known trapping set properties
generate_random_syndrome(H)- generates random error patterns (computes their syndromes using parity check matrices)
Remove the excess debug statements and make them more concise
=#
@testset "Comprehensive BP-OTS Decoder Tests" begin
     # Utility functions
     function create_cycle_matrix(n::Int)
         I_idx = Int[]
         J_idx = Int[]
         for j in 1:n
             push!(I_idx, j)
             push!(J_idx, mod(j,n)+1)
             push!(I_idx, j)
             push!(J_idx, j)
         end
         V = fill(true, length(I_idx))
         return sparse(I_idx, J_idx, V, n, n)
     end
 
     function generate_random_syndrome(H::SparseMatrixCSC{Bool})
         m, n = size(H)
         error = rand(Bool, n)
         return Vector{Bool}(mod.(H * error, 2))
     end
 
     function test_decoder_performance(decoder_type, H, noise_levels, max_iterations)
         results = []
         for noise in noise_levels
             success_count = 0
             total_runs = 100
 
             for _ in 1:total_runs
                 syndrome = generate_random_syndrome(H)
                 decoder = decoder_type(H, noise, max_iterations)
                 reset!(decoder)
                 result, converged = decode!(decoder, syndrome)
                 decoded_syndrome = Bool.(mod.(H * result, 2))
                 
                 if decoded_syndrome == syndrome
                     success_count += 1
                 end
             end
 
             push!(results, (noise=noise, success_rate=success_count/total_runs))
         end
         return results
     end
 
     @testset "Trapping Set Resistance" begin
         # Test with cycle matrices of different sizes
         for n in [4, 8, 16]   # [4, 8, 16]
             H = create_cycle_matrix(n)
             
             # Test BP-OTS with different configurations
             parameter_sets = [
                 (T=3, C=1.0),
                 (T=5, C=2.0),
                 (T=9, C=3.0)
             ]
             
             for (T, C) in parameter_sets
                 bpots_decoder = BPOTSDecoder(H, 0.01, 100; T=T, C=C)
                 reset!(bpots_decoder)
                 
                 # Creates a weight-2 error pattern (known to form trapping sets)
                 error = zeros(Bool, n)
                 error[1:2] .= true
                 syndrome = Vector{Bool}(mod.(H * error, 2))
                 
                 result, converged = decode!(bpots_decoder, syndrome)
                 decoded_syndrome = Bool.(mod.(H * result, 2))
                 
                 @test decoded_syndrome == syndrome
                 println("Cycle Matrix n=$n, T=$T, C=$C: Decoded successfully")
             end
         end
     end
 
     #############################################
     #Not rlly needed anymore
     #############################################
     @testset "Syndrome Matching Robustness" begin
         # Test BPOTS on different types of matrices
         test_matrices = [
             sparse([1,2,3,4], [1,2,3,4], true, 4, 4),  # Diagonal matrix
             create_cycle_matrix(4),                    # Cycle matrix
             sparse([1,1,2,2,3,3], [1,2,2,3,3,4], true, 4, 4)  # Irregular matrix
         ]
         
         noise_levels = [0.01, 0.1, 0.5]
         
         for H in test_matrices
             # Test BP-OTS performance
             performance_results = test_decoder_performance(
                 (H, noise, max_iter) -> BPOTSDecoder(H, noise, max_iter; T=9, C=3.0), 
                 H, noise_levels, 100
             )
             
             println("Performance for matrix:")
             for result in performance_results
                 println("  Noise: $(result.noise), Success Rate: $(result.success_rate)")
                 @test result.success_rate > 0.5  # Expect more than 50% success
             end
         end
     end
 
     @testset "Parameter Sensitivity" begin
         H = create_cycle_matrix(4)
         
         # Test T (threshold) values
         for T in [3, 5, 9, 15]
             bpots = BPOTSDecoder(H, 0.01, 100; T=T, C=3.0)
             reset!(bpots)
             
             syndrome = generate_random_syndrome(H)
             result, converged = decode!(bpots, syndrome)
             decoded_syndrome = Bool.(mod.(H * result, 2))
             
             println("T=$T - Converged: $converged, Syndrome match: $(decoded_syndrome == syndrome)")
             @test decoded_syndrome == syndrome
         end
         
         # Test C values
         for C in [1.0, 2.0, 5.0, 10.0]
             bpots = BPOTSDecoder(H, 0.01, 100; T=9, C=C)
             reset!(bpots)
             
             syndrome = generate_random_syndrome(H)
             result, converged = decode!(bpots, syndrome)
             decoded_syndrome = Bool.(mod.(H * result, 2))
             
             println("C=$C - Converged: $converged, Syndrome match: $(decoded_syndrome == syndrome)")
             @test decoded_syndrome == syndrome
         end
     end

    """
    generate_surface_code_matrix(; d::Int=3)
    Generate a parity check matrix for a surface code with distance `d`.
    """
    # function generate_surface_code_matrix(; d::Int=3)
    #     if d < 2
    #         error("Surface code distance must be at least 2")
    #     end
        
    #     # For a simple surface code, use a lattice with open boundaries
    #     # Number of qubits: (d-1)*d horizontal + d*(d-1) vertical
    #     n_qubits = 2 * d * (d-1)
        
    #     # Number of checks: (d-1)² plaquettes + d² vertices - 1 (dependent check)
    #     n_checks = (d-1)^2 + d^2 - 1
        
    #     # Prepare for sparse matrix construction
    #     I_idx = Int[]  # Row indices
    #     J_idx = Int[]  # Column indices
    #     V = Bool[]     # Values
        
    #     # Horizontal edge at (row, col)
    #     h_edge(row, col) = row*(d-1) + col + 1
        
    #     # Vertical edge at (row, col)
    #     v_edge(row, col) = (d-1)*d + row*d + col + 1
        
    #     # Create plaquette operators (X-type stabilizers)
    #     check_idx = 1
    #     for row in 0:(d-2)
    #         for col in 0:(d-2)
    #             edges = [
    #                 h_edge(row, col),      # Top
    #                 h_edge(row+1, col),    # Bottom
    #                 v_edge(row, col),      # Left
    #                 v_edge(row, col+1)     # Right
    #             ]
                
    #             for edge in edges
    #                 push!(I_idx, check_idx)
    #                 push!(J_idx, edge)
    #                 push!(V, true)
    #             end
    #             check_idx += 1
    #         end
    #     end
        
    #     # Create vertex operators (Z-type stabilizers)
    #     # Skip the last vertex as it's dependent on others
    #     for row in 0:d-1
    #         for col in 0:d-1
    #             # Skip the last vertex
    #             if row == d-1 && col == d-1
    #                 continue
    #             end
                
    #             # Collect edges connected to this vertex
    #             edges = Int[]
                
    #             # Horizontal edges
    #             if col > 0
    #                 push!(edges, h_edge(row, col-1))  # Left
    #             end
    #             if col < d-1
    #                 push!(edges, h_edge(row, col))    # Right
    #             end
                
    #             # Vertical edges
    #             if row > 0
    #                 push!(edges, v_edge(row-1, col))  # Top
    #             end
    #             if row < d-1
    #                 push!(edges, v_edge(row, col))    # Bottom
    #             end
                
    #             # Add entries to sparse matrix
    #             for edge in edges
    #                 push!(I_idx, check_idx)
    #                 push!(J_idx, edge)
    #                 push!(V, true)
    #             end
    #             check_idx += 1
    #         end
    #     end
        
    #     return sparse(I_idx, J_idx, V, n_checks, n_qubits)
    # end

    ##################################################
    #Test only for Toric Code, and make success_rate > 0.9
    ##################################################
    @testset "Toric Code Performance" begin
        # Generate toric code matrix with a smaller distance for testing
        d_test = 3
        println("Generating toric code matrix (d=$d_test)...")
        H_toric = parity_checks_x(Toric(d_test, d_test))
        
        # Test with reduced noise levels and iterations for faster testing
        noise_levels = [0.01, 0.05, 0.1]
        max_iter = 50
        
        matrix_name = "Toric Code"
        H = H_toric
        #println("\nTesting $matrix_name ($(size(H)[1])×$(size(H)[2]))...")
        
        performance_results = test_decoder_performance(
            (H, noise, max_iter) -> BPOTSDecoder(H, noise, max_iter; T=9, C=3.0),
            H, noise_levels, max_iter
        )
        
        for result in performance_results
            #println("  $matrix_name @ Noise $(result.noise): Success Rate $(result.success_rate)")
            @test result.success_rate >= 0.5
        end
    end
end






# OPTIONAL CODE FOR PLOTTING AND CHECKING RESULTS


# ######################

# using SparseArrays
# using Random
# using Plots
# using Measures
# using QuantumClifford.ECC
# import LDPCDecoders: BPOTSDecoder, BeliefPropagationDecoder, decode!, reset!

# # Generate a random error pattern with specified physical error rate
# # Generates random bit-flip errors with probability per for each qubit
# function generate_random_error(n::Int, per::Float64)
#     return Vector{Bool}(rand(n) .< per)
# end

# # Generate syndrome from error pattern
# # Calculates the corresponding syndrome by multiplying the parity check matrix with the error vector
# function generate_syndrome(H::SparseMatrixCSC{Bool}, error::Vector{Bool})
#     return Vector{Bool}(mod.(H * error, 2))
# end

# # Simple function to test a decoder on a single error pattern
# #=
# Takes a decoder instance and an error pattern
# Computes the syndrome from the error
# Runs the decoder on that syndrome
# Verifies if the decoder correctly matched the syndrome
# Returns whether the decoder succeeded and whether it converged
# =#
# function test_decoder(decoder, H::SparseMatrixCSC{Bool}, error::Vector{Bool})
#     syndrome = generate_syndrome(H, error)
#     result, converged = decode!(decoder, syndrome)
#     decoded_syndrome = Vector{Bool}(mod.(H * result, 2))
    
#     # Check if syndrome is matched correctly
#     syndrome_match = (decoded_syndrome == syndrome)
    
#     return syndrome_match, converged
# end

# # Function to evaluate decoder performance at different physical error rates
# #=
# Takes a function that constructs a decoder
# Tests it across multiple physical error rates
# For each error rate, runs many trials (1000 by default)
# Computes the logical error rate as proportion of failed decodings
# Handles the special case of zero errors (important for log plots)
# =#
# function evaluate_performance(decoder_constructor, H::SparseMatrixCSC{Bool}, 
#                              error_rates::Vector{Float64}, 
#                              trials_per_rate::Int=1000, max_iterations::Int=100)
    
#     m, n = size(H)
#     logical_error_rates = Float64[]
    
#     for error_rate in error_rates
#         println("Testing error rate: $error_rate")
#         logical_errors = 0
        
#         for trial in 1:trials_per_rate
#             # Generate random error
#             error = generate_random_error(n, error_rate)
            
#             # Create and reset decoder
#             decoder = decoder_constructor(H, error_rate, max_iterations)
#             reset!(decoder)
            
#             # Test decoder
#             success, _ = test_decoder(decoder, H, error)
            
#             if !success
#                 logical_errors += 1
#             end
            
#             if trial % 100 == 0
#                 println("  Completed $trial/$trials_per_rate trials")
#             end
#         end
        
#         # Calculate logical error rate
#         logical_error_rate = logical_errors / trials_per_rate
        
#         # Handle zero error rate for log plots (set to a very small value)
#         if logical_error_rate == 0
#             logical_error_rate = 1.0 / (10 * trials_per_rate)  # Set to smaller than could be detected
#             println("  Logical error rate: 0 (adjusted to $logical_error_rate for plotting)")
#         else
#             println("  Logical error rate: $logical_error_rate")
#         end
        
#         push!(logical_error_rates, logical_error_rate)
#     end
    
#     return logical_error_rates
# end

# # Compare BP and BP-OTS decoders
# #=
# Takes a code and parameters for testing
# Evaluates both BP and BP-OTS decoders
# Creates a plot comparing their performance
# Returns both the plot and the raw data
# The plotting code uses multiple formatting options to create publication-quality graphs:
# Logarithmic scales for both axes
# Different markers for each decoder (circles for BP, squares for BP-OTS)
# Grid lines for easier reading
# Appropriate coloring (blue for BP, red for BP-OTS)
# =#
# function compare_bp_vs_bpots(code_name::String, H::SparseMatrixCSC{Bool}, 
#                             error_rates::Vector{Float64}, 
#                             trials_per_rate::Int=500, max_iterations::Int=100,
#                             bpots_T::Int=9, bpots_C::Float64=3.0)
    
#     println("Evaluating BP decoder...")
#     bp_constructor = (H, error_rate, max_iter) -> BeliefPropagationDecoder(H, error_rate, max_iter)
#     bp_logical_rates = evaluate_performance(bp_constructor, H, error_rates, trials_per_rate, max_iterations)
    
#     println("Evaluating BP-OTS decoder...")
#     bpots_constructor = (H, error_rate, max_iter) -> 
#         BPOTSDecoder(H, error_rate, max_iter; T=bpots_T, C=bpots_C)
#     bpots_logical_rates = evaluate_performance(bpots_constructor, H, error_rates, trials_per_rate, max_iterations)
    
#     # Fixed figure size
#     default(size=(800, 600), framestyle=:box, margin=10mm)
    
#     # Create plot with improved formatting
#     p = plot(
#         error_rates, bp_logical_rates,
#         label="BP",
#         marker=:circle,
#         markersize=6,
#         linewidth=3,
#         xscale=:log10,
#         yscale=:log10,
#         xlabel="Physical Error Rate",
#         ylabel="Logical Error Rate",
#         title="$code_name: BP vs BP-OTS Comparison",
#         legend=:topleft,
#         grid=true,
#         minorgrid=true,
#         fontfamily="Computer Modern",
#         dpi=300,
#         lw=2
#     )
    
#     plot!(
#         p, error_rates, bpots_logical_rates,
#         label="BP-OTS (T=$bpots_T, C=$bpots_C)",
#         marker=:square,
#         markersize=6,
#         linewidth=3,
#         color=:red
#     )
    
#     # Save the plot
#     savefig(p, "$(lowercase(replace(code_name, " " => "_")))_comparison.png")
    
#     return p, (bp_logical_rates, bpots_logical_rates)
# end

# # Function to run experiments with multiple code distances
# #=
# Tests multiple code distances (3 and 5)
# For each distance, creates the appropriate toric code
# Runs the comparison between decoders
# Creates individual and combined plots
# =#
# function run_toric_comparison()
#     # Parameters
#     error_rates = [0.01, 0.02, 0.05, 0.08, 0.1, 0.15]
#     trials_per_rate = 10000  # More trials for better statistics
#     max_iterations = 100
    
#     # Test different distances
#     distances = [3, 5]  # Add more if computational resources allow
    
#     plots = []
#     results = Dict()
    
#     for d in distances
#         println("\n===== Testing Toric Code (d=$d) =====")
#         H_toric = parity_checks_x(Toric(d, d))
#         p, r = compare_bp_vs_bpots(
#             "Toric Code (d=$d)", 
#             H_toric, 
#             error_rates, 
#             trials_per_rate, 
#             max_iterations
#         )
#         push!(plots, p)
#         results[d] = r
#     end
    
#     # Create combined plot
#     combined_plot = plot(plots..., layout=(length(plots), 1), size=(800, 500*length(plots)))
#     savefig(combined_plot, "toric_code_all_distances.png")
    
#     return combined_plot, results
# end

# # Run for Toric codes
# run_toric_comparison()

# ######################