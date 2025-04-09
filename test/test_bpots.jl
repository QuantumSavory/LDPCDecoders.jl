#BPOTS Tests

using Test
using SparseArrays
using Random
using QuantumClifford.ECC
import LDPCDecoders: BPOTSDecoder, BeliefPropagationDecoder, decode!, reset!

#= 
Tests syndrome matching
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
        
        performance_results = test_decoder_performance(
            (H, noise, max_iter) -> BPOTSDecoder(H, noise, max_iter; T=9, C=3.0),
            H, noise_levels, max_iter
        )
        
        for result in performance_results
            @test result.success_rate >= 0.9
        end
    end
end