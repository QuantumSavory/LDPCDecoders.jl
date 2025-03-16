#BPOTS Tests

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
         for n in [4, 8, 16]
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
    generate_toric_code_matrix(; d::Int=3)
    Generate a parity check matrix for a toric code with distance `d`.
        Just use toric tht exists
        parity_checks_x(Toric(n,n))
        parity_checks_z(Toric(n,n))
        Can implement surface as well (edit(Toric))
    """
    function generate_toric_code_matrix(; d::Int=3)
        # Number of physical qubits is 2 * d^2
        n_qubits = 2 * d^2
        
        # Number of stabilizers (checks) is d^2 for X-type and d^2 for Z-type
        n_checks = 2 * d^2
        
        # Prepare for sparse matrix construction
        I_idx = Int[]  # Row indices
        J_idx = Int[]  # Column indices
        V = Bool[]     # Values
        
        # Create X-type stabilizers (plaquettes)
        for i in 0:(d-1)
            for j in 0:(d-1)
                check_idx = i*d + j + 1
                
                # Four edges of the plaquette
                qubit_h1 = i*d + j + 1
                qubit_h2 = i*d + ((j+1) % d) + 1
                qubit_v1 = d^2 + i*d + j + 1
                qubit_v2 = d^2 + ((i+1) % d)*d + j + 1
                
                # Add entries to sparse matrix
                for qubit in [qubit_h1, qubit_h2, qubit_v1, qubit_v2]
                    push!(I_idx, check_idx)
                    push!(J_idx, qubit)
                    push!(V, true)
                end
            end
        end
        
        # Create Z-type stabilizers (stars)
        for i in 0:(d-1)
            for j in 0:(d-1)
                check_idx = d^2 + i*d + j + 1
                
                # Four edges meeting at a vertex
                qubit_h1 = i*d + j + 1
                qubit_h2 = i*d + ((j-1+d) % d) + 1
                qubit_v1 = d^2 + i*d + j + 1
                qubit_v2 = d^2 + ((i-1+d) % d)*d + j + 1
                
                # Add entries to sparse matrix
                for qubit in [qubit_h1, qubit_h2, qubit_v1, qubit_v2]
                    push!(I_idx, check_idx)
                    push!(J_idx, qubit)
                    push!(V, true)
                end
            end
        end
        
        return sparse(I_idx, J_idx, V, n_checks, n_qubits)
    end

    """
    generate_surface_code_matrix(; d::Int=3)
    Generate a parity check matrix for a surface code with distance `d`.
    """
    function generate_surface_code_matrix(; d::Int=3)
        if d < 2
            error("Surface code distance must be at least 2")
        end
        
        # For a simple surface code, use a lattice with open boundaries
        # Number of qubits: (d-1)*d horizontal + d*(d-1) vertical
        n_qubits = 2 * d * (d-1)
        
        # Number of checks: (d-1)² plaquettes + d² vertices - 1 (dependent check)
        n_checks = (d-1)^2 + d^2 - 1
        
        # Prepare for sparse matrix construction
        I_idx = Int[]  # Row indices
        J_idx = Int[]  # Column indices
        V = Bool[]     # Values
        
        # Horizontal edge at (row, col)
        h_edge(row, col) = row*(d-1) + col + 1
        
        # Vertical edge at (row, col)
        v_edge(row, col) = (d-1)*d + row*d + col + 1
        
        # Create plaquette operators (X-type stabilizers)
        check_idx = 1
        for row in 0:(d-2)
            for col in 0:(d-2)
                edges = [
                    h_edge(row, col),      # Top
                    h_edge(row+1, col),    # Bottom
                    v_edge(row, col),      # Left
                    v_edge(row, col+1)     # Right
                ]
                
                for edge in edges
                    push!(I_idx, check_idx)
                    push!(J_idx, edge)
                    push!(V, true)
                end
                check_idx += 1
            end
        end
        
        # Create vertex operators (Z-type stabilizers)
        # Skip the last vertex as it's dependent on others
        for row in 0:d-1
            for col in 0:d-1
                # Skip the last vertex
                if row == d-1 && col == d-1
                    continue
                end
                
                # Collect edges connected to this vertex
                edges = Int[]
                
                # Horizontal edges
                if col > 0
                    push!(edges, h_edge(row, col-1))  # Left
                end
                if col < d-1
                    push!(edges, h_edge(row, col))    # Right
                end
                
                # Vertical edges
                if row > 0
                    push!(edges, v_edge(row-1, col))  # Top
                end
                if row < d-1
                    push!(edges, v_edge(row, col))    # Bottom
                end
                
                # Add entries to sparse matrix
                for edge in edges
                    push!(I_idx, check_idx)
                    push!(J_idx, edge)
                    push!(V, true)
                end
                check_idx += 1
            end
        end
        
        return sparse(I_idx, J_idx, V, n_checks, n_qubits)
    end

    @testset "Toric and Surface Code Performance" begin
        # Generate matrices with a smaller distance for testing
        d_test = 3
        println("Generating toric code matrix (d=$d_test)...")
        #H_toric = generate_toric_code_matrix(d=d_test)
        H_toric = parity_checks_x(Toric(d_test, d_test))
        
        println("Generating surface code matrix (d=$d_test)...")
        H_surface = generate_surface_code_matrix(d=d_test)
        
        # Test with reduced noise levels and iterations for faster testing
        noise_levels = [0.01, 0.05, 0.1]
        max_iter = 50
        
        for (matrix_idx, H) in enumerate([H_toric, H_surface])
            matrix_name = matrix_idx == 1 ? "Toric Code" : "Surface Code"
            println("\nTesting $matrix_name ($(size(H)[1])×$(size(H)[2]))...")
            
            performance_results = test_decoder_performance(
                (H, noise, max_iter) -> BPOTSDecoder(H, noise, max_iter; T=9, C=3.0),
                H, noise_levels, max_iter
            )
            
            for result in performance_results
                println("  $matrix_name @ Noise $(result.noise): Success Rate $(result.success_rate)")
                @test result.success_rate > 0.5
            end
        end
    end
     
end
