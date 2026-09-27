@testitem "CSS Small-Set-Flip Decoder" begin

using Test
using SparseArrays
using Random
using StableRNGs
import LDPCDecoders: CSSSmallSetFlipDecoder, SmallSetFlipDecoder, decode!, reset!, batchdecode!

    function create_cycle_matrix(n::Int)
        I_idx = Int[]
        J_idx = Int[]
        for j in 1:n
            push!(I_idx, j)
            push!(J_idx, mod(j, n) + 1)
            push!(I_idx, j)
            push!(J_idx, j)
        end
        return sparse(I_idx, J_idx, fill(true, length(I_idx)), n, n)
    end

    @testset "CSS Decoder construction & basic decoding" begin
        Hx = create_cycle_matrix(8)
        Hz = create_cycle_matrix(8)

        decoder = CSSSmallSetFlipDecoder(Hx, Hz, 50; t=2)

        # X error only
        err_x = zeros(Int, 8)
        err_x[1] = 1
        syn_z = Bool.(mod.(Hz * err_x, 2))
        syn_x = zeros(Bool, 8)

        (guess_x, guess_z), conv = decode!(decoder, syn_x, syn_z)
        @test conv
        @test Bool.(mod.(Hz * guess_x, 2)) == syn_z
        @test Bool.(mod.(Hx * guess_z, 2)) == syn_x
    end

    @testset "Combined syndrome vector" begin
        Hx = create_cycle_matrix(8)
        Hz = create_cycle_matrix(8)
        decoder = CSSSmallSetFlipDecoder(Hx, Hz, 50)

        err_x = zeros(Int, 8); err_x[2] = 1
        err_z = zeros(Int, 8); err_z[5] = 1

        syn_x = Bool.(mod.(Hx * err_z, 2))
        syn_z = Bool.(mod.(Hz * err_x, 2))
        combined_syn = vcat(syn_x, syn_z)

        guess_combined, conv = decode!(decoder, combined_syn)
        @test conv
        @test length(guess_combined) == 16
        
        guess_x = guess_combined[1:8]
        guess_z = guess_combined[9:16]
        @test Bool.(mod.(Hz * guess_x, 2)) == syn_z
        @test Bool.(mod.(Hx * guess_z, 2)) == syn_x
    end

    @testset "Reset & batchdecode" begin
        Hx = create_cycle_matrix(8)
        Hz = create_cycle_matrix(8)
        decoder = CSSSmallSetFlipDecoder(Hx, Hz, 50)

        reset!(decoder)

        samples = 3
        syndromes = hcat([begin
            ex = zeros(Int, 8); ex[1] = 1
            ez = zeros(Int, 8); ez[2] = 1
            vcat(Bool.(mod.(Hx * ez, 2)), Bool.(mod.(Hz * ex, 2)))
        end for _ in 1:samples]...)

        errors_buf = zeros(Int, 16, samples)
        guesses, successes = batchdecode!(decoder, syndromes, errors_buf)

        @test all(successes)
        for i in 1:samples
            gx = guesses[1:8, i]
            gz = guesses[9:16, i]
            @test Bool.(mod.(Hz * gx, 2)) == syndromes[9:16, i]
            @test Bool.(mod.(Hx * gz, 2)) == syndromes[1:8, i]
        end
    end

end
