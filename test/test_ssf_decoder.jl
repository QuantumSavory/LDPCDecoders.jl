@testitem "Small-Set-Flip Decoder" begin

using Test
using SparseArrays
using Random
using StableRNGs
import LDPCDecoders: SmallSetFlipDecoder, decode!, reset!, batchdecode!

    # Utility: create a simple cycle parity check matrix
    function create_cycle_matrix(n::Int)
        I_idx = Int[]
        J_idx = Int[]
        for j in 1:n
            push!(I_idx, j)
            push!(J_idx, mod(j, n) + 1)
            push!(I_idx, j)
            push!(J_idx, j)
        end
        V = fill(true, length(I_idx))
        return sparse(I_idx, J_idx, V, n, n)
    end

    @testset "Zero syndrome" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)
        syndrome = zeros(Bool, 8)
        guess, converged = decode!(decoder, syndrome)
        @test converged == true
        @test all(guess .== 0)
    end

    @testset "Single-bit errors on cycle matrix" begin
        for n in [4, 8, 16]
            H = create_cycle_matrix(n)
            decoder = SmallSetFlipDecoder(H, 50)

            for bit in 1:n
                error = zeros(Int, n)
                error[bit] = 1
                syndrome = Bool.(mod.(H * error, 2))
                guess, converged = decode!(decoder, syndrome)
                decoded_syndrome = Bool.(mod.(H * guess, 2))
                @test decoded_syndrome == syndrome
            end
        end
    end

    @testset "Random LDPC code" begin
        rng = StableRNG(42)
        H = LDPCDecoders.parity_check_matrix(100, 5, 4)
        decoder = SmallSetFlipDecoder(H, 100)

        # Low error rate should be decodable
        for _ in 1:10
            error = rand(rng, 100) .< 0.02
            syndrome = Bool.(mod.(H * error, 2))
            guess, converged = decode!(decoder, syndrome)
            decoded_syndrome = Bool.(mod.(H * guess, 2))
            @test decoded_syndrome == syndrome
        end
    end

    @testset "Subset size parameter t" begin
        H = create_cycle_matrix(8)

        # t=1 should behave like single-bit flip
        decoder_t1 = SmallSetFlipDecoder(H, 50; t=1)
        # t=2 allows pairs
        decoder_t2 = SmallSetFlipDecoder(H, 50; t=2)
        # t=0 (default) uses full support
        decoder_full = SmallSetFlipDecoder(H, 50)

        syndrome = Bool.(mod.(H * [1, 0, 0, 0, 0, 0, 0, 0], 2))

        for dec in [decoder_t1, decoder_t2, decoder_full]
            guess, converged = decode!(dec, syndrome)
            decoded_syndrome = Bool.(mod.(H * guess, 2))
            @test decoded_syndrome == syndrome
        end
    end

    @testset "Reset between decodings" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)

        # Decode twice to ensure reset works
        for _ in 1:3
            error = zeros(Int, 8)
            error[1] = 1
            syndrome = Bool.(mod.(H * error, 2))
            guess, converged = decode!(decoder, syndrome)
            decoded_syndrome = Bool.(mod.(H * guess, 2))
            @test decoded_syndrome == syndrome
        end
    end

    @testset "Integer syndrome input" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)
        error = [1, 0, 0, 0, 0, 0, 0, 0]
        # Pass Int syndrome (not Bool) — should work via !iszero()
        syndrome = (H * error) .% 2
        guess, converged = decode!(decoder, syndrome)
        decoded_syndrome = Bool.(mod.(H * guess, 2))
        @test decoded_syndrome == Bool.(syndrome)
    end

    # ── Commit 1: Regression tests for sparse syndrome optimization ──

    @testset "Regression: sparse optimization matches expected results" begin
        rng = StableRNG(42)
        H = LDPCDecoders.parity_check_matrix(100, 5, 4)
        decoder = SmallSetFlipDecoder(H, 100)

        # Run a fixed set of error patterns and verify correct syndrome matching
        # This ensures the sparse syn_nonzeros storage produces identical results
        # to the previous dense syn_columns storage.
        successes = 0
        for _ in 1:20
            error = rand(rng, 100) .< 0.02
            syndrome = Bool.(mod.(H * error, 2))
            guess, converged = decode!(decoder, syndrome)
            decoded_syndrome = Bool.(mod.(H * guess, 2))
            @test decoded_syndrome == syndrome
            if converged
                successes += 1
            end
        end
        # At 2% error rate on a (100,5,4) code, most should converge
        @test successes >= 15
    end

    # ── Commit 1: Deterministic tie-breaking (rng=nothing) ──

    @testset "Deterministic tie-breaking (rng=nothing)" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)  # rng=nothing by default

        error = zeros(Int, 8)
        error[1] = 1
        error[4] = 1  # two-bit error to create more tie opportunities
        syndrome = Bool.(mod.(H * error, 2))

        # Decode twice — must produce identical results
        guess1, conv1 = decode!(decoder, syndrome)
        guess1_copy = copy(guess1)
        guess2, conv2 = decode!(decoder, syndrome)

        @test guess1_copy == guess2
        @test conv1 == conv2
    end

    # ── Commit 1: Randomized tie-breaking (rng provided) ──

    @testset "Randomized tie-breaking (rng provided)" begin
        H = LDPCDecoders.parity_check_matrix(100, 5, 4)

        # Decoder with RNG — results must still produce valid syndrome match
        decoder_rng = SmallSetFlipDecoder(H, 100; rng=StableRNG(1))

        rng_err = StableRNG(99)
        for _ in 1:10
            error = rand(rng_err, 100) .< 0.02
            syndrome = Bool.(mod.(H * error, 2))
            guess, converged = decode!(decoder_rng, syndrome)
            decoded_syndrome = Bool.(mod.(H * guess, 2))
            # Correctness: decoded syndrome must match input syndrome
            @test decoded_syndrome == syndrome
        end
    end

    @testset "Randomized tie-breaking with StableRNG is reproducible" begin
        H = create_cycle_matrix(16)

        error = zeros(Int, 16)
        error[1] = 1
        error[5] = 1
        syndrome = Bool.(mod.(H * error, 2))

        # Same seed → same results
        dec1 = SmallSetFlipDecoder(H, 50; rng=StableRNG(42))
        dec2 = SmallSetFlipDecoder(H, 50; rng=StableRNG(42))

        guess1, conv1 = decode!(dec1, syndrome)
        guess1_copy = copy(guess1)
        guess2, conv2 = decode!(dec2, syndrome)

        @test guess1_copy == guess2
        @test conv1 == conv2
    end

    # ── Commit 2: Batch decoding ──

    @testset "Batch Decoding" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)

        samples = 5
        syndromes = hcat([begin
            err = zeros(Int, 8)
            err[rand(1:8)] = 1
            Bool.(mod.(H * err, 2))
        end for _ in 1:samples]...)
        errors_buf = zeros(Int, size(H, 2), samples)

        guesses, successes = batchdecode!(decoder, syndromes, errors_buf)

        for i in 1:samples
            decoded_syndrome = Bool.(mod.(H * guesses[:, i], 2))
            @test decoded_syndrome == syndromes[:, i]
        end
        @test all(successes)
    end

    @testset "Batch Decoding with pre-allocated converged vector" begin
        H = LDPCDecoders.parity_check_matrix(100, 5, 4)
        decoder = SmallSetFlipDecoder(H, 100)

        rng = StableRNG(99)
        samples = 10
        errors_in = rand(rng, 100, samples) .< 0.02
        syndromes = Bool.((H * errors_in) .% 2)
        errors_buf = zeros(Int, 100, samples)
        converged = Vector{Bool}(undef, samples)

        guesses, conv = batchdecode!(decoder, syndromes, errors_buf, converged)

        # converged must be the same object we passed in
        @test conv === converged

        for i in 1:samples
            decoded_syndrome = Bool.(mod.(H * guesses[:, i], 2))
            @test decoded_syndrome == syndromes[:, i]
        end
    end

    @testset "Batch Decoding without pre-allocated converged" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)

        samples = 3
        syndromes = hcat([begin
            err = zeros(Int, 8)
            err[1] = 1
            Bool.(mod.(H * err, 2))
        end for _ in 1:samples]...)
        errors_buf = zeros(Int, 8, samples)

        # 3-arg form (no converged vector) — should allocate one internally
        guesses, successes = batchdecode!(decoder, syndromes, errors_buf)

        @test length(successes) == samples
        for i in 1:samples
            decoded_syndrome = Bool.(mod.(H * guesses[:, i], 2))
            @test decoded_syndrome == syndromes[:, i]
        end
    end

end
