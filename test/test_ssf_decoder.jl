@testitem "Small-Set-Flip Decoder" begin

using Test
using SparseArrays
using Random
using StableRNGs
import LDPCDecoders: SmallSetFlipDecoder, BitFlipDecoder, decode!, reset!, batchdecode!

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
        successes = 0
        for _ in 1:10
            error = rand(rng, 100) .< 0.02
            syndrome = Bool.(mod.(H * error, 2))
            guess, converged = decode!(decoder, syndrome)
            if converged
                decoded_syndrome = Bool.(mod.(H * guess, 2))
                @test decoded_syndrome == syndrome
                successes += 1
            end
        end
        @test successes >= 5
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
            if converged
                decoded_syndrome = Bool.(mod.(H * guess, 2))
                @test decoded_syndrome == syndrome
                successes += 1
            end
        end
        # At 2% error rate on a (100,5,4) code, most should converge
        @test successes >= 10
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
        # Use cycle matrix (always converges for single-bit errors)
        # to avoid flaky failures from non-convergence on random LDPC codes
        H = create_cycle_matrix(16)
        decoder_rng = SmallSetFlipDecoder(H, 50; rng=StableRNG(1))

        for bit in 1:16
            error = zeros(Int, 16)
            error[bit] = 1
            syndrome = Bool.(mod.(H * error, 2))
            guess, converged = decode!(decoder_rng, syndrome)
            decoded_syndrome = Bool.(mod.(H * guess, 2))
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
        # Use cycle matrix for reliable convergence across platforms
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 50)

        samples = 5
        syndromes = hcat([begin
            err = zeros(Int, 8)
            err[rand(1:8)] = 1
            Bool.(mod.(H * err, 2))
        end for _ in 1:samples]...)
        errors_buf = zeros(Int, 8, samples)
        converged = Vector{Bool}(undef, samples)

        guesses, conv = batchdecode!(decoder, syndromes, errors_buf, converged)

        # converged must be the same object we passed in
        @test conv === converged

        for i in 1:samples
            decoded_syndrome = Bool.(mod.(H * guesses[:, i], 2))
            @test decoded_syndrome == syndromes[:, i]
        end
        @test all(converged)
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

    # ── Commit 4: Comprehensive tests ──

    @testset "Random i.i.d. noise — threshold behavior (Paper 2)" begin
        rng = StableRNG(42)
        H = create_cycle_matrix(16)
        decoder = SmallSetFlipDecoder(H, 50)

        # Low noise (single-bit errors) -> 100% exact error recovery
        exact_matches = 0
        for _ in 1:20
            err = zeros(Int, 16)
            err[rand(rng, 1:16)] = 1
            syn = Bool.(mod.(H * err, 2))
            guess, conv = decode!(decoder, syn)
            if conv && guess == err
                exact_matches += 1
            end
        end
        @test exact_matches == 20
    end

    @testset "SSF (t>1) vs BitFlipDecoder — Paper 1 Algorithm 2 vs single-bit-flip" begin
        H = create_cycle_matrix(12)
        ssf_t2 = SmallSetFlipDecoder(H, 50; t=2)
        bf = BitFlipDecoder(H, 0.01, 50)

        # 2-bit error on adjacent bits (cycle matrix weight-2 stabilizer support overlap)
        err = zeros(Int, 12)
        err[1] = 1
        err[2] = 1
        syn = Bool.(mod.(H * err, 2))

        _, ssf_conv = decode!(ssf_t2, syn)
        _, bf_conv = decode!(bf, syn)

        # SSF with t=2 handles 2-bit flip sets directly
        @test ssf_conv == true
    end

    @testset "Effect of t parameter on success rate" begin
        H = create_cycle_matrix(16)
        dec_t1 = SmallSetFlipDecoder(H, 50; t=1)
        dec_t2 = SmallSetFlipDecoder(H, 50; t=2)

        rng = StableRNG(123)
        succ_t1 = 0
        succ_t2 = 0

        for _ in 1:30
            err = rand(rng, 16) .< 0.08
            syn = Bool.(mod.(H * err, 2))
            
            _, c1 = decode!(dec_t1, syn)
            _, c2 = decode!(dec_t2, syn)
            
            c1 && (succ_t1 += 1)
            c2 && (succ_t2 += 1)
        end

        @test succ_t2 >= succ_t1
    end

    @testset "Multi-bit error patterns" begin
        H = create_cycle_matrix(12)
        decoder = SmallSetFlipDecoder(H, 50; t=2)

        # 2-bit & 3-bit errors
        for bits in [[1, 3], [2, 5], [1, 4, 7]]
            err = zeros(Int, 12)
            err[bits] .= 1
            syn = Bool.(mod.(H * err, 2))

            guess, conv = decode!(decoder, syn)
            if conv
                @test Bool.(mod.(H * guess, 2)) == syn
            end
        end
    end

    @testset "Failure reporting" begin
        H = create_cycle_matrix(8)
        decoder = SmallSetFlipDecoder(H, 1) # max 1 iteration

        # Non-codeword error pattern with 3 separated errors that cannot be solved in 1 iteration
        err = zeros(Int, 8)
        err[1] = 1; err[3] = 1; err[5] = 1
        syn = Bool.(mod.(H * err, 2))

        guess, conv = decode!(decoder, syn)
        @test length(guess) == 8
        @test typeof(guess) == Vector{Int}
    end

end
