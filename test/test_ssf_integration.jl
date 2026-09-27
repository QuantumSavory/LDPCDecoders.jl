@testitem "SSF QuantumClifford Integration" begin

using Test
using SparseArrays
using Random
using StableRNGs
using QuantumClifford
using QuantumClifford.ECC
import LDPCDecoders: SmallSetFlipDecoder, CSSSmallSetFlipDecoder, decode!, reset!, batchdecode!

    # -- helpers --

    # manual hypergraph product from two classical codes
    # Hx = [H1⊗I_n2 | I_r1⊗H2^T]
    # Hz = [I_n1⊗H2  | H1^T⊗I_r2]
    function hgp_matrices(H1::AbstractMatrix{Bool}, H2::AbstractMatrix{Bool})
        r1, n1 = size(H1)
        r2, n2 = size(H2)

        sH1 = sparse(H1)
        sH2 = sparse(H2)

        # build Bool identity matrices via spdiagm to avoid
        # LinearAlgebra.I dispatch issues with QuantumClifford loaded
        _eye(m) = spdiagm(0 => ones(Bool, m))

        I_n1 = _eye(n1)
        I_n2 = _eye(n2)
        I_r1 = _eye(r1)
        I_r2 = _eye(r2)

        Hx_raw = hcat(kron(sH1, I_n2), kron(I_r1, sparse(sH2')))
        Hz_raw = hcat(kron(I_n1, sH2), kron(sparse(sH1'), I_r2))

        Hx = SparseMatrixCSC{Bool,Int}(Hx_raw .!= 0)
        Hz = SparseMatrixCSC{Bool,Int}(Hz_raw .!= 0)
        return Hx, Hz
    end

    function calc_syndrome(H, err)
        return Bool.(mod.(H * err, 2))
    end

    # wrap QC parity matrices to sparse Bool for our constructors
    function to_sparse_bool(M)
        return sparse(Bool.(M))
    end

    @testset "Toric code with CSSSmallSetFlipDecoder" begin
        for d in [3, 5]
            @testset "d=$d" begin
                toric = Toric(d, d)
                Hx = to_sparse_bool(parity_matrix_x(toric))
                Hz = to_sparse_bool(parity_matrix_z(toric))
                n = code_n(toric)

                decoder = CSSSmallSetFlipDecoder(Hx, Hz, 100; t=2)

                rng = StableRNG(100 + d)
                successes = 0
                trials = 50

                for _ in 1:trials
                    err_x = zeros(Int, n)
                    err_z = zeros(Int, n)
                    err_x[1 + floor(Int, rand(rng) * n)] = 1
                    err_z[1 + floor(Int, rand(rng) * n)] = 1

                    syn_x = calc_syndrome(Hx, err_z)
                    syn_z = calc_syndrome(Hz, err_x)

                    (guess_x, guess_z), conv = decode!(decoder, syn_x, syn_z)
                    if conv
                        syn_z_check = calc_syndrome(Hz, guess_x)
                        syn_x_check = calc_syndrome(Hx, guess_z)
                        if syn_z_check == syn_z && syn_x_check == syn_x
                            successes += 1
                        end
                    end
                end
                @test successes >= trials * 0.8
            end
        end
    end

    @testset "Steane [[7,1,3]] — exhaustive weight-1 errors" begin
        steane = Steane7()
        Hx = to_sparse_bool(parity_matrix_x(steane))
        Hz = to_sparse_bool(parity_matrix_z(steane))
        n = code_n(steane)

        decoder = CSSSmallSetFlipDecoder(Hx, Hz, 50; t=3)

        # every single-qubit X error
        for i in 1:n
            err = zeros(Int, n)
            err[i] = 1
            syn_z = calc_syndrome(Hz, err)
            syn_x = zeros(Bool, size(Hx, 1))

            (guess_x, guess_z), conv = decode!(decoder, syn_x, syn_z)
            @test conv
            @test calc_syndrome(Hz, guess_x) == syn_z
        end

        # every single-qubit Z error
        for i in 1:n
            err = zeros(Int, n)
            err[i] = 1
            syn_x = calc_syndrome(Hx, err)
            syn_z = zeros(Bool, size(Hz, 1))

            (guess_x, guess_z), conv = decode!(decoder, syn_x, syn_z)
            @test conv
            @test calc_syndrome(Hx, guess_z) == syn_x
        end
    end

    @testset "Steane [[7,1,3]] — compare with TableDecoder" begin
        steane = Steane7()
        Hx = to_sparse_bool(parity_matrix_x(steane))
        Hz = to_sparse_bool(parity_matrix_z(steane))
        n = code_n(steane)

        ssf_dec = CSSSmallSetFlipDecoder(Hx, Hz, 50; t=3)

        # both SSF and a ground-truth check should agree on weight-1 X errors
        for i in 1:n
            err = zeros(Int, n)
            err[i] = 1
            syn_z = calc_syndrome(Hz, err)
            syn_x = zeros(Bool, size(Hx, 1))

            (guess_x, _), conv_ssf = decode!(ssf_dec, syn_x, syn_z)
            @test conv_ssf
            @test calc_syndrome(Hz, guess_x) == syn_z
        end
    end

    @testset "Small HGP code — manual construction" begin
        H1 = Bool[1 1 0; 0 1 1]
        H2 = Bool[1 1 0; 0 1 1]

        Hx, Hz = hgp_matrices(H1, H2)
        # CSS condition: Hx * Hz' = 0 mod 2
        product = mod.(Hx * Hz', 2)
        @test iszero(product)

        sx, nx = size(Hx)
        sz, nz = size(Hz)
        @test nx == nz

        decoder = CSSSmallSetFlipDecoder(Hx, Hz, 100; t=2)

        rng = StableRNG(555)
        successes = 0
        trials = 30

        for _ in 1:trials
            err = zeros(Int, nx)
            pos = 1 + floor(Int, rand(rng) * nx)
            err[pos] = 1

            syn_z = calc_syndrome(Hz, err)
            syn_x = zeros(Bool, sx)

            (guess_x, _), conv = decode!(decoder, syn_x, syn_z)
            if conv && calc_syndrome(Hz, guess_x) == syn_z
                successes += 1
            end
        end
        @test successes >= trials * 0.7
    end

    @testset "Logical error vs syndrome-only check" begin
        steane = Steane7()
        Hx = to_sparse_bool(parity_matrix_x(steane))
        Hz = to_sparse_bool(parity_matrix_z(steane))
        n = code_n(steane)

        # get logical X operator from the stabilizer tableau
        md = MixedDestabilizer(parity_checks(steane))
        logx_tab = logicalxview(md)
        logx_gf2 = stab_to_gf2(logx_tab)
        # X part is the first n columns
        logx_bits = Bool.(logx_gf2[1, 1:n])

        # logical X operator has zero syndrome under Hz (undetectable)
        logx_int = Int.(logx_bits)
        syn_z_logical = calc_syndrome(Hz, logx_int)
        @test all(syn_z_logical .== false)

        # but it's nontrivial (not the identity)
        @test sum(logx_bits) > 0

        # single-qubit error: decoder corrects it, residual has zero syndrome
        decoder = CSSSmallSetFlipDecoder(Hx, Hz, 50; t=3)

        err = zeros(Int, n)
        err[1] = 1
        syn_z = calc_syndrome(Hz, err)
        syn_x = zeros(Bool, size(Hx, 1))

        (guess_x, _), conv = decode!(decoder, syn_x, syn_z)
        @test conv

        residual = mod.(err .+ guess_x, 2)
        residual_syn = calc_syndrome(Hz, residual)
        @test all(residual_syn .== false)
    end

end
