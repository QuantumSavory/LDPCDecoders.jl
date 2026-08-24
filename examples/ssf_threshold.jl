# SSF Threshold Estimation Script
#
# Reproduces the Monte Carlo threshold estimation from:
#   Grospellier & Krishna, "Numerical study of hypergraph product codes"
#   arXiv:1810.03681 (2019)
#
# Estimates the logical error rate of the SSF decoder across a sweep of
# physical error rates and code sizes. The threshold p_th is visible as
# the crossing point of the resulting curves.
#
# Usage:
#   julia --project=. examples/ssf_threshold.jl
#   julia --project=. examples/ssf_threshold.jl --samples 5000 --outfile results.csv
#
# The script constructs Toric-code-like CSS codes of increasing distance
# (via manual HGP construction from cyclic repetition codes), runs SSF
# decoding under i.i.d. depolarizing noise, and records logical error rates.
#
# For HGP codes from biregular random graphs (Paper 3 families), one would
# need the full QuantumClifford.ECC.GeneralizedHyperGraphProduct constructor
# with Hecke.jl — here we use the simpler Toric-code pathway as a
# self-contained demonstration.

using LDPCDecoders
using SparseArrays
using LinearAlgebra
using Random
using Printf
using DelimitedFiles

# ============================================================================
# CLI argument parsing
# ============================================================================

function parse_args(args)
    defaults = Dict(
        :samples   => 1000,
        :outfile   => "ssf_threshold_results.csv",
        :seed      => 42,
        :max_iters => 100,
        :t         => 2,
    )
    i = 1
    while i <= length(args)
        if args[i] == "--samples"
            defaults[:samples] = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--outfile"
            defaults[:outfile] = args[i+1]; i += 2
        elseif args[i] == "--seed"
            defaults[:seed] = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--max-iters"
            defaults[:max_iters] = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--t"
            defaults[:t] = parse(Int, args[i+1]); i += 2
        else
            @warn "Unknown argument: $(args[i])"
            i += 1
        end
    end
    return defaults
end

# ============================================================================
# Toric code construction (manual HGP from cyclic repetition code)
# ============================================================================

"""
    build_cyclic_repetition_code(L)

Build the (L × L) parity check matrix for the length-L cyclic repetition code.
H[i,i] = 1, H[i, (i mod L) + 1] = 1 for i = 1..L.
This is the classical code underlying the toric code HGP construction.
"""
function build_cyclic_repetition_code(L::Int)
    H = spzeros(Bool, L, L)
    for i in 1:L
        H[i, i] = true
        H[i, mod1(i + 1, L)] = true
    end
    return H
end

"""
    build_toric_hgp(L)

Build Hx and Hz for an [[2L², 2, L]] toric code via the hypergraph product
of two cyclic repetition codes of length L.

Returns (Hx, Hz, n, k) where n = 2L² is the number of qubits and k = 2.
"""
function build_toric_hgp(L::Int)
    h = build_cyclic_repetition_code(L)  # L × L sparse Bool
    I_L = sparse(Bool.(Matrix(I, L, L)))

    # HGP: Hx = [h ⊗ I_L | I_L ⊗ h^T]
    #       Hz = [I_L ⊗ h  | h^T ⊗ I_L]
    Hx = hcat(kron(h, I_L), kron(I_L, sparse(h')))
    Hz = hcat(kron(I_L, h), kron(sparse(h'), I_L))

    # convert to Bool sparse
    Hx = SparseMatrixCSC{Bool, Int}(Hx)
    Hz = SparseMatrixCSC{Bool, Int}(Hz)

    n = size(Hx, 2)
    # rank() needs numeric matrix, not Bool
    k = n - rank(Matrix{Float64}(Hx)) - rank(Matrix{Float64}(Hz))  # should be 2 for toric code

    return Hx, Hz, n, k
end

# ============================================================================
# Logical operator construction for toric code
# ============================================================================

"""
    build_toric_logicals(L)

Return (Lx, Lz) logical operator matrices for the [[2L², 2, L]] toric code.
Each is a (2 × 2L²) sparse Bool matrix. The two X logicals are horizontal and
vertical loops; the two Z logicals are the dual loops.

A logical error is detected when (error ⊕ correction) has non-trivial overlap
with the logical operators (i.e. the residual anticommutes with at least one
logical operator from the complementary type).
"""
function build_toric_logicals(L::Int)
    n = 2 * L^2

    # The first L² qubits correspond to horizontal edges (h ⊗ I block).
    # The second L² qubits correspond to vertical edges (I ⊗ h block).
    # Index a horizontal edge as (row, col) → (row-1)*L + col, 1-indexed.
    # Index a vertical edge as L² + (row-1)*L + col.

    # X logical 1: horizontal loop across row 1 (all horizontal edges in row 1)
    lx1 = spzeros(Bool, 1, n)
    for col in 1:L
        lx1[1, col] = true  # horizontal edge (1, col)
    end

    # X logical 2: vertical loop down column 1 (all vertical edges in col 1)
    lx2 = spzeros(Bool, 1, n)
    for row in 1:L
        lx2[1, L^2 + (row - 1) * L + 1] = true  # vertical edge (row, 1)
    end

    # Z logical 1: vertical loop down column 1 (vertical edges)
    lz1 = spzeros(Bool, 1, n)
    for row in 1:L
        lz1[1, L^2 + (row - 1) * L + 1] = true
    end

    # Z logical 2: horizontal loop across row 1 (horizontal edges)
    lz2 = spzeros(Bool, 1, n)
    for col in 1:L
        lz2[1, col] = true
    end

    Lx = vcat(lx1, lx2)
    Lz = vcat(lz1, lz2)

    return Lx, Lz
end

# ============================================================================
# Logical error check
# ============================================================================

"""
    is_logical_error(residual_x, residual_z, Lx, Lz)

Check if the residual error (error ⊕ correction) is a non-trivial logical operator.

For CSS codes, a logical X error is detected when the X-residual anticommutes
with at least one Z logical (i.e., Lz · residual_x ≠ 0 mod 2), and similarly
for Z errors.

Returns true if a logical error occurred.
"""
function is_logical_error(residual_x::AbstractVector, residual_z::AbstractVector,
                          Lx::AbstractMatrix, Lz::AbstractMatrix)
    # X-type logical error: Lz · residual_x mod 2 has any nonzero entry
    x_logical = any(Bool.(mod.(Lz * residual_x, 2)) .!= 0)
    # Z-type logical error: Lx · residual_z mod 2 has any nonzero entry
    z_logical = any(Bool.(mod.(Lx * residual_z, 2)) .!= 0)
    return x_logical || z_logical
end

# ============================================================================
# Monte Carlo simulation
# ============================================================================

"""
    run_threshold_sweep(distances, p_range, samples, max_iters, t, seed)

Run Monte Carlo threshold estimation for SSF decoder on toric codes.

For each code distance L and physical error rate p:
1. Generate i.i.d. X and Z errors with probability p.
2. Compute syndromes.
3. Decode with CSSSmallSetFlipDecoder.
4. Check for logical errors (not just syndrome matching).

Returns a vector of named tuples (distance, p, logical_error_rate, decode_failure_rate, samples).
"""
function run_threshold_sweep(distances::Vector{Int}, p_range, samples::Int,
                             max_iters::Int, t::Int, seed::Int)
    results = NamedTuple{(:distance, :p, :logical_error_rate, :decode_failure_rate, :samples),
                         Tuple{Int, Float64, Float64, Float64, Int}}[]

    for L in distances
        @printf("Building toric code L=%d ...\n", L)
        Hx, Hz, n, k = build_toric_hgp(L)
        Lx, Lz = build_toric_logicals(L)

        # Verify CSS commutativity: Hx * Hz' = 0 mod 2
        comm_check = mod.(Hx * Hz', 2)
        @assert nnz(comm_check) == 0 "CSS commutativity violated for L=$L"

        decoder = CSSSmallSetFlipDecoder(Hx, Hz, max_iters; t=t)

        @printf("  n=%d, k=%d, stabilizer checks: %d (X) + %d (Z)\n",
                n, k, size(Hx, 1), size(Hz, 1))

        for p in p_range
            rng = MersenneTwister(seed)
            logical_errors = 0
            decode_failures = 0

            for trial in 1:samples
                # Generate i.i.d. X and Z errors independently
                error_x = rand(rng, n) .< p
                error_z = rand(rng, n) .< p

                syn_x = Bool.(mod.(Hx * error_z, 2))  # X-checks detect Z-errors
                syn_z = Bool.(mod.(Hz * error_x, 2))  # Z-checks detect X-errors

                LDPCDecoders.reset!(decoder)
                (guess_x, guess_z), converged = decode!(decoder, syn_x, syn_z)

                if !converged
                    decode_failures += 1
                    logical_errors += 1  # count decode failure as logical error
                    continue
                end

                # Compute residual: error ⊕ correction
                residual_x = Int.(error_x) .⊻ Int.(guess_x)
                residual_z = Int.(error_z) .⊻ Int.(guess_z)

                if is_logical_error(residual_x, residual_z, Lx, Lz)
                    logical_errors += 1
                end
            end

            ler = logical_errors / samples
            dfr = decode_failures / samples
            push!(results, (distance=L, p=p, logical_error_rate=ler,
                            decode_failure_rate=dfr, samples=samples))

            @printf("  L=%2d  p=%.4f  logical_err=%.4f  decode_fail=%.4f\n",
                    L, p, ler, dfr)
        end
    end

    return results
end

# ============================================================================
# Output
# ============================================================================

function write_results(results, outfile)
    header = ["distance" "p" "logical_error_rate" "decode_failure_rate" "samples"]

    data = hcat(
        [r.distance for r in results],
        [r.p for r in results],
        [r.logical_error_rate for r in results],
        [r.decode_failure_rate for r in results],
        [r.samples for r in results],
    )

    open(outfile, "w") do io
        writedlm(io, vcat(header, data), ',')
    end
    @printf("Results written to %s\n", outfile)
end

function print_summary(results, distances)
    println("\n" * "="^60)
    println("THRESHOLD ESTIMATION SUMMARY")
    println("="^60)
    println()
    println("Look for the crossing point of logical error rate curves")
    println("across different code distances — that's the threshold p_th.")
    println()

    # For each pair of adjacent distances, find approximate crossing
    for i in 1:length(distances)-1
        L1, L2 = distances[i], distances[i+1]
        r1 = filter(r -> r.distance == L1, results)
        r2 = filter(r -> r.distance == L2, results)

        # Find where the curves cross (L2 goes from better to worse relative to L1)
        crossings = Float64[]
        for j in 1:min(length(r1), length(r2))-1
            diff_j   = r2[j].logical_error_rate   - r1[j].logical_error_rate
            diff_j1  = r2[j+1].logical_error_rate - r1[j+1].logical_error_rate
            if diff_j * diff_j1 < 0  # sign change → crossing
                # Linear interpolation
                p1, p2 = r1[j].p, r1[j+1].p
                frac = abs(diff_j) / (abs(diff_j) + abs(diff_j1))
                p_cross = p1 + frac * (p2 - p1)
                push!(crossings, p_cross)
            end
        end

        if !isempty(crossings)
            @printf("  Crossing L=%d vs L=%d:  p_th ≈ %.4f\n", L1, L2, mean(crossings))
        else
            @printf("  Crossing L=%d vs L=%d:  not found (try wider p range or more samples)\n", L1, L2)
        end
    end
    println()
end

# simple mean as we don't want to pull in Statistics.jl
mean(x) = sum(x) / length(x)

# ============================================================================
# Main
# ============================================================================

function main()
    config = parse_args(ARGS)

    distances = [3, 5, 7, 9]
    p_range = range(0.02, 0.12; length=15)

    @printf("SSF Threshold Estimation\n")
    @printf("========================\n")
    @printf("Code family:   Toric codes (HGP from cyclic repetition codes)\n")
    @printf("Distances:     %s\n", string(distances))
    @printf("Error rates:   %.3f to %.3f (%d points)\n", first(p_range), last(p_range), length(p_range))
    @printf("Samples/point: %d\n", config[:samples])
    @printf("SSF t:         %d\n", config[:t])
    @printf("Max iters:     %d\n", config[:max_iters])
    @printf("RNG seed:      %d\n", config[:seed])
    @printf("Output:        %s\n\n", config[:outfile])

    results = run_threshold_sweep(distances, p_range, config[:samples],
                                  config[:max_iters], config[:t], config[:seed])

    write_results(results, config[:outfile])
    print_summary(results, distances)
end

main()
