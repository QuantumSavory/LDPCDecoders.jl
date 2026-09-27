"""
    CSSSmallSetFlipDecoder <: AbstractDecoder

CSS-aware Small-Set Flip decoder that decodes X-errors and Z-errors independently.

Holds two [`SmallSetFlipDecoder`](@ref) instances: `x_decoder` (using `Hz`) for X-errors
and `z_decoder` (using `Hx`) for Z-errors, following the CSS convention where Z-checks
detect X-errors and vice versa.

# Example
```jldoctest
julia> using LDPCDecoders, SparseArrays, StableRNGs

julia> H = sparse(Bool[1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1]);

julia> dec = CSSSmallSetFlipDecoder(H, H, 50; t=2);

julia> syn_x = zeros(Bool, 4); syn_z = Bool.(mod.(H * [1,0,0,0], 2));

julia> (gx, gz), conv = decode!(dec, syn_x, syn_z);

julia> conv
true
```
"""
struct CSSSmallSetFlipDecoder <: AbstractDecoder
    x_decoder::SmallSetFlipDecoder
    z_decoder::SmallSetFlipDecoder
end

"""
    CSSSmallSetFlipDecoder(Hx, Hz, max_iters; t=0, rng=nothing)

Construct a CSS SSF decoder from X and Z parity check matrices.
`Hx` is used for the `z_decoder` (detects Z-errors) and `Hz` for the `x_decoder` (detects X-errors).
"""
function CSSSmallSetFlipDecoder(Hx::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, Hz::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, max_iters::Int; t::Int=0, rng::Union{Nothing, AbstractRNG}=nothing)
    x_decoder = SmallSetFlipDecoder(Hz, max_iters; t=t, rng=rng)
    z_decoder = SmallSetFlipDecoder(Hx, max_iters; t=t, rng=rng)
    return CSSSmallSetFlipDecoder(x_decoder, z_decoder)
end

"""
    reset!(decoder::CSSSmallSetFlipDecoder)

Reset both internal sub-decoders.
"""
function reset!(decoder::CSSSmallSetFlipDecoder)
    reset!(decoder.x_decoder)
    reset!(decoder.z_decoder)
    return decoder
end

"""
    decode!(decoder::CSSSmallSetFlipDecoder, syn_x, syn_z)

Decode X-errors and Z-errors independently. Returns `((err_x, err_z), converged)`.
"""
function decode!(decoder::CSSSmallSetFlipDecoder, syn_x::AbstractVector, syn_z::AbstractVector)
    err_x, conv_x = decode!(decoder.x_decoder, syn_z)
    err_z, conv_z = decode!(decoder.z_decoder, syn_x)
    return (err_x=err_x, err_z=err_z), conv_x && conv_z
end

"""
    decode!(decoder::CSSSmallSetFlipDecoder, syndrome)

Decode from a combined syndrome vector `[σ_x; σ_z]`. Returns `([ê_x; ê_z], converged)`.
"""
function decode!(decoder::CSSSmallSetFlipDecoder, syndrome::AbstractVector)
    sx = decoder.z_decoder.s
    sz = decoder.x_decoder.s
    @assert length(syndrome) == sx + sz "Syndrome length mismatch: expected $(sx + sz), got $(length(syndrome))"
    
    syn_x = @view syndrome[1:sx]
    syn_z = @view syndrome[sx+1:end]
    
    (err_x, err_z), conv = decode!(decoder, syn_x, syn_z)
    return vcat(err_x, err_z), conv
end

