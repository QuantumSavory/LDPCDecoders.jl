struct CSSSmallSetFlipDecoder <: AbstractDecoder
    x_decoder::SmallSetFlipDecoder
    z_decoder::SmallSetFlipDecoder
end

function CSSSmallSetFlipDecoder(Hx::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, Hz::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, max_iters::Int; t::Int=0, rng::Union{Nothing, AbstractRNG}=nothing)
    x_decoder = SmallSetFlipDecoder(Hz, max_iters; t=t, rng=rng)
    z_decoder = SmallSetFlipDecoder(Hx, max_iters; t=t, rng=rng)
    return CSSSmallSetFlipDecoder(x_decoder, z_decoder)
end

function reset!(decoder::CSSSmallSetFlipDecoder)
    reset!(decoder.x_decoder)
    reset!(decoder.z_decoder)
    return decoder
end

function decode!(decoder::CSSSmallSetFlipDecoder, syn_x::AbstractVector, syn_z::AbstractVector)
    err_x, conv_x = decode!(decoder.x_decoder, syn_z)
    err_z, conv_z = decode!(decoder.z_decoder, syn_x)
    return (err_x=err_x, err_z=err_z), conv_x && conv_z
end

function decode!(decoder::CSSSmallSetFlipDecoder, syndrome::AbstractVector)
    sx = decoder.z_decoder.s
    sz = decoder.x_decoder.s
    @assert length(syndrome) == sx + sz "Syndrome length mismatch: expected $(sx + sz), got $(length(syndrome))"
    
    syn_x = @view syndrome[1:sx]
    syn_z = @view syndrome[sx+1:end]
    
    (err_x, err_z), conv = decode!(decoder, syn_x, syn_z)
    return vcat(err_x, err_z), conv
end
