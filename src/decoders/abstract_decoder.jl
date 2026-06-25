"""
    AbstractDecoder

Abstract supertype for all LDPC decoders in this package.

Concrete decoders (e.g., [`BeliefPropagationDecoder`](@ref), [`BitFlipDecoder`](@ref))
subtype `AbstractDecoder` and must implement [`decode!`](@ref).
"""
abstract type AbstractDecoder end

"""
    batchdecode!(decoder::AbstractDecoder, syndromes::AbstractMatrix, errors::AbstractMatrix, [converged::AbstractVector{Bool}])

Decode each column of `syndromes` with `decoder`, writing the results into the
corresponding column of `errors`.

This is the generic fallback for any `AbstractDecoder` that implements
[`decode!`](@ref). Individual decoders can override it with a faster batch
strategy.

# Arguments
- `decoder`: Any concrete `AbstractDecoder`.
- `syndromes`: Matrix whose columns are the syndromes to decode.
- `errors`: Pre-allocated matrix (same column count); overwritten with results.
- `converged`: Optional pre-allocated boolean vector. When provided, the function runs with zero allocations.

Returns `(errors, converged)`. The `converged[i]` is `true` if the decoder found
an error estimate whose syndrome matches `syndromes[:, i]` within the iteration
limit, and `false` otherwise.
"""
function batchdecode!(decoder::AbstractDecoder, syndromes::AbstractMatrix, errors::AbstractMatrix, converged::AbstractVector{Bool})
    @assert size(syndromes, 2) == size(errors, 2)
    @assert size(syndromes, 2) == length(converged)

    @views for i in axes(syndromes, 2)
        guess, conv = decode!(decoder, syndromes[:, i])
        converged[i] = conv
        errors[:, i] .= guess
    end

    return errors, converged
end

function batchdecode!(decoder::AbstractDecoder, syndromes::AbstractMatrix, errors::AbstractMatrix)
    num_trials = size(syndromes, 2)
    converged = Vector{Bool}(undef, num_trials)
    return batchdecode!(decoder, syndromes, errors, converged)
end
