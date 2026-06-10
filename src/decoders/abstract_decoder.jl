"""
    AbstractDecoder

Abstract supertype for all LDPC decoders in this package.

Concrete decoders (e.g., [`BeliefPropagationDecoder`](@ref), [`BitFlipDecoder`](@ref))
subtype `AbstractDecoder` and must implement [`decode!`](@ref).
"""
abstract type AbstractDecoder end

"""
    batchdecode!(decoder::AbstractDecoder, syndromes, errors)

Decode every column of `syndromes` with `decoder`, writing the guesses into
`errors` column-by-column.

This generic fallback works for any `AbstractDecoder` that implements
[`decode!`](@ref).  Individual decoders may override it if a more efficient
batch strategy exists.

# Arguments
- `decoder`: Any concrete `AbstractDecoder`.
- `syndromes`: Matrix whose columns are the syndromes to decode.
- `errors`: Pre-allocated matrix (same column count); overwritten with guesses.

Returns `(errors, converged)` where `converged` is a `Vector{Bool}`.
"""
function batchdecode!(decoder::AbstractDecoder, syndromes, errors)
    @assert size(syndromes, 2) == size(errors, 2)
    num_trials = size(syndromes, 2)
    converged = Vector{Bool}(undef, num_trials)

    for i in axes(syndromes, 2)
        guess, conv = decode!(decoder, syndromes[:, i])
        converged[i] = conv
        errors[:, i] = guess
    end

    return errors, converged
end
