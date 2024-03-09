using SparseArrays

struct BeliefPropagationScratchSpace
  "Log probabilities"
  log_probabs::Vector{Float64}

  "Channel probabilities"
  channel_probs::Vector{Float64}

  "Bit to check message matrix"
  bit_2_check::Matrix{Float64}

  "Check to bit message matrix"
  check_2_bit::Matrix{Float64}

  "predicted error for each bit"
  err::Vector{Float64}
end

function BeliefPropagationScratchSpace(n, s, per)
  return BeliefPropagationScratchSpace(zeros(n), fill(per, n), zeros(s, n), zeros(s, n), zeros(n))
end

struct BeliefPropagationDecoder <: AbstractDecoder
  "Physical error rate"
  per::Float64

  "Number of max iterations of Belief propagation decoder"
  max_iters::Int

  "Num of parity checks i.e. number of rows of parity check matrix"
  s::Int

  "Num of bits in the code i.e number of columns of parity check matrix"
  n::Int

  "Sparse form of the parity check matrix"
  sparse_H::SparseArrays.SparseMatrixCSC{Bool,Int}

  "Sparse form of the transform of the parity check matrix"
  sparse_HT::SparseArrays.SparseMatrixCSC{Bool,Int}

  "Scratch space for the decoder"
  scratch::BeliefPropagationScratchSpace
end

function BeliefPropagationDecoder(H, per::Float64, max_iters::Int)
  s, n = size(H)
  sparse_H = sparse(H)
  sparse_HT = sparse(H')
  scratch = BeliefPropagationScratchSpace(n, s, per)
  return BeliefPropagationDecoder(per, max_iters, s, n, sparse_H, sparse_HT, scratch)
end

"""
Function to reset the scratch space for Belief propagation decoding

# Arguments
* `bp_setup`: The Belief Propagation Scratch space
* `bp_decoder`: The belief propagation decoder configuration

# Examples
```jldoctest
julia> decoder = BeliefPropagationDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> scratch = BeliefPropagationScratchSpace(zeros(decoder.n), fill(decoder.per, decoder.n),
                zeros(decoder.s, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.n));

julia> reset!(scratch, decoder);
````
"""
function reset!(bp_decoder::BeliefPropagationDecoder)
  bp_setup = bp_decoder.scratch
  bp_setup.log_probabs .= 0.0
  bp_setup.channel_probs .= bp_decoder.per
  bp_setup.bit_2_check .= 0.0
  bp_setup.check_2_bit .= 0.0
  bp_setup.err .= 0.0
  bp_decoder
end

# TODO: Belief decoder type, syndrome, error -> should be modified

"""
Function to decode given the parity check matrix, syndrome and error

# Arguments
* `decoder`: The Belief Propagation Decoder configuration
* `syndrome`: Syndrome that is taken as input
* `errors`: Predefined error array that this function manipulates

# Examples
```jldoctest
julia> H = LDPCDecoders.parity_check_matrix(1000, 10, 9);

julia> decoder = BeliefPropagationDecoder(H, 0.01, 100);

julia> error = rand(1000) .< 0.01;

julia> syndrome = (H * error) .% 2;

julia> guess, success = decode!(decoder, syndrome, zeros(1000));

julia> error == guess
true
```
"""
function decode!(decoder::BeliefPropagationDecoder, syndrome::AbstractVector, error::AbstractVector) # TODO check if casting to bitarrays helps with performance -- if it does, set up warnings to the user for cases where they have not done the casting
  reset!(decoder)
  success = syndrome_decode!(decoder, decoder.scratch, syndrome) # TODO: Delete syndrome_decode! and just move its content here. There is no point in this indirection.
  return decoder.scratch.err, success
end

"""
Function to batch decode given the parity check matrix, syndrome and error

Scratch space allocations are done once and re-used for better performance

# Arguments
* `decoder`: The Belief Propagation decoder configuration
* `syndromes`: Syndrome matrix that contains batch syndromes
* `errors`: Predefined error matrix that this function manipulates

# Examples
```jldoctest
julia> decoder = BeliefPropagationDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> samples = 100

julia> errors = rand(1000,samples) .< 0.01;

julia> syndromes = zeros(900, samples);

julia> syndromes = H*errors .% 2

julia> guesses, successes = batchdecode!(decoder, syndromes, zero(errors));

julia> sum((guesses[:,i] == errors[:,i] for i in 1:samples)) > 0.995*samples
"""
function batchdecode!(decoder::BeliefPropagationDecoder, syndromes, errors)
  @assert size(syndromes, 2) == size(errors, 2)
  num_trials::Int = size(syndromes, 2)
  success::AbstractVector{Bool} = zeros(num_trials)

  for i in axes(syndromes, 2)
    reset!(decoder)
    success[i] = syndrome_decode!(decoder, decoder.scratch, syndromes[:, i])
    errors[:, i] = decoder.scratch.err
  end

  return errors, success
end
