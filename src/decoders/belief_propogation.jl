using SparseArrays

include("abstract_decoder.jl")

struct BeliefPropagationSetup
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


struct BeliefPropagationDecoder{T} <: AbstractDecoder
  "Parity check matrix"
  H::T

  "Physical error rate"
  per::Float64

  "Number of max iterations of Belief propagation decoder"
  max_iters::Int
end

struct BeliefPropagationDecoderSparse{T} <: AbstractDecoder
  "Parity check matrix"
  H::T

  "Physical error rate"
  per::Float64

  "Number of max iterations of Belief propagation decoder"
  max_iters::Int
  
  "Num of parity checks i.e. number of rows of parity check matrix"
  s::Int
  
  "Num of bits in the code i.e number of columns of parity check matrix"
  n::Int
  
  sparse_H::SparseArrays.SparseMatrixCSC{Bool,Int}
  
  sparse_HT::SparseArrays.SparseMatrixCSC{Bool,Int}
end

"""
Function to reset the scratch space for Belief propagation decoding

# Arguments
* `bp_setup`: The Belief Propagation Scratch space
* `bp_decoder`: The belief propagation decoder configuration

# Examples
julia> decode!(bpd, syndrome, errors)
"""
function reset!(bp_setup::BeliefPropagationSetup, bp_decoder::BeliefPropagationDecoder)
  bp_setup.log_probabs .= 0.0
  bp_setup.channel_probs .= bp_decoder.per
  bp_setup.bit_2_check .= 0.0
  bp_setup.check_2_bit .= 0.0
  bp_setup.err .= 0.0
end

# TODO: Belief decoder type, syndrome, error -> should be modified 

"""
Function to decode given the parity check matrix, syndrome and error

# Arguments
* `decoder`: The Belief Propagation Decoder configuration
* `syndrome`: Syndrome that is taken as input
* `errors`: Predefined error array that this function manipulates

# Examples
```julia
julia> decode!(bpd, syndrome, error)
```
"""
function decode!(decoder::BeliefPropagationDecoder, syndrome::BitArray{1}, error::BitArray{1})
  s, n = size(decoder.H)
  
  sparse_decoder = BeliefPropagationDecoderSparse{typeof(decoder.H)}(decoder.H, decoder.per, decoder.max_iters, s, n, sparse(decoder.H), sparse(decoder.H'))
  setup = BeliefPropagationSetup(zeros(n), fill(decoder.per, n), zeros(s,n), zeros(s,n), zeros(n))

  reset!(setup, decoder)
  success = syndrome_decode!(sparse_decoder, setup, syndrome)
  
  for i in eachindex(setup.err)
    error[i] = setup.err[i]
  end

  return success
end


"""
Function to batch decode given the parity check matrix, syndrome and error

Scratch space allocations are done once and re-used for better performance

# Arguments
* `decoder`: The Belief Propagation decoder configuration
* `syndromes`: Syndrome matrix that contains batch syndromes
* `errors`: Predefined error matrix that this function manipulates

# Examples
```julia
julia> batchdecode!(bpd, syndromes, errors)
```
"""
function batchdecode!(decoder::BeliefPropagationDecoder, syndromes::BitArray{2}, errors::BitArray{2})

  @assert size(syndromes, 1) == size(errors, 1)
  num_trials::Int = size(syndromes, 1)

  s, n = size(decoder.H)
  sparse_decoder = BeliefPropagationDecoderSparse{typeof(decoder.H)}(decoder.H, decoder.per, decoder.max_iters, s, n, sparse(decoder.H), sparse(decoder.H'))
  setup = BeliefPropagationSetup(zeros(n), fill(decoder.per, n), zeros(s,n), zeros(s,n), zeros(n))

  success::BitArray{1} = zeros(num_trials)

  for i in axes(syndromes, 1)
    reset!(setup, decoder)
    success[i] = syndrome_decode!(sparse_decoder, setup, syndromes[i, :])

    errors[i, :] = setup.err
  end

  return success
end
