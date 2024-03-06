using SparseArrays

include("abstract_decoder.jl")

# TODO: Better nomenclature
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
end

function BeliefPropagationDecoder(H, per::Float64, max_iters::Int)
  s, n = size(H)
  sparse_H = sparse(H)
  sparse_HT = sparse(H')
  return BeliefPropagationDecoder(per, max_iters, s, n, sparse_H, sparse_HT)
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
function reset!(bp_setup::BeliefPropagationScratchSpace, bp_decoder::BeliefPropagationDecoder)
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
```jldoctest
julia> decoder = BeliefPropagationDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> error = rand(1000) .< 0.01;

julia> syndrome::BitArray{1} = (H * error) .% 2;

julia> decoded_error::BitArray{1} = zeros(1000);

julia> decode!(decoder, syndrome, decoded_error)
true
```
"""

function decode!(decoder::BeliefPropagationDecoder, syndrome::BitArray{1}, error::BitArray{1})

  setup = BeliefPropagationScratchSpace(zeros(decoder.n), fill(decoder.per, decoder.n), 
            zeros(decoder.s, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.n))

  reset!(setup, decoder)
  success = syndrome_decode!(decoder, setup, syndrome)
  copyto!(error, setup.err)
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
```jldoctest
julia> decoder = BeliefPropagationDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> errors::BitArray{2} = rand(Base.Float64, (10,1000)) .< 0.01;

julia> syndromes::BitArray{2} = zeros(10, 900);

julia> for i in axes(errors, 1)
             syn = (H * errors[i, :]) .% 2;
             syndromes[i, :] = syn;
           end  

julia> batchdecode!(decoder, syndromes, errors) < 0.9
true
"""
function batchdecode!(decoder::BeliefPropagationDecoder, syndromes::BitArray{2}, errors::BitArray{2})

  @assert size(syndromes, 1) == size(errors, 1)
  num_trials::Int = size(syndromes, 1)

  setup = BeliefPropagationScratchSpace(zeros(decoder.n), fill(decoder.per, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.n))

  success::BitArray{1} = zeros(num_trials)

  for i in axes(syndromes, 1)
    reset!(setup, decoder)
    success[i] = syndrome_decode!(decoder, setup, syndromes[i, :])
    errors[i, :] = setup.err
  end

  return success
end