using SparseArrays

struct BitFlipScratchSpace
  "Predicted error for each bit"
  err::Vector{Int64}

  "Votes in favor to flip the bit"
  votes::Vector{Int64}

  "Guess syndrome"
  syn::Vector{Int64}

  "Error checks"
  error_checks::Vector{Int64}
end

function BitFlipScratchSpace(s, n)
  return BitFlipScratchSpace(zeros(n), zeros(n), zeros(s), zeros(s))
end

struct BitFlipDecoder <: AbstractDecoder
  "Physical error rate"
  per::Float64

  "Number of max iterations of Iterative bit flip decoder"
  max_iters::Int

  "Num of parity checks i.e. number of rows of parity check matrix"
  s::Int

  "Num of bits in the code i.e number of columns of parity check matrix"
  n::Int

  "Sparse form of the parity check matrix"
  sparse_H::SparseArrays.SparseMatrixCSC{Bool,Int}

  "Sparse form of the transform of the parity check matrix"
  sparse_HT::SparseArrays.SparseMatrixCSC{Bool,Int}

  "Non-zero rowvals of the transform of sparse parity check matrix"
  rowsT::Vector{Int}

  "Scratch space for the decoder"
  scratch::BitFlipScratchSpace
end

function BitFlipDecoder(H, per::Float64, max_iters::Int)
  s, n = size(H)
  scratch = BitFlipScratchSpace(s, n)
  sparse_H = sparse(H)
  sparse_HT = sparse(H')
  rowsT = rowvals(sparse_HT)
  return BitFlipDecoder(per, max_iters, s, n, sparse_H, sparse_HT, rowsT, scratch)
end


"""
Function to reset the scratch space for bit flip decoding algorithm

# Arguments
* `bf_decoder`: The Bit Flip decoder configuration

# Examples
```jldoctest
julia> decoder = BitFlipDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> reset!(decoder);
````
"""
function reset!(bf_decoder::BitFlipDecoder)
  scratch = bf_decoder.scratch
  scratch.err .= 0.0
  scratch.votes .= 0.0
end

"""
Function to decode the given syndrome 

# Arguments
* `decoder`: The Bit Flip Decoder configuration
* `syndrome`: Syndrome that is taken as input
* `errors`: Predefined error array that this function manipulates

# Examples
```jldoctest
julia> H = LDPCDecoders.parity_check_matrix(1000, 10, 9);

julia> decoder = BitFlipDecoder(H, 0.01, 100);

julia> error = rand(1000) .< 0.01;

julia> syndrome = (H * error) .% 2;

julia> guess, success = decode!(decoder, syndrome);

julia> error == guess
true
```
"""
function decode!(decoder::BitFlipDecoder, syndrome::AbstractVector)
  reset!(decoder)
  setup = decoder.scratch
  converged = false

  for iter in 1:decoder.max_iters
    setup.syn .= (decoder.sparse_H * setup.err) .% 2

    if setup.syn == syndrome
      converged = true
      break
    end

    setup.error_checks .= (setup.syn .!= syndrome)
    
    for i in 1:decoder.s
      if setup.error_checks[i] == 1
        for k::Int in nzrange(decoder.sparse_HT, i)
          j = decoder.rowsT[k]
          setup.votes[j] += decoder.sparse_H[i, j]
        end
      else
        for k::Int in nzrange(decoder.sparse_HT, i)
          j = decoder.rowsT[k]
          setup.votes[j] -= decoder.sparse_H[i, j]
        end
      end
    end

    max_votes = maximum(setup.votes)
    if max_votes >= 0
      max_idxs = findall(setup.votes .== max_votes)
      max_idx = rand(max_idxs)
      setup.err[max_idx] = 1 - setup.err[max_idx]
    else
      converged = true
      break 
    end
  end

  return setup.err, converged
end

"""
Function to batch decode given the parity check matrix, syndrome and errors

Scratch space allocations are done once and re-used for better performance

# Arguments
* `decoder`: The Bit flip decoder configuration
* `syndromes`: Syndrome matrix that contains batch syndromes
* `errors`: Predefined error matrix that this function manipulates

# Examples
```jldoctest
julia> decoder = BitFlipDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> samples = 100

julia> errors = rand(1000,samples) .< 0.01;

julia> syndromes = zeros(900, samples);

julia> syndromes = H*errors .% 2

julia> guesses, successes = batchdecode!(decoder, syndromes, zero(errors));

julia> sum((guesses[:,i] == errors[:,i] for i in 1:samples)) > 0.995*samples
"""
function batchdecode!(decoder::BitFlipDecoder, syndromes, errors)
  @assert size(syndromes, 2) == size(errors, 2)
  num_trials::Int = size(syndromes, 2)
  converged::AbstractVector{Bool} = zeros(num_trials)

  for i in axes(syndromes, 2)
    reset!(decoder)
    guess, conv = decode!(decoder, syndromes[:, i])
    # guess, conv = syndrome_it_decode(decoder.sparse_H, syndromes[:, i], decoder.max_iters, decoder.scratch.err, decoder.scratch.votes)
    converged[i] = conv
    errors[:, i] = guess
  end

  return errors, converged
end
