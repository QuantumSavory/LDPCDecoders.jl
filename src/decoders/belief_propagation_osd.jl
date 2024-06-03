using LinearAlgebra

using SparseArrays

struct BPOSDScratchSpace
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

function BPOSDScratchSpace(n, s, per)
  return BPOSDScratchSpace(zeros(n), fill(per, n), zeros(s, n), zeros(s, n), zeros(n))
end

struct BPOSDDecoder <: AbstractDecoder
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
  scratch::BPOSDScratchSpace
end

function BPOSDDecoder(H, per::Float64, max_iters::Int)
  s, n = size(H)
  sparse_H = sparse(H)
  sparse_HT = sparse(H')
  scratch = BPOSDScratchSpace(n, s, per)
  return BPOSDDecoder(per, max_iters, s, n, sparse_H, sparse_HT, scratch)
end

"""
Function to reset the scratch space for Belief propagation decoding

# Arguments
* `bp_setup`: The Belief Propagation Scratch space
* `bp_osd_decoder`: The belief propagation decoder configuration

# Examples
```jldoctest
julia> decoder = BPOSDDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> scratch = BPOSDScratchSpace(zeros(decoder.n), fill(decoder.per, decoder.n),
                zeros(decoder.s, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.n));

julia> reset!(scratch, decoder);
````
"""
function reset!(bp_osd_decoder::BPOSDDecoder)
  bp_setup = bp_osd_decoder.scratch
  bp_setup.log_probabs .= 0.0
  bp_setup.channel_probs .= bp_osd_decoder.per
  bp_setup.bit_2_check .= 0.0
  bp_setup.check_2_bit .= 0.0
  bp_setup.err .= 0.0
  bp_osd_decoder
end

# TODO: Belief decoder type, syndrome, error -> should be modified

"""
Function to decode given the parity check matrix, syndrome and error

# Arguments
* `decoder`: The Belief Propagation OSD Decoder configuration
* `syndrome`: Syndrome that is taken as input
* `errors`: Predefined error array that this function manipulates

# Examples
```jldoctest
julia> H = LDPCDecoders.parity_check_matrix(1000, 10, 9);

julia> decoder = BPOSDDecoder(H, 0.01, 100);

julia> error = rand(1000) .< 0.01;

julia> syndrome = (H * error) .% 2;

julia> guess, success = decode!(decoder, syndrome);

julia> error == guess
true
```
"""
function decode!(decoder::BPOSDDecoder, syndrome::AbstractVector) # TODO check if casting to bitarrays helps with performance -- if it does, set up warnings to the user for cases where they have not done the casting
  bp_scratch = BeliefPropagationScratchSpace(decoder.n, decoder.s, decoder.per)
  bp_decoder = BeliefPropagationDecoder(
      decoder.per, decoder.max_iters, decoder.s, 
      decoder.n, decoder.sparse_H, decoder.sparse_HT, 
      bp_scratch
    )
  
  setup = decoder.scratch
  guess_error, converged = decode!(bp_decoder, syndrome)
  
  if !converged  
    dense_H = Nemo.matrix(Nemo.GF(2), collect(decoder.sparse_H))
    rank_H = LinearAlgebra.rank(dense_H)
    sorted_cols = sortperm(setup.log_probabs, alg=QuickSort, rev=true)
    trim_cols = sorted_cols[1:rank_H]
    trim_H = dense_H[:, trim_cols]
    syndrome_gf2 = Nemo.matrix(Nemo.GF(2), decoder.s,1,syndrome)
    _, trim_error = Nemo.can_solve_with_solution(trim_H, syndrome_gf2, side=:right)
    trim_error_vec = trim_error .== 1
    
    setup.err .= 0
    setup.err[trim_cols] .= trim_error_vec
    
    syndrome_decoded = (decoder.sparse_H * setup.err) .% 2 
    converged = all(syndrome_decoded .== syndrome)
  else
    setup.err = guess_error
  end

  return bp_decoder.scratch.err, converged
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
julia> decoder = BPOSDDecoder(LDPCDecoders.parity_check_matrix(1000, 10, 9), 0.01, 100);

julia> samples = 100

julia> errors = rand(1000,samples) .< 0.01;

julia> syndromes = zeros(900, samples);

julia> syndromes = H*errors .% 2

julia> guesses, successes = batchdecode!(decoder, syndromes, zero(errors));
`
julia> sum((guesses[:,i] == errors[:,i] for i in 1:samples)) > 0.995*samples
"""
function batchdecode!(decoder::BPOSDDecoder, syndromes, errors)
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

function bp_osd0_decode(pcm, pcmT, syndrome, maxIters, channelProbs, b2c, c2b, logProbabs, error)
  pcmRank = rank(pcm)
  # println(typeof(pcm))
  pcm = sparse(pcm)
  pcmT = sparse(pcmT)  
  decodedError, converged, probDecisions = syndrome_decode(pcm, pcmT, syndrome, maxIters, channelProbs, b2c, c2b, logProbabs, error)
  
  # If converged, return decoded_error with belief propogation
  if converged 
    return decodedError, converged
  end

  
  # println("Rank of the pcm matrix is : $pcmRank")
  numChecks, numBits = size(pcm)

  sortedCols = sortperm(probDecisions)
  sortedColsTruncated = sortedCols[1:pcmRank]
  # println("The size of sortedColsTruncated : ")
  # display(size(sortedColsTruncated))
  
  pcmTruncated = pcm[:,sortedColsTruncated]
  # println("The size of pcmTruncated : ")
  # display(size(pcmTruncated))
  
  # println(typeof(pcmTruncated))
  # println(typeof(syndrome))
  errorS = BitArray(pcmTruncated) \ syndrome
  
  errorOsd0 = zeros(numBits)

  for idx in eachindex(errorS)
    if errorS[idx] == 1
      errorOsd0[sortedColsTruncated[idx]] = 1
    end
  end

  osd0Converged = false
  decodedSyndrome = (pcm * errorOsd0) .% 2
    if all(decodedSyndrome .== syndrome)
      osd0Converged = true
      return Bool.(error), osd0Converged
    end
  
  if converged != osd0Converged
    println("OSD helped !")
  end
  return errorOsd0, osd0Converged

end