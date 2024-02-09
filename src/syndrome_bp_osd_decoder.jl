using LinearAlgebra

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