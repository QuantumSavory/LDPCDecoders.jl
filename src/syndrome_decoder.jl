using SparseArrays

function syndrome_decode(pcm, pcmT, syndrome, max_iters, channel_probs, b2c, c2b, log_probabs, error)
  
  # Get size of Parity check matrix
  m, n = size(pcm)
  rows = rowvals(pcm)
  rowsT = rowvals(pcmT)

  # Initiliase bit to check messages
  for j in 1:n
    for k in nzrange(pcm, j)
      i = rows[k]
      b2c[i, j] = (channel_probs[j]) / (1 - channel_probs[j])
    end
  end

  converged = false
  for iter in 1:max_iters

    # Check to bit messages
    for i in 1:m
      temp = (-1) ^ syndrome[i]
      # @inbounds for j in 1:n
      for k in nzrange(pcmT, i)
          j = rowsT[k]
          c2b[i,j] = temp
          temp *= 2 / (1 + b2c[i,j]) - 1
      end

      temp = 1.0
      # @inbounds for j=n:-1:1
      for k in reverse(nzrange(pcmT, i))
          j = rowsT[k]
          c2b[i,j] *= temp
          c2b[i,j] = (1 - c2b[i,j]) / (1 + c2b[i,j])
          temp *= 2/(1 + b2c[i,j]) - 1
      end
    end

    
    # Bit to check messages
    for j in 1:n
      temp = channel_probs[j] / (1 - channel_probs[j])

      for k in nzrange(pcm, j)
        i = rows[k]
          b2c[i,j] = temp
          temp *= c2b[i,j]
          if isnan(temp)
            temp = 1.0
          end
      end

      log_probabs[j] = log(1 / temp)
      if temp >= 1
        error[j] = 1
      else 
        error[j] = 0
      end

      temp = 1.0
      for k in reverse(nzrange(pcm, j))
          i = rows[k]
          b2c[i,j] *= temp
          temp *= c2b[i,j]
          if isnan(temp)
            temp = 1.0
          end
        # end
      end
    end

    syndrome_decoded = (pcm * error) .% 2
    if all(syndrome_decoded .== syndrome)
      converged = true
      return Bool.(error), converged, log_probabs
    end
  end
  
  return Bool.(error), converged, log_probabs

end