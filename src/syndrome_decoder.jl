function syndrome_decode(pcm, syndrome, error_rate, max_iters)
  
  # Get size of Parity check matrix
  m, n = size(pcm)

  # Initialise channel probabilities
  channel_probs = fill(error_rate, n)
  
  # Initialise messages between checks, bits array
  # 3rd dimension: 1 -> bit to check, 2-> check to bit
  msg = zeros(m, n, 2)

  # Initialize log probabilities
  log_probabs = zeros(n)

  # Initiliase bit to check messages
  for j in 1:n
    for i in 1:m
      if pcm[i,j] == 1
        msg[i, j, 1] = (channel_probs[j]) / (1 - channel_probs[j])
      end
    end
  end

  error = zeros(Bool,n)
  converged = 0
  for iter in 1:max_iters

    # Check to bit messages
    for i in 1:m
      temp = (-1) ^ syndrome[i]
      for j in 1:n
        if pcm[i,j] == 1
          msg[i,j,2] = temp
          temp *= 2 / (1 + msg[i,j,1]) - 1
        end
      end

      temp = 1.0
      for j=n:-1:1
        if pcm[i,j] == 1
          msg[i,j,2] *= temp
          msg[i,j,2] = (1 - msg[i,j,2]) / (1 + msg[i,j,2])
          temp *= 2/(1 + msg[i,j,1]) - 1
        end
      end
    end

    # display(msg)
    # Bit to check messages
    for j in 1:n
      temp = channel_probs[j] / (1 - channel_probs[j])

      for i in 1:m
        if pcm[i,j] == 1
          msg[i,j,1] = temp
          temp *= msg[i,j,2]
          if isnan(temp)
            temp = 1.0
          end
        end
      end

      log_probabs[j] = log(1 / temp)
      if temp >= 1
        error[j] = 1
      else 
        error[j] = 0
      end

      temp = 1.0
      for i=m:-1:1
        if pcm[i,j] == 1
          msg[i,j,1] *= temp
          temp *= msg[i,j,2]
          if isnan(temp)
            temp = 1.0
          end
        end
      end
    end

    syndrome_decoded = (pcm * error) .% 2
    if all(syndrome_decoded .== syndrome)
      converged = true
      return error, converged
    end
  end
  
  return Bool.(error), converged

end