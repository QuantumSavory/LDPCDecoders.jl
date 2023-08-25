function syndrome_it_decode(pcm, syndrome, max_iters::Int, error, votes)
  
  # Get size of parity check matrix 
  ## n is codeword length, m is number of parity check equations
  m, n = size(pcm)

  # Initialise error and success
  # error = zeros(Bool, n)
  success = false

  # Iterative decoding
  for iter in 1:max_iters
    curr_syn = (pcm * error) .% 2
    # println("Curr syn : $curr_syn, Syndrome : $syndrome")

    # Check if the syndrome matches
    if curr_syn == syndrome
      success = true
      break
    end

    # Initialise votes
    # votes = zeros(Int, n)
    # votes = vec(votes)

    error_checks = .!(curr_syn .== syndrome)
    
    for i in 1:m
      if error_checks[i] == 1
        votes += pcm[i, :]
      else
        votes -= pcm[i, :]
      end
    end

    max_idx = argmax(votes)
    error[max_idx] = 1 - error[max_idx]
  end

  return error, success
end