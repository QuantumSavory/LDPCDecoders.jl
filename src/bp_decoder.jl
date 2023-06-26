using LinearAlgebra


# Bi-Symmetric Channel LLR (Bit flip)
function llr(error_rate, bit)
  return ((-1) ^ bit) * log((1-error_rate)/error_rate)
end

# Calculate Phi
function phi(message)
  return -log(tanh(message * 0.5))
end

# Calculate estimate
function estimate(message_c2v, num_vs, num_cs, vllr, parity_check_matrix)
  estimate = zeros(num_vs)
  for v in 1:num_vs
    estimate[v] += vllr[v]
    for c in 1:num_cs
      if parity_check_matrix[c, v] == 1
        estimate[v] += message_c2v[c, v]
      end
    end
  end 
  return estimate
end

# Belief Propagation Decoder
function bp_decode(received_message, parity_check_matrix, error_rate, max_iterations=100)
  
  num_checks, num_bits = size(parity_check_matrix)
  num_edges = sum(parity_check_matrix)

  num_cs, num_vs = size(parity_check_matrix)
  
  # Initialize messages
  message_c2v = zeros(num_checks, num_bits)
  message_v2c = zeros(num_checks, num_bits)


  # Messages received by the variables 
  
  ## For every variable node, there should be initialised llr
  vllr = zeros(num_bits)
  for i in 1:num_bits
    vllr[i] = llr(error_rate, received_message[i])
  end 
  

  # Initial message received by check nodes
  for check in 1:num_cs
    for v in 1:num_vs
      if parity_check_matrix[check, v] == 1
        sum = 0
      
        for c in 1:num_cs
          if c != check && parity_check_matrix[c, v] == 1
            sum += message_c2v[c, v]  
          end
        end

        message_v2c[check, v] = sum + vllr[v]  
      end
    end  
  end

  local decoded
  local syndrome

  for iter in 1:max_iterations
    
    # Message receive by variable nodes
    for variable in 1:num_vs 
      for c in 1:num_cs
        if parity_check_matrix[c, variable] == 1 

          sum = 0
          sgn = 1
          for v in 1:num_vs 
            if parity_check_matrix[c, v] == 1 && v != variable 
              sgn *= sign(message_v2c[c, v])
              sum += phi(abs(message_v2c[c,v]))
            end
          end
          message_c2v[c, variable] = sgn * phi(sum)
        end
      end
    end

    # Message receive by check nodes
    for check in 1:num_cs
      for v in 1:num_vs
        if parity_check_matrix[check, v] == 1
          sum = 0
        
          for c in 1:num_cs
            if c != check && parity_check_matrix[c, v] == 1
              sum += message_c2v[c, v]  
            end
          end
  
          message_v2c[check, v] = sum + vllr[v]  
        end
      end  
    end


    llr = estimate(message_c2v, num_vs, num_cs, vllr, parity_check_matrix)
    decoded = llr .< 0
    syndrome = (parity_check_matrix * decoded) .% 2
    
    if iszero(syndrome)
      break
    end

    
  end

  return decoded, iszero(syndrome)
end
