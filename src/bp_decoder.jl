# Bi-Symmetric Channel LLR (Bit flip)
function llr(bit, error_rate)
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

function send_variable_message!(parity_check_matrix, message_c2v, message_v2c, num_cs, num_vs, vllr)
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
end

function send_check_message!(parity_check_matrix, message_c2v, message_v2c, num_cs, num_vs)
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
end

function initialise_checks!(pcm, message_c2v, syndrome, num_vs, num_cs)
  for check in 1:num_cs
    msg = syndrome[check]
    for variable in 1:num_vs
      if pcm[check, variable] == 1
        message_c2v = msg
      end
    end
  end
end


# Belief Propagation Decoder
function bp_decode(parity_check_matrix, received_message, error_rate, max_iterations=100)

  num_checks, num_bits = size(parity_check_matrix)
  num_edges = sum(parity_check_matrix)

  num_cs, num_vs = size(parity_check_matrix)

  # Initialize messages
  message_v2c = zeros(num_checks, num_bits)

  # Initialize check messages
  message_c2v = zeros(num_checks, num_bits)
  initialise_checks!(parity_check_matrix, message_c2v, syndrome, num_vs, num_cs)

  # Intialize llr for variable nodes
  vllr = llr.(received_message, error_rate)

  # Send message from variable to check nodes
  send_variable_message!(parity_check_matrix, message_c2v, message_v2c, num_checks, num_bits, vllr)

  local decoded
  local syndrome

  for iter in 1:max_iterations

    # Send Message from check to variable nodes
    send_check_message!(parity_check_matrix, message_c2v, message_v2c, num_checks, num_bits)

    # Send message from variable to check nodes
    send_variable_message!(parity_check_matrix, message_c2v, message_v2c, num_checks, num_bits, vllr)

    # Estimate LLR and calculate Syndrome to check H.c+e = 0
    llr = estimate(message_c2v, num_vs, num_cs, vllr, parity_check_matrix)
    decoded = llr .< 0
    syndrome = (parity_check_matrix * decoded) .% 2

    # Check if decoded is a codeword
    if iszero(syndrome)
      break
    end

  end

  return decoded, iszero(syndrome)
end
