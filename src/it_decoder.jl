# Bit flipping algorithm
# Iterative decoder
function it_decode(parity_check_matrix, received_message, iterations)

  num_checks, num_bits = size(parity_check_matrix)
  
  syndrome = (received_message * transpose(parity_check_matrix)) .% 2 

  # No of iterations
  for iter in 1:iterations
    if iszero(syndrome)
      break
    end

    votes = zeros(num_bits, 2)
    for check in 1:num_checks
      indices = zeros(num_bits)
      sum = 0
      for bit in 1:num_bits
        if parity_check_matrix[check, bit] == 1
          sum += parity_check_matrix[check, bit] * received_message[bit]
          indices[num_bits] = 1
        end
      end

      sum = sum % 2
      for bit in 1:num_bits
        if indices[bit] == 1
          votes[bit, sum+1] += 1
        end
      end
    end

    for bit in 1:num_bits
      if votes[bit, 1] < votes[bit, 2]
        received_message[bit] = (received_message[bit] + 1) % 2
      end
    end

    syndrome = (received_message * transpose(parity_check_matrix)) .% 2 
  end
  
  return received_message
end

