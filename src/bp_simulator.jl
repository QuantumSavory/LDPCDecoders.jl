function bp_simulate(parity_check_matrix, generator_matrix, error_rate, max_trials)

  # Get size of parity check matrix
  num_checks, num_bits = size(parity_check_matrix)
  # generator_matrix = parity_to_generator(parity_check_matrix)
  message_size = num_bits - num_checks
  display(generator_matrix')

  success = 0
  for i in 1:max_trials
    println("******************* Iteration number : ", i)

    # Generate error based on error rate
    message = bitrand(message_size)
    println("message = ", message)
    # code = (generator_matrix * message) .% 2
    code = (transpose(message) * generator_matrix) .% 2
    println("code = ", code)
    received_message = code
    println("received message = ", received_message)

    error = rand(num_bits)
    for j in 1:num_bits

      if error[j] < error_rate
        error[j] = 1
      else
        error[j] = 0
      end
    end


    error = Int.(error)
    println("error = ", error)
    received_message = vec(code) .⊻ error
    println("received message2 = ", received_message)
    display(received_message)

    # syndrome = (pcm * error) .% 2

    # Belief propogation decoder
    decoded_message, success = bp_decode(received_message, parity_check_matrix, 0.1)
    println("decoded message = ", decoded_message)
    display(decoded_message)
    if decoded_message == code
      success += 1
    end

    # n = size(decoded_message)
    # correct = true
    # for k in 1:n
    #   if decoded_message[i] != code[i]
    #     correct = false
    #   end
    # end

    # if correct
    #   success += 1
    # end



  end
  println("The success rate of the decoder is ")
  println(success/max_trials)


end
