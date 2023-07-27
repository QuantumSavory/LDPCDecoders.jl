using LinearAlgebra
using Random
using Statistics


function syndrome_simulate(parity_check_matrix, physical_error_rate, max_trials)

  # Get size of parity check matrix
  num_checks, num_bits = size(parity_check_matrix)
  

  success = 0
  logical_error_rates = []
  for i in 1:max_trials
    println("Trial number : ", i)
   
    # Generate error based on error rate
    error = zeros(Bool, num_bits)
    bits = rand(1:num_bits, Int(physical_error_rate * num_bits))
    for bit in bits
      error[bit] = 1
    end
  
    # println("error = ", error)
    syndrome = (parity_check_matrix * error) .% 2
    # println("syndrome = ", syndrome)

    # Belief propogation decoder
    decoded_error, decoded = syndrome_decode(parity_check_matrix, syndrome, physical_error_rate, 100)
    # println("decoded error = ", decoded_error, " success = ", decoded)

    # if all(Int.(decoded_error) .== error)
    #   success += 1
    # end

    ler = sum(decoded_error .!= error) / num_bits
    logical_error_rates = [logical_error_rates; ler]
    
  end
  
  mean_ler = mean(logical_error_rates)
  println("The average logical error rate is : ", mean_ler)
  return mean_ler

end

