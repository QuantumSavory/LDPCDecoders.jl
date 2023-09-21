using LinearAlgebra
using Random
using Statistics


function syndrome_it_simulate(parity_check_matrix, physical_error_rate, max_trials)

  # Get size of parity check matrix
  num_checks, num_bits = size(parity_check_matrix)
  
  suc = 0
  logical_error_rates = []
  @info "Simulating for $physical_error_rate for $max_trials trials"
  tenths = floor(max_trials/10)

  # Initalization
  err = zeros(Bool, num_bits)
  votes = zeros(Int, num_bits)
  curr_syn = zeros(Bool, num_checks)
  error_checks = zeros(Bool, num_checks)

  for i in 1:max_trials
   
    # Generate error 
    error = rand(num_bits) .< physical_error_rate
  
    # Calculate syndrome
    syndrome = (parity_check_matrix * error) .% 2

    # Iterative bit flip decoder
    decoded_error, decoded = syndrome_it_decode(parity_check_matrix, syndrome, 100, copy(err), copy(votes))
    
    if decoded == true
      suc += 1
    end

    if i % tenths == 0
      @info "Completed $i trials"
    end
  end
  
  
  ler = (max_trials - suc) / max_trials
  println("The logical error rate is : ", ler)
  return ler

end

