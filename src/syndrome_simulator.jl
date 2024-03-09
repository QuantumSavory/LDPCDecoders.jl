function syndrome_simulate(parity_check_matrix, physical_error_rate, max_trials)

  # Get size of parity check matrix
  num_checks, num_bits = size(parity_check_matrix)

  parity_check_matrix_T = sparse(parity_check_matrix')

  suc = 0

  @info "Simulating for $physical_error_rate for $max_trials trials"
  tenths = floor(max_trials/10)

  # Initialization
  log_probabs = zeros(num_bits)
  channel_probabs = fill(physical_error_rate, num_bits)
  b2c = zeros(num_checks, num_bits)
  c2b = zeros(num_checks, num_bits)
  err = zeros(num_bits)

  for i in 1:max_trials

    # Generate error based on error rate
    error = rand(num_bits) .< physical_error_rate

    syndrome = (parity_check_matrix * error) .% 2

    # Belief propagation decoder
    decoded_error, decoded = syndrome_decode(parity_check_matrix, parity_check_matrix_T, syndrome, 10, copy(channel_probabs), copy(b2c), copy(c2b), copy(log_probabs), copy(err))

    if decoded == true
      suc += 1
    end

    if i % tenths == 0
      @info "Completed $i trials"
    end
  end

  # mean_ler = mean(logical_error_rates)
  ler = (max_trials - suc) / max_trials
  println("The logical error rate is : ", ler)
  return ler

end
