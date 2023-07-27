using DelimitedFiles
function simulate_bp(parity_check_matrix, physical_error_rates,  max_trials, output_file_path)
  outputs = []
  
  for per in physical_error_rates
    @info "Simulation for per: " + per
    ler = syndrome_simulate(parity_check_matrix, per, max_trials)
    outputs = [outputs; [per ler]]
  end

  @info "Writing to output file: " + output_file_path
  writedlm(output_file_path, outputs)
end