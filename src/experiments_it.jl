using DelimitedFiles
using PyPlot

function simulate_it(parity_check_matrix, physical_error_rates, max_trials, output_file_path)
  outputs = []
  
  for per in physical_error_rates
    @info "Simulation for per: $per"
    ler = syndrome_it_simulate(parity_check_matrix, per, max_trials)
    append!(outputs, ler)
  end

  to_save = ["IT_syndrome_decoder"; physical_error_rates'; outputs']
  @info "Writing to output file: $output_file_path"
  writedlm(output_file_path, to_save)
end
