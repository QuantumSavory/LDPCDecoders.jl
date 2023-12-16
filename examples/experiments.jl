using DelimitedFiles
using PyPlot
using LDPCDecoders

function simulate_bp(parity_check_matrix, physical_error_rates,  max_trials, output_file_path)
  outputs = []

  for per in physical_error_rates
    @info "Simulation for per: $per"
    ler = syndrome_simulate(parity_check_matrix, per, max_trials)
    append!(outputs, ler)
  end

  to_save = ["BP_syndrome_decoder"; physical_error_rates'; outputs']
  @info "Writing to output file: $output_file_path"
  writedlm(output_file_path, to_save)
end

function plot_per_vs_ler(file_path)

  # Read data from file
  data = readdlm(file_path)

  title = data[1,1]
  pers = data[2, :]
  lers = data[3, :]

  # pers = log.(data[2, :])
  # lers = log.(data[3, :])

  plt.loglog(pers, lers, label="$title")
  plt.loglog(pers, pers, "k--")
  # plt.title(title)
  plt.legend()
  plt.xlabel("Physical error rate")
  plt.ylabel("Logical error rate")
  plt.show()

end
