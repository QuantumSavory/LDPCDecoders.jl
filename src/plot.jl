using PyPlot

function plot_per_vs_ler(file_path, log_scale = false, save_file_path = "")
  
  # Read data from file
  data = readdlm(file_path)

  title = data[1,1]
  pers = data[2, :]
  lers = data[3, :]

  if log_scale
    pers = log10.(pers)
    lers = log10.(lers)
    # plt.loglog(pers, lers, label="$title")
    # plt.loglog(pers, pers, "k--")
  end

  plt.plot(pers, lers, label="$title")
  plt.plot(pers, pers, "k--")
  # plt.title(title)
  plt.legend()
  plt.xlabel("Physical error rate")
  plt.ylabel("Logical error rate")
  if isempty(save_file_path)  
    plt.show()
  else
    plt.savefig(save_file_path)
  end
  

end