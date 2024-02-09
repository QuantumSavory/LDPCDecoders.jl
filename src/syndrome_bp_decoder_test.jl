using NPZ
using SparseArrays

# Hx = npzread("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_X_rankX120_rankZ179_minWtX2_minWtZ2.npz")
# Hz = npzread("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_Z_rankX120_rankZ179_minWtX2_minWtZ2.npz")

# Hx = npzread("/Users/krishnapg/Desktop/umass/masters-project/steane_hx.npz")
# Hz = npzread("/Users/krishnapg/Desktop/umass/masters-project/steane_hz.npz")

# println(typeof(Hx))
# rankHx = rank(Hx)
# println("The rank of the Hx matrix is : $rankHx")
# Hx = BitArray(Hx)
# Hx = Hx["arr_0"]
# Hz = Hz["arr_0"]
# display(Hx["arr_0"])
# display(Hz)
# physicalErrorRates = collect(0.0:0.0025:0.015)
physicalErrorRates = range(0.001, 0.2, 10)

H = npzread("/Users/krishnapg/Desktop/umass/masters-project/comparisons/pcm.npy" )
H = BitArray(H)
H = sparse(H)
filePathBP = "/Users/krishnapg/Desktop/umass/masters-project/comparisons-2/bp-classical.txt" 
filePathBPOSD = "/Users/krishnapg/Desktop/umass/masters-project/comparisons-2/bp-osd-classical.txt" 
# file_path = "/Users/krishnapg/Desktop/umass/masters-project/Hx_simulation.txt"

# HxFilePath = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-Hx.txt"
# HxSaveFilePath = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-plot-Hx.png"

# HzFilePath = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-Hz.txt"
# HzSaveFilePath = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-plot-Hz.png"

# HxFilePathBpOsd = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-bp-osd-Hx.txt"
# HxSaveFilePathBpOsd = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-plot-bp-osd-Hx.png"

# HzFilePathBpOsd = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-bp-osd-Hz.txt"
# HzSaveFilePathBpOsd = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-plot-bp-osd-Hz.png"

# Hx = Float16.(Hx)
# Hx = sparse(Hx)
# simulate_bp(Hx, physicalErrorRates, 100, HxFilePath)
# plot_per_vs_ler(HxFilePath)
# plot_per_vs_ler(HxFilePath, HxSaveFilePath)

# Hz = sparse(Hz)
# simulate_bp(Hz, physicalErrorRates, 10000, HzFilePath)
# # plot_per_vs_ler(HzFilePath)
# plot_per_vs_ler(HzFilePath, HzSaveFilePath)

# simulate_bp_osd(Hx, physicalErrorRates, 10000, HxFilePathBpOsd)
# plot_per_vs_ler(HxFilePathBpOsd)

# simulate_bp_osd(H, physicalErrorRates, 1000, filePathBPOSD)
# simulate_bp(H, physicalErrorRates, 10000, filePathBP)
# plot_per_vs_ler(HzFilePathBpOsd)
plot_per_vs_ler(filePathBP, true)
