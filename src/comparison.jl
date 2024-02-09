using LDPCDecoders
using NPZ
using SparseArrays

# Hx = npzread("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_X_rankX120_rankZ179_minWtX2_minWtZ2.npz")
# Hz = npzread("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_Z_rankX120_rankZ179_minWtX2_minWtZ2.npz")

H = parity_check_matrix(1000, 10, 9)
filePathPCM = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/pcm.npy" 
npzwrite(filePathPCM, H)

physicalErrorRates = range(0.001, 0.2, 10)
maxTrials = 5000


#################### Classical ####################
filePathIterative = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/iterative-classical.txt" 
# simulate_it(H, physicalErrorRates, maxTrials, filePathIterative)

filePathBP = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/bp-classical.txt" 
H = sparse(H)
# simulate_bp(H, physicalErrorRates, maxTrials, filePathIterative)

filePathBPOSD = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/bp-osd-classical.txt" 
simulate_bp_osd(H, physicalErrorRates, maxTrials, filePathIterative)

#################### Quantum ####################
Hx = npzread("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_X_rankX120_rankZ179_minWtX2_minWtZ2.npz")
Hz = npzread("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_Z_rankX120_rankZ179_minWtX2_minWtZ2.npz")

#### Hx
H = Hx 
filePathIterative = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/iterative-hx.txt" 
# simulate_it(H, physicalErrorRates, maxTrials, filePathIterative)

H = sparse(H)
filePathBP = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/bp-hx.txt" 
simulate_bp(H, physicalErrorRates, maxTrials, filePathIterative)

filePathBPOSD = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/bp-osd-hx.txt" 
simulate_bp_osd(H, physicalErrorRates, maxTrials, filePathIterative)

#### Hz
H = Hz
filePathIterative = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/iterative-hz.txt" 
# simulate_it(H, physicalErrorRates, maxTrials, filePathIterative)

H = sparse(H)
filePathBP = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/bp-hz.txt" 
simulate_bp(H, physicalErrorRates, maxTrials, filePathIterative)

filePathBPOSD = "/Users/krishnapg/Desktop/umass/masters-project/comparisons/bp-osd-hz.txt" 
simulate_bp_osd(H, physicalErrorRates, maxTrials, filePathIterative)