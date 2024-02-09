import numpy as np
from ldpc import bposd_decoder
import matplotlib.pyplot as plt
from bposd.css import css_code

Hx = np.load("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_X_rankX120_rankZ179_minWtX2_minWtZ2.npz")
rank = np.linalg.matrix_rank(Hx)
# print("The rank of the matrix is : ", rank)

physicalErrorRates = np.arange(0.00, 0.015, 0.0025)
# print(physicalErrorRates)
maxTrials = 1000

Hx = np.load("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_X_rankX120_rankZ179_minWtX2_minWtZ2.npz")
Hz = np.load("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_Z_rankX120_rankZ179_minWtX2_minWtZ2.npz")
code = css_code(Hx, Hz)
lers = []

# for row in Hx:
#   print(sum(row))


# juliaFilePath = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-bp-osd-Hz.txt"

# with open(juliaFilePath, 'r') as file:
#    data = file.readlines()
   
# plt.figure()
# title = data[0].strip().split('\t')[0]
# pers2 = [float(x) for x in data[1].strip().split('\t')[:-1]]
# lers2 = [float(x) for x in data[2].strip().split('\t')[:-1]]
# print(physicalErrorRates)
# print(pers2)
# print(lers2)

# plt.plot(physicalErrorRates, lers2, label="Julia - LDPCDecoders")
# plt.plot(physicalErrorRates, physicalErrorRates, "k--")


# for per in physicalErrorRates:
#   nchecks, nbits = np.shape(Hz)
#   print(nchecks, nbits)
  
#   bpd=bposd_decoder(
#       code.hz,#the parity check matrix
#       error_rate=per,
#       channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
#       max_iter=1000, #the maximum number of iterations for BP)
#       bp_method="ms",
#       ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
#       osd_method="osd0", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
#       osd_order=7 #the osd search depth
#       )
  
#   suc = 0
#   print(f"Running simulation for per: {per}")
#   for i in range(maxTrials):

#       err = np.random.rand(code.N)
#       # print("length of err: ", len(err))
#       err = err < per
#       syn = code.hz @ err % 2
#       syn = np.array(syn)
#       # print("length of syn: ", len(syn))
#       bpd.decode(syn)

#       residual = (bpd.osdw_decoding + err) % 2
#       a=(code.lz @ residual % 2).any()
#       # if np.array_equal(bpd.osdw_decoding,err): 
#       #     suc+=1
#       if not a:
#          suc += 1


#       if i % 100 == 99:
#           print(f"Completed {i+1} trials")

#   ler = (maxTrials - suc) / maxTrials
#   print(f"The logical error rate is: {ler}")
#   lers.append(ler)


# # plt.title(title)

# plt.xlabel("Physical error rate")
# plt.ylabel("Logical error rate")

# plt.plot(physicalErrorRates, lers, label="Python - ldpc")
# # plt.plot(physicalErrorRates, physicalErrorRates, '--')
# plt.legend()
# plt.title("Python vs Julia decoder")
# plt.show()






