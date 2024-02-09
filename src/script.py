import numpy as np
from ldpc import bposd_decoder, bp_decoder
import matplotlib.pyplot as plt
from bposd.css import css_code
import csv

## Import the Hx, Hz matrices
Hx = np.load("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_X_rankX120_rankZ179_minWtX2_minWtZ2.npz")
Hz = np.load("/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/codes_for_hardware_test 2/1_ra1_rb2_Z_rankX120_rankZ179_minWtX2_minWtZ2.npz")
code = css_code(Hx, Hz)
# rank = np.linalg.matrix_rank(Hx)

# Define physical error rates
# physicalErrorRates = np.logspace(np.log10(0.0001), np.log10(0.5), 10) # logspaced pers
physicalErrorRates = np.linspace(0.001, 0.2, 10) # linearly spaces pers



maxTrials = 10000
lers = []

H = np.load("/Users/krishnapg/Desktop/umass/masters-project/comparisons/pcm.npy")
# juliaFilePath = "/Users/krishnapg/Desktop/umass/masters-project/pauli-frame/output/output-bp-osd-Hz.txt"
# juliaFilePath = "/Users/krishnapg/Desktop/umass/masters-project/steane_simulation.txt"
# with open(juliaFilePath, 'r') as file:
#    data = file.readlines()
   
# plt.figure()
# title = data[0].strip().split('\t')[0]
# pers2 = [float(x) for x in data[1].strip().split('\t')[:-1]]
# lers2 = [float(x) for x in data[2].strip().split('\t')[:-1]]
# print(physicalErrorRates)
# print(pers2)
# print(lers2)

logPers = np.log10(physicalErrorRates)
# logPers2 = np.log10(pers2)
# logLers2 = np.log10(lers2)
# plt.plot(logPers2, logLers2, label="Julia - LDPCDecoders")
plt.plot(logPers, logPers, "k--")

for per in physicalErrorRates:
  nchecks, nbits = np.shape(H)
  print(nchecks, nbits)
  
  # bpd=bposd_decoder(
  #     code.hz,#the parity check matrix
  #     error_rate=per,
  #     max_iter=Hz.shape[1], #the maximum number of iterations for BP)
  #     bp_method=1,
  #     )

  bpd=bp_decoder(
      H,#the parity check matrix
      error_rate=per,
      max_iter=10, #the maximum number of iterations for BP)
      bp_method=2,
      )
  
  suc = 0
  print(f"Running simulation for per: {per}")
  for i in range(maxTrials):

      err = np.random.rand(nbits)
      # print("length of err: ", len(err))
      err = err < per
      err = err.astype(int)
      # err[err>=per] = int(0)
      # err[err>0] = int(1)
      

      # err2 = np.random.rand(code.N)
      # err2 = err2 < per 
      # print(err)
      
      # err2 = err2.astype(int)
      # print("err2", err2)
      syn = H @ err % 2
      syn = np.array(syn)
      # print("length of syn: ", len(syn))
      decoding = bpd.decode(syn)

    #   residual = (bpd.osdw_decoding + err) % 2
    #   a=(code.lz @ residual % 2).any()
      

      if abs(decoding - err).sum() == 0: 
          suc+=1
      # if not a:
      #    suc += 1

      if i % 100 == 99:
          print(f"Completed {i+1} trials")

  ler = (maxTrials - suc) / maxTrials
  print(f"The logical error rate is: {ler}")
  lers.append(ler)

pythonHzFilePath = "/Users/krishnapg/Desktop/umass/masters-project/comparisons-2/python-classical-bp.txt"
logLers = np.log10(lers)
with open(pythonHzFilePath, 'w') as file: 
   csvwriter = csv.writer(file, delimiter=' ')

   csvwriter.writerow(['Python'])
   csvwriter.writerow(physicalErrorRates)
   csvwriter.writerow(lers)


# plt.title(title)

plt.xlabel("Physical error rate")
plt.ylabel("Logical error rate")
plt.plot(logPers, logLers, label="Python - ldpc")
# plt.plot(physicalErrorRates, physicalErrorRates, '--')
plt.legend()
plt.title("Python vs Julia decoder")
plt.show()






