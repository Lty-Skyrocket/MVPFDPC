import numpy as np




lnc1 = np.loadtxt(r"Expression/lnc-expression-p.txt", dtype=float)
mir1 = np.loadtxt(r"Expression/mi-expression-p.txt", dtype=float)
lnc2 = np.loadtxt(r"Function/lnc-function-p.txt", dtype=float)
mir2 = np.loadtxt(r"Function/mi-function-p.txt", dtype=float)
lnc3 = np.loadtxt(r"Sequence/lnc-seq-p-weight.txt", dtype=float)
mir3 = np.loadtxt(r"Sequence/mi-seq-p.txt", dtype=float)




lnc = (lnc1 + lnc2 + lnc3 )/3

mir = (mir1 + mir2 + mir3 )/3

np.savetxt(r'lnc.txt', lnc, delimiter='\t', fmt='%.9f')
np.savetxt(r'mir.txt', mir, delimiter='\t', fmt='%.9f')
