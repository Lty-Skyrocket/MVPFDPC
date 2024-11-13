import numpy as np



lnc1 = np.loadtxt(r"lnc.txt", dtype=float)
mir1 = np.loadtxt(r"mir.txt", dtype=float)
lnc2 = np.loadtxt(r"GIP/lncRNA_GIP.txt", dtype=float)
mir2 = np.loadtxt(r"GIP/miRNA_GIP.txt", dtype=float)



lnc_GIP = 0.1*lnc1 + 0.9*lnc2

mir_GIP = 0.1*mir1 + 0.9*mir2

np.savetxt(r'lnc_integration.txt', lnc_GIP, delimiter='\t', fmt='%.9f')
np.savetxt(r'mir_integration.txt', mir_GIP, delimiter='\t', fmt='%.9f')
