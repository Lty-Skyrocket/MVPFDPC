import numpy as np




Sim1 = np.loadtxt(r"Sequence/lnc_seq_similarity.txt", dtype=float)
Sim2 = np.loadtxt(r"Sequence/mi_seq_similarity.txt", dtype=float)
Inv1 = np.loadtxt(r"Sequence/invalid_lnc_seq.txt", dtype=int)





# 避免NaN元素
Sim1[Inv1, :] = 0
Sim1[:, Inv1] = 0



# 计算均值
mean_lnc = np.mean(Sim1)
mean_mi = np.mean(Sim2)

# 避免NaN元素
Sim1[Inv1, :] = mean_lnc
Sim1[:, Inv1] = mean_lnc



#归一化

def normalize(in_img):
    # all elements divide the max
    A = np.tril(in_img)
    for i in range(A.shape[1]):
        A[i:, i] = A[i:, i] / A[i, i]
    C = np.tril(A, k=-1)
    out_img = C.T + A
    return out_img

Sim1=normalize(Sim1);
Sim2=normalize(Sim2);

def set_negative_to_zero(matrix):
    matrix[matrix < 0] = 0
    return matrix

Sim1 = set_negative_to_zero(Sim1)


np.savetxt(r'lnc-seq-p-weight.txt', Sim1, delimiter='\t', fmt='%.9f')
np.savetxt(r'mi-seq-p.txt', Sim2, delimiter='\t', fmt='%.9f')

# np.savetxt(r'lnc-expression-p.txt', Sim1, delimiter='\t', fmt='%.9f')
# np.savetxt(r'mi-expression-p.txt', Sim2, delimiter='\t', fmt='%.9f')

