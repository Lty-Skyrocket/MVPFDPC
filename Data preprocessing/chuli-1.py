import numpy as np



Sim1 = np.loadtxt(r"Expression/lnc_expression_similarity.txt", dtype=float)
Sim2 = np.loadtxt(r"Expression/mi_expression_similarity.txt", dtype=float)
Inv1 = np.loadtxt(r"Expression/invalid_lnc_expression.txt", dtype=int)
Inv2 = np.loadtxt(r"Expression/invalid_mi_expression.txt", dtype=int)

# Sim1 = np.loadtxt(r"Function/lnc_function_similarity.txt", dtype=float)
# Sim2 = np.loadtxt(r"Function/mi_function_similarity.txt", dtype=float)
# Inv1 = np.loadtxt(r"Function/invalid_lnc_function.txt", dtype=int)
# Inv2 = np.loadtxt(r"Function/invalid_mi_function.txt", dtype=int)


def subtract_one(matrix):
    return matrix - 1

Inv1 = subtract_one(Inv1)
Inv2 = subtract_one(Inv2)



# 避免NaN元素
Sim1[Inv1, :] = 0
Sim1[:, Inv1] = 0
Sim2[Inv2, :] = 0
Sim2[:, Inv2] = 0


# 计算均值
mean_lnc = np.mean(Sim1)
mean_mi = np.mean(Sim2)

# 避免NaN元素
Sim1[Inv1, :] = mean_lnc
Sim1[:, Inv1] = mean_lnc
Sim2[Inv2, :] = mean_mi
Sim2[:, Inv2] = mean_mi

np.savetxt(r'lnc-expression-p.txt', Sim1, delimiter='\t', fmt='%.9f')
np.savetxt(r'mi-expression-p.txt', Sim2, delimiter='\t', fmt='%.9f')
#
# np.savetxt(r'lnc-function-p.txt', Sim1, delimiter='\t', fmt='%.9f')
# np.savetxt(r'mi-function-p.txt', Sim2, delimiter='\t', fmt='%.9f')