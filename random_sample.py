import numpy as np





Y1= np.loadtxt(r"MVPF_association.txt",dtype=float)

SM = np.loadtxt(r"lnc_integration.txt", dtype=float)
MM = np.loadtxt(r"mir_integation.txt", dtype=float)


selected_columns_sorted6 = np.loadtxt(r"selected_columns_sorted.txt", dtype=float)
selected_rows_sorted6 = np.loadtxt(r"selected_rows_sorted.txt", dtype=float)
# 选取的列数和行数
# l = 200
# d = 200
#生成子矩阵
# 随机选取l列
# selected_columns = np.random.choice(Y1.shape[1], l, replace=False)
# selected_columns_sorted = sorted(selected_columns)  # 对列索引进行排序
# np.savetxt(r'selected_columns_sorted.txt', selected_columns_sorted, delimiter='\t', fmt='%d')
#
# # 随机选取d行
# selected_rows = np.random.choice(Y1.shape[0], d, replace=False)
# selected_rows_sorted = sorted(selected_rows)  # 对列索引进行排序
# np.savetxt(r'selected_rows_sorted.txt', selected_rows_sorted, delimiter='\t', fmt='%d')



# selected_columns_sorted = []
#
# for x in selected_columns_sorted6:
#     selected_columns_sorted.append(x)
# selected_rows_sorted = []
#
# for x in selected_rows_sorted6:
#     selected_rows_sorted.append(x)
def generate_submatrix(matrix, l, d,selected_columns_sorted,selected_rows_sorted):


    submatrix_1 = matrix[:, selected_columns_sorted]



    submatrix_2 = matrix[selected_rows_sorted, :]
    # print(submatrix_1)
    # print(submatrix_2)
    # # step1:处理miRNA相似性

    NSM = np.zeros((l, l))

    for i in range(l):
        for j in range(i, l):
            NSM[i][j] = MM[selected_columns_sorted[i]][selected_columns_sorted[j]]
            if i != j:
                NSM[j][i] = MM[selected_columns_sorted[i]][selected_columns_sorted[j]]  # 对称位置也需要填充

    # step2:处理SM相似性
    NMM = np.zeros((d, d))
    for i in range(d):
        for j in range(i, d):
            NMM[i][j] = SM[selected_rows_sorted[i]][selected_rows_sorted[j]]
            if i != j:
                NMM[j][i] = SM[selected_rows_sorted[i]][selected_rows_sorted[j]]  # 对称位置也需要填充

    return submatrix_1,submatrix_2,NSM,NMM,

def intersection(d,l,selected_rows_sorted,selected_columns_sorted):
    NJX = np.zeros((d, l))
    for i in range(d):
        for j in range(l):
            NJX[i][j] = Y1[selected_rows_sorted[i]][selected_columns_sorted[j]]

    # return NJX
    return NJX









