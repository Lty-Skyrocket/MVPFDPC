

import numpy as np
#import pandas as pd
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import matplotlib.pyplot as plt
import copy
import random_sample



# SM = np.loadtxt(r"lnc.txt", dtype=float)
# MM = np.loadtxt(r"mir.txt", dtype=float)

lncRNA = np.loadtxt(r"lnc_integation.txt", dtype=float)
miRNA = np.loadtxt(r"mir_integation.txt", dtype=float)




# Y1= np.loadtxt(r"interaction.txt",dtype=float)
Y1= np.loadtxt(r"MVPF_association.txt",dtype=float)
lncRNA_miRNA_k = np.loadtxt(r"known.txt",dtype=int)
lncRNA_miRNA_uk = np.loadtxt(r"unknown.txt",dtype=int)
selected_columns_sorted6 = np.loadtxt(r"selected_columns_sorted.txt", dtype=int)
selected_rows_sorted6 = np.loadtxt(r"selected_rows_sorted.txt", dtype=int)

# # 随机取796行
# row_rand_array = np.arange(SM_miRNA_uk.shape[0])
# np.random.shuffle(row_rand_array)
# SM_miRNA_uk_5118.txt = SM_miRNA_uk[row_rand_array[0:5118]]
# np.savetxt(r'SM_miRNA_uk_5118.txt', SM_miRNA_uk_5118.txt, delimiter='\t', fmt='%d')

#分治矩阵
# 选取的列数和行数
lty = 250  #列数
dty = 500  #行数
# 随机选取l列
# selected_columns = np.random.choice(Y1.shape[1], lty, replace=False)
# selected_columns_sorted = sorted(selected_columns)  # 对列索引进行排序
# np.savetxt(r'selected_columns_sorted.txt', selected_columns_sorted, delimiter='\t', fmt='%d')

# # # # 随机选取d行
# selected_rows = np.random.choice(Y1.shape[0], dty, replace=False)
# selected_rows_sorted = sorted(selected_rows)  # 对列索引进行排序
# np.savetxt(r'selected_rows_sorted.txt', selected_rows_sorted, delimiter='\t', fmt='%d')
# #

# 将数组转换为列表
selected_columns_sorted = selected_columns_sorted6.flatten().tolist()
selected_rows_sorted = selected_rows_sorted6.flatten().tolist()

def run_MC_1(P):


    # 生成子矩阵
    submatrix_1, submatrix_2, NSM, NMM = random_sample.generate_submatrix(P, lty, dty,selected_columns_sorted,selected_rows_sorted)


    # # 选取的列数和行数
    # l = 200
    # d = 200
    # # 生成子矩阵
    # submatrix_1, submatrix_2, NSM, NMM = random_sample.generate_submatrix(P, l, d)
    #经过矩阵分解后的两个新子矩阵
    submatrix_n1, submatrix_n2 = run_MC_2(submatrix_1, submatrix_2, NSM, NMM)


    # print(submatrix_n1, submatrix_n2)
    #生成交叉矩阵并且SVD求伪逆
    intersection_matrix = random_sample.intersection(dty,lty,selected_rows_sorted,selected_columns_sorted)

    # 使用 SVD 对矩阵 A 进行奇异值分解
    U, sigma, Vt = np.linalg.svd(intersection_matrix)
    r3 = 50
    Wt = np.zeros([r3,r3])
    for i in range(0,r3):
        Wt[i][i]=sigma[i]
    U= U[:, 0:r3]
    Vt= Vt[0:r3,:]

    # 计算伪逆的步骤
    # 首先将 sigma 转换为对角矩阵
    Sigma_plus = np.zeros((Wt.shape[0], Wt.shape[1]))
    # Sigma[:min(intersection_matrix.shape[0], intersection_matrix.shape[1]), :min(intersection_matrix.shape[0], intersection_matrix.shape[1])] = np.diag(sigma)
    # Sigma_plus=np.zeros_like(Sigma)

    # 对非零奇异值取倒数
    for i in range (Wt.shape[0]):
        if Wt[i][i]!=0:
            Sigma_plus[i][i] = 1/Wt[i][i]   # 求解伪逆矩阵
    A_plus = np.dot(Vt.T, np.dot(Sigma_plus, U.T))




    #计算最终结果
    LTY = np.dot(submatrix_n1,np.dot(A_plus,submatrix_n2))

    return LTY

#矩阵分解

def TCMF(alpha, beta,gamma, Y, maxiter,A,B,C,SM,MM):

    iter0=1
    while True:

        a = np.dot(Y,B)+beta*np.dot(SM,A)
        b = np.dot(np.transpose(B),B)+alpha*C+beta*np.dot(np.transpose(A),A)
        c = np.dot(np.transpose(Y),A)+gamma*np.dot(MM,B)
        d = np.dot(np.transpose(A), A) + alpha * C + gamma * np.dot(np.transpose(B), B)

        A = np.dot(a,np.linalg.inv(b))
        B = np.dot(c, np.linalg.inv(d))

        if iter0 >= maxiter:

            #print('reach maximum iteration!')
            break
        iter0 = iter0 + 1
    Y= np.dot(A,np.transpose(B))
    Y_recover = Y
    return Y_recover






def run_MC_2(Y,L,NSM, NMM):
    maxiter = 500
    alpha = 0.30
    beta = 0.0001
    gamma = 35

#SVD
#对矩阵X_2的奇异值分解

    U, S, V = np.linalg.svd(Y)
    S1=np.sqrt(S)
    r1 = 30
    Wt1 = np.zeros([r1,r1])
    for i in range(0,r1):
        Wt1[i][i]=S1[i]
    U= U[:, 0:r1]
    V= V[0:r1,:]
    A = np.dot(U,Wt1)
    B1 = np.dot(Wt1,V)
    B=np.transpose(B1)
    C = np.zeros([r1,r1])
    for i in range(0,r1):
        C[i][i] = 1

#对矩阵X_2的奇异值分解
    U, S, V = np.linalg.svd(L)
    S2 = np.sqrt(S)
    r2 = 30
    Wt2 = np.zeros([r2, r2])
    for i in range(0, r2):
        Wt2[i][i] = S2[i]
    U = U[:, 0:r2]
    V = V[0:r2, :]
    A1 = np.dot(U, Wt2)
    B2 = np.dot(Wt2, V)
    B1 = np.transpose(B2)
    C1= np.zeros([r2, r2])
    for i in range(0, r2):
        C1[i][i] = 1
    Y = TCMF(alpha, beta,gamma,Y, maxiter,A,B,C,SM,NSM)

    L = TCMF(alpha, beta, gamma, L, maxiter, A1, B1, C1, NMM, MM)


    return Y,L







def main():
    roc_sum, time = 0, 0
    all_fpr, all_tpr, all_auc = [], [], []
    # pr_sum, time = 0, 0
    kf = KFold(n_splits=5, shuffle=True, random_state=4444)#//初始化kfold
    for train_index,test_index in kf.split(lncRNA_miRNA_k):
        X_2 = copy.deepcopy(Y1)

        for index in test_index:
            X_2[lncRNA_miRNA_k[index, 0] , lncRNA_miRNA_k[index, 1] ] = 0


        M_1 = run_MC_1(X_2)
        Label = np.zeros(lncRNA_miRNA_uk.shape[0] + test_index.size)
        Score = np.zeros(lncRNA_miRNA_uk.shape[0] + test_index.size)
        i , j = 0 , 0
        for s_index in test_index:
            Label[i] = 1
            Score[i] = M_1[lncRNA_miRNA_k[s_index,0],lncRNA_miRNA_k[s_index,1]]

            i = i + 1
        for i in range(test_index.size, lncRNA_miRNA_uk.shape[0] + test_index.size):
            Score[i] = M_1[lncRNA_miRNA_uk[j,0],lncRNA_miRNA_uk[j,1]]
            j = j + 1
        fpr, tpr, thersholds = roc_curve(y_true=Label, y_score=Score, drop_intermediate=False)
        roc_auc = auc(fpr, tpr)
        roc_sum = roc_sum + roc_auc
        time += 1
        s=roc_sum/time
        print(f'Fold {time}: AUC = {roc_auc:.4f}, Cumulative AUC = {s:.4f}')

        # Store values for ROC curve
        all_fpr.append(fpr)
        all_tpr.append(tpr)
        all_auc.append(roc_auc)

    # Plot ROC curves for each fold
    plt.figure(figsize=(8, 6))

    for i in range(len(all_fpr)):
        plt.plot(all_fpr[i], all_tpr[i], label=f'ROC fold {i + 1} (AUC = {all_auc[i]:.4f})', linestyle='-',linewidth=1)
    # 先找到最小的长度
    min_length = min(len(fpr) for fpr in all_fpr)

    # 对每个子数组进行截断或插值，使它们具有相同的长度
    all_fpr_fixed = [np.interp(np.linspace(0, 1, min_length), fpr, fpr) for fpr in all_fpr]
    all_tpr_fixed = [np.interp(np.linspace(0, 1, min_length), tpr, tpr) for tpr in all_tpr]

    # 然后计算平均值
    mean_fpr = np.mean(all_fpr_fixed, axis=0)
    mean_tpr = np.mean(all_tpr_fixed, axis=0)
    mean_auc = np.mean(all_auc)
    np.savetxt(r'mean_fpr.txt', mean_fpr, delimiter='\t', fmt='%.9f')
    np.savetxt(r'mean_tpr.txt', mean_tpr, delimiter='\t', fmt='%.9f')


    plt.plot(fpr, tpr, label=f'Mean ROC (AUC = {mean_auc:.4f})', linestyle='-')
    # plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve for 5-Fold CV')
    plt.legend(loc='lower right')
    # 保存图像到文件
    plt.savefig('roc_5fold.png')
    plt.show()




if __name__ == "__main__":
    main()

