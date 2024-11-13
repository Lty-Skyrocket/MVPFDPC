import numpy as np


lncRNA2 = np.loadtxt(r'lncRNA_GIP.txt', dtype=float)
miRNA2 = np.loadtxt(r'miRNA_GIP.txt', dtype=float)

lncRNA3 = np.loadtxt(r'lnc-expression-p.txt', dtype=float)
miRNA3 = np.loadtxt(r'mi-expression-p.txt', dtype=float)

lncRNA4 = np.loadtxt(r'lnc-function-p.txt', dtype=float)
miRNA4 = np.loadtxt(r'mi-function-p.txt', dtype=float)

lncRNA5 = np.loadtxt(r'lnc-seq-p.txt', dtype=float)
miRNA5 = np.loadtxt(r'mi-seq-p.txt', dtype=float)


Y = np.loadtxt(r'interaction.txt',dtype=float)


#小分子投影空间
def fMVP_miRNA(Y, miRNA):

    m, n = Y.shape
    MVP_sm = np.zeros((m, n))



    for i in range(m):
        for j in range(n):  # Adjusted loop to only iterate up to i
            dot_product = np.dot(miRNA[i],Y[:, j])
            mod_j = np.linalg.norm(Y[:, j])
            MVP_sm[i,j] = dot_product / mod_j


    return MVP_sm

#lncRNA投影空间
def fMVP_MM(Y, MM):

    m, n = Y.shape
    MVP_mm = np.zeros((m, n))



    for i in range(m):
        for j in range(n):  # Adjusted loop to only iterate up to i
            dot_product = np.dot(Y[i],MM[:, j])
            mod_j = np.linalg.norm(Y[i])
            MVP_mm[i,j] = dot_product / mod_j



    return MVP_mm



def main():
    m, n = Y.shape

    MVP_association2 = np.zeros((m, n))
    MVP_association3 = np.zeros((m, n))
    MVP_association4 = np.zeros((m, n))
    MVP_association5 = np.zeros((m, n))


    # miRNA投影空间

    MVP_sm2 = fMVP_miRNA(Y, miRNA2)
    MVP_sm3 = fMVP_miRNA(Y, miRNA2)
    MVP_sm4 = fMVP_miRNA(Y, miRNA2)
    MVP_sm5 = fMVP_miRNA(Y, miRNA2)


    #lncRNA投影空间
    MVP_mm2 = fMVP_MM(Y, lncRNA2)
    MVP_mm3 = fMVP_MM(Y, lncRNA2)
    MVP_mm4 = fMVP_MM(Y, lncRNA2)
    MVP_mm5 = fMVP_MM(Y, lncRNA2)




    #对关联矩阵的预处理：
    for i in range(m):
        for j in range(n):  # Adjusted loop to only iterate up to i
            MVP_association2[i,j] = (MVP_sm2[i,j]+MVP_mm2[i,j])/(np.linalg.norm(miRNA2[i])+np.linalg.norm(lncRNA2[:, j]))




    # # #对关联矩阵的预处理：
    for i in range(m):
        for j in range(n):  # Adjusted loop to only iterate up to i
            MVP_association3[i,j] = (MVP_sm3[i,j]+MVP_mm3[i,j])/(np.linalg.norm(miRNA3[i])+np.linalg.norm(lncRNA3[:, j]))


    # # #对关联矩阵的预处理：
    for i in range(m):
        for j in range(n):  # Adjusted loop to only iterate up to i
            MVP_association4[i,j] = (MVP_sm4[i,j]+MVP_mm4[i,j])/(np.linalg.norm(miRNA4[i])+np.linalg.norm(lncRNA4[:, j]))
    # # #对关联矩阵的预处理：
    for i in range(m):
        for j in range(n):  # Adjusted loop to only iterate up to i
            MVP_association5[i,j] = (MVP_sm5[i,j]+MVP_mm5[i,j])/(np.linalg.norm(miRNA5[i])+np.linalg.norm(lncRNA5[:, j]))


    MVP_association66 = np.multiply(MVP_association2 , MVP_association3)
    MVP_association77 = np.multiply(MVP_association4, MVP_association5)
    MVP_association = np.multiply(MVP_association66, MVP_association77)
    MVP_association = MVP_association/4

    #
    for i in range(m):
        for j in range(n):
            if(Y[i,j]==0):
                Y[i,j]=MVP_association[i,j]



    np.savetxt(r'MVPF_association.txt', Y, delimiter='\t', fmt='%.9f')

    return Y





if __name__ == "__main__":

    MVPF_association = main()
    print("预处理后的关联矩阵:")
    print(MVPF_association)




