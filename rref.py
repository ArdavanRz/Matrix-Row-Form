import numpy as np
from sympy import *
from scipy.linalg import null_space
from numpy.linalg import matrix_rank
def rref(mat,precision=0,GJ=False):
    m,n = mat.shape
    p,t = precision, 1e-1**precision
    A = around(mat.astype(float).copy(),decimals=p )
    if GJ:
        A = hstack((A,identity(n)))
    pcol = -1 #pivot colum
    for i in xrange(m):
        pcol += 1
        if pcol >= n : break
        #pivot index
        pid = argmax( abs(A[i:,pcol]) )
        #Row exchange
        A[i,:],A[pid+i,:] = A[pid+i,:].copy(),A[i,:].copy()
        #pivot with given precision
        while pcol < n and abs(A[i,pcol]) < t:
            #pivot index
            pid = argmax( abs(A[i:,pcol]) )
            #Row exchange
            A[i,:],A[pid+i,:] = A[pid+i,:].copy(),A[i,:].copy()
            pcol += 1
        if pcol >= n : break
        pivot = float(A[i,pcol])
        for j in xrange(m):
            if j == i: continue
            mul = float(A[j,pcol])/pivot
            A[j,:] = around(A[j,:] - A[i,:]*mul,decimals=p)
        A[i,:] /= pivot
        A[i,:] = around(A[i,:],decimals=p)

    if GJ:
        return A[:,:n].copy(),A[:,n:].copy()
    else:
        return A
#############################################################################################
run = true
c = 0
while run == true :
    print("enter 1 for rref matrix solver or 2 for null_space and rank finder")
    c = input()
    choice = int(c)
    if choice == 1:
        print("enter the total number of variables present in one line including the right hand")
        m = input()
        columns = int(m)
        print("enter the total number of equations")
        n = input()
        rows = int(n)
        print("then, the size of your augmented matrix is", rows, "by", columns, )
        c_counter = 0
        r_counter = 0
        mx = []
        ##############################################################################################
        while r_counter < rows:
            print("enter varabile multipliers and the right hand result for row number", r_counter + 1, "\n")
            while c_counter < columns:
                if c_counter + 1 < columns:
                    print("enter multiplier for X", c_counter + 1, "if not present enter zero")
                    x = input()
                    mx.append(int(x))
                    c_counter = c_counter + 1
                if c_counter + 1 == columns:
                    print("enter right hand result")
                    x = input()
                    mx.append(int(x))
                    c_counter = 0
                    break
            r_counter = r_counter + 1
        #############################################################################################
        np.asarray(mx)
        print("your augmented matrix would look like this\n")
        print(np.resize(mx, (rows, columns)))
        #############################################################################################
        newmx = Matrix(np.resize(mx, (rows, columns)))
        mx_rref = newmx.rref()
        print("The Row echelon form of matrix M and the pivot columns : {}".format(mx_rref))

    if choice == 2 :
            print("how many columns in your matrix?")
            m = input()
            columns = int(m)
            print("how many rows in your matrix?")
            n = input()
            rows = int(n)
            c_counter = 0
            r_counter = 0
            mx = []
            ##############################################################################################
            while r_counter < rows:
                print("enter the numbers for row number", r_counter + 1, "\n")
                while c_counter < columns:
                    if c_counter + 1 < columns:
                        print("enter data", c_counter + 1, "if not present enter zero")
                        x = input()
                        mx.append(int(x))
                        c_counter = c_counter + 1
                    if c_counter + 1 == columns:
                        print("enter data", c_counter + 1, "if not present enter zero")
                        x = input()
                        mx.append(int(x))
                        c_counter = 0
                        break
                r_counter = r_counter + 1
            np.asarray(mx)
            ns = null_space(np.resize(mx, (rows, columns)))
            print("the null_space of your matrix is :",ns,)
            newmx = np.matrix(np.resize(mx, (rows, columns)))
            print("the rank of your matrix is :",np.linalg.matrix_rank(newmx),)









