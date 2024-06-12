import time
import numpy as np
from scipy.linalg import eigvals

def Eliminacion_Gaussiana(A,b):
    n = len(b)
    x = np.zeros(n)
    for k in range(0,n-1):
        for i in range(k+1,n):
            lam=A[i,k]/(A[k,k])
            A[i,k:n]=A[i,k:n]-lam*A[k,k:n]
            b[i]=b[i]-lam*b[k]
    for k in range(n-1,-1,-1):
        x[k] = (b[k]-np.dot(A[k,k+1:n],x[k+1:n]))/(A[k,k])
    print("La solucion es: ",x)
    return x

def Gauss_seidel(a, b, xo, tol=1e-6):
    D = np.diag(np.diag(a))
    L = D - np.tril(a)
    U = D - np.triu(a)
    Tg = np.dot(np.linalg.inv(D - L), U)
    Cg = np.dot(np.linalg.inv(D - L), b)
    lam, vec = np.linalg.eig(Tg)
    radio = max(abs(lam))
    if radio < 1:
        x1 = np.dot(Tg, xo) + Cg
        iteraciones = 1
        while max(np.abs((x1 - xo))) > tol:
            xo = x1
            x1 = np.dot(Tg, xo) + Cg
            iteraciones += 1
        return x1
    else:
        print("El sistema iterativo no converge a la solucion unica del sistema")

def Gauss_sum(A,B,X0,tol):
    n=len(B)
    norm = 2
    cont = 0
    M = 50
    X1 = np.zeros(n)
    while norm>= tol or cont >M:
        for i in range(n):
            aux = 0
            for j in range(n):
                if i != j:
                    aux = aux-A[i,j]*X0[j]
            X1[i]=(B[i]+aux)/A[i,i]
            print(X1,X0)
            norm = np.max(np.abs(X1-X0))
            X0=X1
