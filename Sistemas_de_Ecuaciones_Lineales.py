import time
import numpy as np
from scipy.linalg import eigvals

def Gauss_seidel(a,b,xo,tol=1e-6):
  D=np.diag(np.diag(a))
  L=D-np.tril(a)
  U=D-np.triu(a)
  Tg=np.dot(np.linalg.inv(D-L),U)
  Cg=np.dot(np.linalg.inv(D-L),b)
  lam,vec=np.linalg.eig(Tg)
  radio=max(abs(lam))
  if radio<1:
    x1=np.dot(Tg,xo)+Cg
    iteraciones = 1
    while max(np.abs((x1 - xo)))>tol:
        xo=x1
        x1=np.dot(Tg,xo)+Cg
        iteraciones+=1
    return x1
  else:
    print("El sistema iterativo no converge a la solucion unica del sistema")


def Gauss_seidel_time(a,b,xo,tol=1e-6):
  D=np.diag(np.diag(a))
  L=D-np.tril(a)
  U=D-np.triu(a)
  Tg=np.dot(np.linalg.inv(D-L),U)
  Cg=np.dot(np.linalg.inv(D-L),b)
  lam,vec=np.linalg.eig(Tg)
  radio=max(abs(lam))
  if radio<1:
    x1=np.dot(Tg,xo)+Cg
    iteraciones = 1
    start_time = time.time()
    while max(np.abs((x1 - xo)))>tol:
        xo=x1
        x1=np.dot(Tg,xo)+Cg
        iteraciones+=1
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"El tiempo de computo fue de {execution_time}")
    return x1
  else:
    print("El sistema iterativo no converge a la solucion unica del sistema")

def Gauss_sum_time(A, B, X0, tol=1e-6):
    n = len(B)
    norm = tol + 1
    cont = 0
    max_iter = 1000
    X1 = np.zeros(n)
    start_time = time.time()
    while norm >= tol and cont < max_iter:
        for i in range(n):
            aux = 0
            for j in range(n):
                if i != j:
                    aux = aux - A[i, j] * X0[j]
            X1[i] = (B[i] + aux) / A[i, i]
        norm = np.max(np.abs(X1 - X0))
        X0 = X1.copy()
        cont += 1
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"El tiempo de computo fue de {execution_time}")
    if cont == max_iter:
        print("El método no convergió en el número máximo de iteraciones.")
    return f"La solución del sistema es: {X1}, con {cont} iteracciones"

def Gauss_sum_errores(A, B, X0, tol=1e-6):
    n = len(B)
    norm = tol + 1
    max_iter = 1000
    X1 = np.zeros(n)
    errores_sumas = []
    iteracciones_sumas = 1
    while norm >= tol and iteracciones_sumas < max_iter:
        for i in range(n):
            aux = 0
            for j in range(n):
                if i != j:
                    aux = aux - A[i, j] * X0[j]
            X1[i] = (B[i] + aux) / A[i, i]
        norm = np.max(np.abs(X1 - X0))
        errores_sumas.append(norm)
        X0 = X1.copy()
        iteracciones_sumas += 1
    if iteracciones_sumas == max_iter:
        print("El método no convergió en el número máximo de iteraciones.")
    return X1, iteracciones_sumas, errores_sumas

def Gauss_seidel_errores(a,b,x0,tol=1e-6):
  D=np.diag(np.diag(a))
  L=D-np.tril(a)
  U=D-np.triu(a)
  Tg=np.dot(np.linalg.inv(D-L),U)
  Cg=np.dot(np.linalg.inv(D-L),b)
  lam,vec=np.linalg.eig(Tg)
  radio=max(abs(lam))
  if radio<1:
    x1=np.dot(Tg,x0)+Cg
    errores_seidel = [np.max(np.abs(x1 - x0))]
    iteracciones_seidel = 1
    while max(np.abs((x1 - x0)))>tol:
        x0=x1
        x1=np.dot(Tg,x0)+Cg
        errores_seidel.append(np.max(np.abs(x1 - x0)))
        iteracciones_seidel+=1
    return x1, iteracciones_seidel, errores_seidel
  else:
    print("El sistema iterativo no converge a la solucion unica del sistema")

