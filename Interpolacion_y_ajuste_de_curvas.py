import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def Gauss_s(a, b, xo, tol):
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
        while (max(np.abs((x1 - xo))) > tol):
            xo = x1
            x1 = np.dot(Tg, xo) + Cg
            iteraciones += 1
        return x1
    else:
        print("El sistema iterativo no converge a la solucion unica del sistema")


def Pol_simple(x, y):
    n = len(x)
    M_p = np.zeros([n, n])
    x0 = np.zeros(n)
    for i in range(n):
        M_p[i, 0] = 1
        for j in range(1, n):
            M_p[i, j] = M_p[i, j - 1] * x[i]
    a_i = Gauss_s(M_p, y, x0, 1e-6)
    return a_i


def Lagrange(xd, yd):
    n = len(xd)
    P = 0

    for i in range(n):
        L = 1  # Polinomio base L_i(x)
        for j in range(n):
            if j != i:
                L *= (x - xd[j]) / (xd[i] - xd[j])  # Construcción de L_i(x)
        P += L * yd[i]  # Sumar el término L_i(x) * y_i al polinomio P
    Poly = sp.expand(P)
    return Poly


def min_c(x, y):
    Sx = sum(x)
    Sf = sum(y)
    Sx2 = sum((x ** 2))
    Sxy = sum((x * y))
    n = len(x)
    a0 = (Sf * Sx2 - Sx * Sxy) / (n * Sx2 - Sx ** 2)
    a1 = (n * Sxy - Sx * Sf) / (n * Sx2 - Sx ** 2)
    return a0, a1
