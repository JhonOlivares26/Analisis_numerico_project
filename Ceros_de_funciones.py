import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from math import factorial
import math

x = sp.symbols('x')

def biseccion(f, a, b, tol):  # tolerancia
    if f(a) * f(b) > 0:
        print("no se cumple el teorema en ", [a, b])
    else:
        i = 0
        while np.abs(b - a) > tol:
            i +=1
            c = (a + b) / 2
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c
            print(c)
        print("iteraciones: ",i)
    return c

def newton(f, x0, tol):
    x = sp.symbols('x')
    df = sp.diff(f, x)
    NewT = x - f / df
    NewT = sp.lambdify(x, NewT)
    x1 = NewT(x0)
    i=0
    while abs(x1 - x0) > tol:
        i+=1
        x0 = x1
        x1 = NewT(x0)
        print(x1)
    print("iteraciones: ", i)
    return x1

def posicion_falsa(f, a, b, tol):
    c = 0
    if f(a) * f(b) > 0:
        print("No cumple el teorema")
    else:
        i=0
        c = a - f(a) * (a - b) / (f(a) - f(b))
        while abs(f(c)) > tol:
            i+=1
            c = a - f(a) * (a - b) / (f(a) - f(b))
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c
            # print(abs(b - a), abs(f(c)))
        print("iteraciones: ",i)
        print(f"La solución es {c}")
    return c

def Secante(f, x0, xi, tol):
    S = 0
    i = 0
    while abs(xi - x0)>tol:
        S = xi - f(xi)*(x0-xi)/(f(x0)-f(xi))
        x0 = xi
        xi = S
        i += 1
    print(f'La solución es: {S} con {i} iteraciones fueron:')
    return S
