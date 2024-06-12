import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from math import factorial
import math
import pandas as pd
x = sp.symbols('x')

def S_Taylor(f,x0,n):
    P=0
    for k in range(n+1):
        df = sp.diff(f,x,k)
        dfx0 = df.subs(x,x0)
        P = P+dfx0*(x-x0)**k/factorial(k)
    return P

def Cota_t(f,x0,xp,n):
    m = min(x0,xp)
    M = max(x0,xp)
    u = np.linspace(m,M,500)
    df = sp.diff(f,x,n+1)
    df = sp.lambdify(x,df)
    Mc = np.max(np.abs(df(u)))
    return Mc*np.abs((xp-x0)**(n+1)/factorial(n+1))