import sympy as sp
import numpy as np
import tkinter as tk
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

class CerosDeFuncionesApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Ceros de Funciones")
        root.geometry("600x300")

        # Crear campos de entrada para los parámetros
        self.f_entry = tk.Entry(root)
        self.a_entry = tk.Entry(root)
        self.b_entry = tk.Entry(root)
        self.tol_entry = tk.Entry(root)
        self.x0_entry = tk.Entry(root)
        self.xi_entry = tk.Entry(root)

        # Crear botones para ejecutar las funciones
        self.biseccion_button = tk.Button(root, text="Bisección", command=self.run_biseccion)
        self.newton_button = tk.Button(root, text="Newton", command=self.run_newton)
        self.posicion_falsa_button = tk.Button(root, text="Posición Falsa", command=self.run_posicion_falsa)
        self.secante_button = tk.Button(root, text="Secante", command=self.run_secante)

        # Posicionar los campos de entrada y los botones
        self.f_entry.pack()
        self.a_entry.pack()
        self.b_entry.pack()
        self.tol_entry.pack()
        self.x0_entry.pack()
        self.xi_entry.pack()
        self.biseccion_button.pack()
        self.newton_button.pack()
        self.posicion_falsa_button.pack()
        self.secante_button.pack()

    def run_biseccion(self):
        # Obtener los valores de los campos de entrada
        f = self.f_entry.get()
        a = float(self.a_entry.get())
        b = float(self.b_entry.get())
        tol = float(self.tol_entry.get())

        # Ejecutar la función biseccion y mostrar los resultados
        c = biseccion(f, a, b, tol)
        print(c)

    def run_newton(self):
        # Obtener los valores de los campos de entrada
        f = self.f_entry.get()
        x0 = float(self.x0_entry.get())
        tol = float(self.tol_entry.get())

        # Ejecutar la función newton y mostrar los resultados
        x1 = newton(f, x0, tol)
        print(x1)

    def run_posicion_falsa(self):
        # Obtener los valores de los campos de entrada
        f = self.f_entry.get()
        a = float(self.a_entry.get())
        b = float(self.b_entry.get())
        tol = float(self.tol_entry.get())

        # Ejecutar la función posicion_falsa y mostrar los resultados
        c = posicion_falsa(f, a, b, tol)
        print(c)

    def run_secante(self):
        # Obtener los valores de los campos de entrada
        f = self.f_entry.get()
        x0 = float(self.x0_entry.get())
        xi = float(self.xi_entry.get())
        tol = float(self.tol_entry.get())

        # Ejecutar la función Secante y mostrar los resultados
        S = Secante(f, x0, xi, tol)
        print(S)
