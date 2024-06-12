import tkinter as tk
import sympy as sp
import numpy as np
import math

x = sp.symbols('x')

def biseccion(f, a, b, tol):
    if f(a) * f(b) > 0:
        return "No se cumple el teorema en [{}, {}]".format(a, b)
    else:
        i = 0
        while np.abs(b - a) > tol:
            i += 1
            c = (a + b) / 2
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c
        return "Resultado con biseccion: {} en {} iteraciones".format(c, i)

def newton(f_expr, x0, tol):
    df_expr = sp.diff(f_expr, x)
    f = sp.lambdify(x, f_expr, "numpy")
    df = sp.lambdify(x, df_expr, "numpy")
    i = 0
    while True:
        try:
            x1 = x0 - f(x0) / df(x0)
        except ZeroDivisionError:
            return "Error: división por cero durante la evaluación en x = {}".format(x0)
        if abs(x1 - x0) < tol:
            break
        x0 = x1
        i += 1
    return "Resultado con newton: {} en {} iteraciones".format(x1, i)

def posicion_falsa(f, a, b, tol):
    if f(a) * f(b) > 0:
        return "No se cumple el teorema en [{}, {}]".format(a, b)
    else:
        i = 0
        while True:
            c = a - f(a) * (b - a) / (f(b) - f(a))
            if abs(f(c)) < tol:
                break
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c
            i += 1
        return "Resultado con posicion falsa: {} en {} iteraciones".format(c, i)

def secante(f, x0, x1, tol):
    i = 0
    while True:
        try:
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        except ZeroDivisionError:
            return "Error: división por cero durante la evaluación en x = {}".format(x1)
        if abs(x2 - x1) < tol:
            break
        x0, x1 = x1, x2
        i += 1
    return "Resultado con secante: {} en {} iteraciones".format(x2, i)

class CerosDeFuncionesApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Ceros de Funciones")
        root.geometry("600x600")

        self.f_entry = tk.Entry(root, width=50)
        self.a_entry = tk.Entry(root)
        self.b_entry = tk.Entry(root)
        self.tol_entry = tk.Entry(root)
        self.x0_entry = tk.Entry(root)
        self.xi_entry = tk.Entry(root)

        tk.Label(root, text="Función (en términos de x):").pack()
        self.f_entry.pack()
        tk.Label(root, text="a (para Biseccion y Falsa posicion):").pack()
        self.a_entry.pack()
        tk.Label(root, text="b (para Biseccion y Falsa posicion):").pack()
        self.b_entry.pack()
        tk.Label(root, text="Tolerancia:").pack()
        self.tol_entry.pack()
        tk.Label(root, text="x0 (para Newton y Secante):").pack()
        self.x0_entry.pack()
        tk.Label(root, text="xi (para Secante):").pack()
        self.xi_entry.pack()

        self.biseccion_button = tk.Button(root, text="Bisección", command=self.run_biseccion)
        self.newton_button = tk.Button(root, text="Newton", command=self.run_newton)
        self.posicion_falsa_button = tk.Button(root, text="Posición Falsa", command=self.run_posicion_falsa)
        self.secante_button = tk.Button(root, text="Secante", command=self.run_secante)

        self.biseccion_button.pack(pady=5)
        self.newton_button.pack(pady=5)
        self.posicion_falsa_button.pack(pady=5)
        self.secante_button.pack(pady=5)

        self.result_text = tk.Text(root, height=10, width=70)
        self.result_text.pack(pady=10)

    def run_biseccion(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f = self.get_function()
            a = float(self.a_entry.get())
            b = float(self.b_entry.get())
            tol = float(self.tol_entry.get())
            result = biseccion(f, a, b, tol)
        except ValueError as e:
            result = f"Error de valor: {e}"
        self.result_text.insert(tk.END, result)

    def run_newton(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f_expr = self.get_function_expr()
            x0 = float(self.x0_entry.get())
            tol = float(self.tol_entry.get())
            result = newton(f_expr, x0, tol)
        except ValueError as e:
            result = f"Error de valor: {e}"
        self.result_text.insert(tk.END, result)

    def run_posicion_falsa(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f = self.get_function()
            a = float(self.a_entry.get())
            b = float(self.b_entry.get())
            tol = float(self.tol_entry.get())
            result = posicion_falsa(f, a, b, tol)
        except ValueError as e:
            result = f"Error de valor: {e}"
        self.result_text.insert(tk.END, result)

    def run_secante(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f = self.get_function()
            x0 = float(self.x0_entry.get())
            xi = float(self.xi_entry.get())
            tol = float(self.tol_entry.get())
            result = secante(f, x0, xi, tol)
        except ValueError as e:
            result = f"Error de valor: {e}"
        self.result_text.insert(tk.END, result)

    def get_function(self):
        f_str = self.f_entry.get()
        f = eval(f"lambda x: {f_str}", {"x": x, "sp": sp, "np": np, "math": math})
        return f

    def get_function_expr(self):
        f_str = self.f_entry.get()
        f_expr = sp.sympify(f_str, {"x": x, "sp": sp, "np": np, "math": math})
        return f_expr