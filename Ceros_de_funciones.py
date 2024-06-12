import tkinter as tk
from tkinter import ttk, messagebox
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
        
        self.create_widgets()

    def create_widgets(self):
        self.input_frame = ttk.LabelFrame(self.root, text="Informacion de uso")
        self.input_frame.pack(padx=10, pady=10, fill="x")

        ttk.Label(self.input_frame, text="La funcion se debe ingresar con los datos completa y las incognitas en terminos de X (ejemplo: 2*x + 4*x), \n para ingresar exponenciales, logaritmos,etc, se recomienda Sympy con el prefijo sp(ej sp.exp(x))").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(self.input_frame, text="El programa para calcular con newton,para el valor X0 se calcula el promedio que se ingrese en los campos a y b,\n para secante, el valor X0 y X1 son los valores ingresados en a y b respectivamente").grid(row=2, column=0, padx=5, pady=5, sticky="e") 

        self.input_frame = ttk.LabelFrame(self.root, text="Datos de Entrada")
        self.input_frame.pack(padx=10, pady=10, fill="x")

        ttk.Label(self.input_frame, text="Función (en términos de x):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.f_entry = ttk.Entry(self.input_frame, width=50)
        self.f_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Valor de a:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.a_entry = ttk.Entry(self.input_frame)
        self.a_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Valor de b:").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.b_entry = ttk.Entry(self.input_frame)
        self.b_entry.grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Tolerancia:").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.tol_entry = ttk.Entry(self.input_frame)
        self.tol_entry.grid(row=3, column=1, padx=5, pady=5)

        self.button_frame = ttk.Frame(self.root)
        self.button_frame.pack(padx=10, pady=10, fill="x")

        self.biseccion_button = ttk.Button(self.button_frame, text="Bisección", command=self.run_biseccion)
        self.biseccion_button.grid(row=0, column=0, padx=5, pady=5)

        self.newton_button = ttk.Button(self.button_frame, text="Newton", command=self.run_newton)
        self.newton_button.grid(row=0, column=1, padx=5, pady=5)

        self.posicion_falsa_button = ttk.Button(self.button_frame, text="Posición Falsa", command=self.run_posicion_falsa)
        self.posicion_falsa_button.grid(row=0, column=2, padx=5, pady=5)

        self.secante_button = ttk.Button(self.button_frame, text="Secante", command=self.run_secante)
        self.secante_button.grid(row=0, column=3, padx=5, pady=5)

        self.output_frame = ttk.LabelFrame(self.root, text="Resultados")
        self.output_frame.pack(padx=10, pady=10, fill="both", expand=True)

        self.result_text = tk.Text(self.output_frame, height=10, width=70)
        self.result_text.pack(padx=5, pady=5, fill="both", expand=True)

    def run_biseccion(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f = self.get_function()
            a = float(self.a_entry.get())
            b = float(self.b_entry.get())
            tol = float(self.tol_entry.get())
            result = biseccion(f, a, b, tol)
        except Exception as e:
            result = f"Error: {str(e)}"
        self.result_text.insert(tk.END, result)

    def run_newton(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f_expr = self.get_function_expr()
            x0 = (float(self.a_entry.get()) + float(self.b_entry.get())) / 2
            tol = float(self.tol_entry.get())
            result = newton(f_expr, x0, tol)
        except Exception as e:
            result = f"Error: {str(e)}"
        self.result_text.insert(tk.END, result)

    def run_posicion_falsa(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f = self.get_function()
            a = float(self.a_entry.get())
            b = float(self.b_entry.get())
            tol = float(self.tol_entry.get())
            result = posicion_falsa(f, a, b, tol)
        except Exception as e:
            result = f"Error: {str(e)}"
        self.result_text.insert(tk.END, result)

    def run_secante(self):
        self.result_text.delete(1.0, tk.END)
        try:
            f = self.get_function()
            x0 = float(self.a_entry.get())
            x1 = float(self.b_entry.get())
            tol = float(self.tol_entry.get())
            result = secante(f, x0, x1, tol)
        except Exception as e:
            result = f"Error: {str(e)}"
        self.result_text.insert(tk.END, result)

    def get_function(self):
        f_str = self.f_entry.get()
        f = eval(f"lambda x: {f_str}", {"x": x, "sp": sp, "np": np, "math": math})
        return f

    def get_function_expr(self):
        f_str = self.f_entry.get()
        f_expr = sp.sympify(f_str, {"x": x, "sp": sp, "np": np, "math": math})
        return f_expr