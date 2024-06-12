import tkinter as tk
from tkinter import ttk, messagebox
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from math import factorial

x = sp.symbols('x')

def S_Taylor(f, x0, n):
    P = 0
    for k in range(n + 1):
        df = sp.diff(f, x, k)
        dfx0 = df.subs(x, x0)
        P = P + dfx0 * (x - x0)**k / factorial(k)
    return P

def Cota_t(f, x0, xp, n):
    m = min(x0, xp)
    M = max(x0, xp)
    u = np.linspace(m, M, 500)
    df = sp.diff(f, x, n + 1)
    df = sp.lambdify(x, df)
    Mc = np.max(np.abs(df(u)))
    return Mc * np.abs((xp - x0)**(n + 1) / factorial(n + 1))

class SeriesTaylorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Serie de Taylor App")
        
        self.create_widgets()

    def create_widgets(self):
        self.input_frame = ttk.LabelFrame(self.root, text="Datos de Entrada")
        self.input_frame.pack(padx=10, pady=10, fill="x")

        ttk.Label(self.input_frame, text="Función f(x):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.func_entry = ttk.Entry(self.input_frame, width=50)
        self.func_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Valor de x0:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.x0_entry = ttk.Entry(self.input_frame)
        self.x0_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Grado del Polinomio:").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.n_entry = ttk.Entry(self.input_frame)
        self.n_entry.grid(row=2, column=1, padx=5, pady=5)

        self.calculate_button = ttk.Button(self.input_frame, text="Calcular", command=self.calculate)
        self.calculate_button.grid(row=3, columnspan=2, pady=10)

        self.output_frame = ttk.LabelFrame(self.root, text="Resultados")
        self.output_frame.pack(padx=10, pady=10, fill="both", expand=True)

        self.poly_label = ttk.Label(self.output_frame, text="", wraplength=600)
        self.poly_label.pack(pady=10)

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.output_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def calculate(self):
        try:
            f_str = self.func_entry.get()
            x0 = float(self.x0_entry.get())
            n = int(self.n_entry.get())

            f = sp.sympify(f_str, locals={"ln": sp.ln})
            
            poly = S_Taylor(f, x0, n)
            self.poly_label.config(text=f"Polinomio de Taylor: {sp.pretty(poly)}")

            self.ax.clear()
            f_lambdified = sp.lambdify(x, f, modules=['numpy'])
            poly_lambdified = sp.lambdify(x, poly, modules=['numpy'])

            x_vals = np.linspace(x0 - 10, x0 + 10, 400)
            y_vals = f_lambdified(x_vals)
            y_poly_vals = poly_lambdified(x_vals)

            self.ax.plot(x_vals, y_vals, label='Función Original')
            self.ax.plot(x_vals, y_poly_vals, label=f'Polinomio de Taylor (grado {n})', linestyle='--')
            self.ax.legend()
            self.ax.set_title("Aproximación con Polinomio de Taylor")

            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", str(e))