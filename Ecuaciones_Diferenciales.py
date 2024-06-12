import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox
import sympy as sp

def RungeKutta(f, a, b, h, y0):
    n = int((b - a) / h)
    t = np.linspace(a, b, n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        k1 = h * np.array(f(t[i], y[i])).astype(float)
        k2 = h * np.array(f(t[i] + h / 2, y[i] + k1 / 2)).astype(float)
        k3 = h * np.array(f(t[i] + h / 2, y[i] + k2 / 2)).astype(float)
        k4 = h * np.array(f(t[i] + h, y[i] + k3)).astype(float)
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return t, y

def Euler(f, a, b, h, y0):
    n = int((b - a) / h)
    t = np.linspace(a, b, n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        y[i + 1] = y[i] + h * np.array(f(t[i], y[i])).astype(float)
    return t, y

class DifferentialEquationsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Ecuaciones Diferenciales")
        self.root.geometry("1000x600")  # Ajustar el tamaño de la ventana

        # Crear un marco para los controles
        control_frame = ttk.LabelFrame(self.root, text="Datos de Entrada")
        control_frame.pack(padx=10, pady=10, fill="both", expand=True)

        # Crear un marco para la gráfica
        self.graph_frame = ttk.LabelFrame(self.root, text="Gráfica")
        self.graph_frame.pack(padx=10, pady=10, fill="both", expand=True)

        # Crear etiquetas y campos de entrada para los parámetros
        ttk.Label(control_frame, text="Función f(t, y):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.f_entry = ttk.Entry(control_frame, width=50)
        self.f_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Valor inicial t (a):").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.t0_entry = ttk.Entry(control_frame)
        self.t0_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Valor final t (b):").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.tf_entry = ttk.Entry(control_frame)
        self.tf_entry.grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Tamaño de paso h:").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.h_entry = ttk.Entry(control_frame)
        self.h_entry.grid(row=3, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Valor inicial y0 (separados por comas si son múltiples):").grid(row=4, column=0, padx=5, pady=5, sticky="e")
        self.y0_entry = ttk.Entry(control_frame)
        self.y0_entry.grid(row=4, column=1, padx=5, pady=5)

        # Crear botones para ejecutar las funciones
        self.runge_kutta_button = ttk.Button(control_frame, text="Runge Kutta", command=self.run_runge_kutta)
        self.runge_kutta_button.grid(row=5, column=0, padx=5, pady=10, sticky="ew")

        self.euler_button = ttk.Button(control_frame, text="Euler", command=self.run_euler)
        self.euler_button.grid(row=5, column=1, padx=5, pady=10, sticky="ew")

        # Crear un widget Text para mostrar los resultados
        self.result_text = tk.Text(self.graph_frame, height=10, width=60)
        self.result_text.pack(side=tk.LEFT, fill="both", expand=True)

        # Inicializar la variable del lienzo de la gráfica
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill="both", expand=True)

    def plot_solution(self, t_vals, y_vals, method_name):
        # Limpiar la gráfica
        self.ax.clear()

        # Graficar la solución
        self.ax.plot(t_vals, y_vals, label=method_name)
        self.ax.set_xlabel('t')
        self.ax.set_ylabel('y')
        self.ax.set_title(f'Solución usando {method_name}')
        self.ax.legend()

        # Actualizar el lienzo de tkinter
        self.canvas.draw()

    def run_runge_kutta(self):
        try:
            # Obtener los valores de los campos de entrada
            f_str = self.f_entry.get()
            t0 = float(self.t0_entry.get())
            tf = float(self.tf_entry.get())
            h = float(self.h_entry.get())
            y0 = list(map(float, self.y0_entry.get().split(',')))

            # Convertir la cadena de texto en una función
            t, y = sp.symbols('t y')  # Define los símbolos
            f_expr = sp.sympify(f_str)  # Convierte la cadena en una expresión sympy
            f_lambda = sp.lambdify((t, y), f_expr, modules=['numpy'])  # Convierte la expresión en una función lambda

            # Ejecutar la función RungeKutta y mostrar los resultados
            t_vals, y_vals = RungeKutta(f_lambda, t0, tf, h, y0)
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"t: {t_vals}\n\ny: {y_vals}\n")
            self.plot_solution(t_vals, y_vals, "Runge Kutta")

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def run_euler(self):
        try:
            # Obtener los valores de los campos de entrada
            f_str = self.f_entry.get()
            t0 = float(self.t0_entry.get())
            tf = float(self.tf_entry.get())
            h = float(self.h_entry.get())
            y0 = list(map(float, self.y0_entry.get().split(',')))

            # Convertir la cadena de texto en una función
            t, y = sp.symbols('t y')
            f_expr = sp.sympify(f_str)
            f_lambda = sp.lambdify((t, y), f_expr, 'numpy')

            # Ejecutar la función Euler y mostrar los resultados
            t_vals, y_vals = Euler(f_lambda, t0, tf, h, y0)
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Resultados de Euler:\nt: {t_vals}\ny: {y_vals}\n")
            self.plot_solution(t_vals, y_vals, "Euler")

        except Exception as e:
            messagebox.showerror("Error", str(e))



