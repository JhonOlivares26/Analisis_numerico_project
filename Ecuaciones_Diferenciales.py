import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
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
        root.geometry("1200x600")  # Ajustar el tamaño de la ventana

        # Crear un marco para los controles
        control_frame = tk.Frame(root, width=300)
        control_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False)

        # Crear un marco para la gráfica
        self.graph_frame = tk.Frame(root, width=700)
        self.graph_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Crear etiquetas y campos de entrada para los parámetros
        tk.Label(control_frame, text="Función (en términos de t y y):").pack()
        self.f_entry = tk.Entry(control_frame)
        self.f_entry.pack()

        tk.Label(control_frame, text="Valor inicial a:").pack()
        self.a_entry = tk.Entry(control_frame)
        self.a_entry.pack()

        tk.Label(control_frame, text="Valor final b:").pack()
        self.b_entry = tk.Entry(control_frame)
        self.b_entry.pack()

        tk.Label(control_frame, text="Tamaño de paso h:").pack()
        self.h_entry = tk.Entry(control_frame)
        self.h_entry.pack()

        tk.Label(control_frame, text="Valor inicial y0 (separados por comas si son múltiples):").pack()
        self.y0_entry = tk.Entry(control_frame)
        self.y0_entry.pack()

        # Crear botones para ejecutar las funciones
        self.runge_kutta_button = tk.Button(control_frame, text="Runge Kutta", command=self.run_runge_kutta)
        self.runge_kutta_button.pack()

        self.euler_button = tk.Button(control_frame, text="Euler", command=self.run_euler)
        self.euler_button.pack()

        # Crear un widget Text para mostrar los resultados
        self.result_text = tk.Text(control_frame, height=10)
        self.result_text.pack()

        # Inicializar la variable del lienzo de la gráfica
        self.canvas = None

    def plot_solution(self, t_vals, y_vals, method_name):
        # Limpiar el marco de la gráfica
        for widget in self.graph_frame.winfo_children():
            widget.destroy()

        # Crear una nueva figura de matplotlib
        fig, ax = plt.subplots()
        ax.plot(t_vals, y_vals, label=method_name)
        ax.set_xlabel('t')
        ax.set_ylabel('y')
        ax.set_title(f'Solución usando {method_name}')
        ax.legend()

        # Crear un lienzo de tkinter para la figura de matplotlib
        self.canvas = FigureCanvasTkAgg(fig, master=self.graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def run_runge_kutta(self):
        # Limpiar el contenido del widget de texto
        self.result_text.delete(1.0, tk.END)

        # Obtener los valores de los campos de entrada
        f_str = self.f_entry.get()
        a = float(self.a_entry.get())
        b = float(self.b_entry.get())
        h = float(self.h_entry.get())
        y0 = list(map(float, self.y0_entry.get().split(',')))

        # Convertir la cadena de texto en una función
        t, y = sp.symbols('t y')
        f_expr = sp.sympify(f_str)
        f_lambda = sp.lambdify((t, y), f_expr, 'numpy')

        # Ejecutar la función RungeKutta y mostrar los resultados
        t_vals, y_vals = RungeKutta(f_lambda, a, b, h, y0)
        self.result_text.insert(tk.END, f"Resultados de Runge Kutta:\nt: {t_vals}\ny: {y_vals}\n")
        self.plot_solution(t_vals, y_vals, "Runge Kutta")

    def run_euler(self):
        # Limpiar el contenido del widget de texto
        self.result_text.delete(1.0, tk.END)

        # Obtener los valores de los campos de entrada
        f_str = self.f_entry.get()
        a = float(self.a_entry.get())
        b = float(self.b_entry.get())
        h = float(self.h_entry.get())
        y0 = list(map(float, self.y0_entry.get().split(',')))

        # Convertir la cadena de texto en una función
        t, y = sp.symbols('t y')
        f_expr = sp.sympify(f_str)
        f_lambda = sp.lambdify((t, y), f_expr, 'numpy')

        # Ejecutar la función Euler y mostrar los resultados
        t_vals, y_vals = Euler(f_lambda, a, b, h, y0)
        self.result_text.insert(tk.END, f"Resultados de Euler:\nt: {t_vals}\ny: {y_vals}\n")
        self.plot_solution(t_vals, y_vals, "Euler")








