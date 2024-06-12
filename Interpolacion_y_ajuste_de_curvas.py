import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox
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
    return np.flip(a_i)

def Lagrange(xd, yd):
    n = len(xd)
    P = 0
    x = sp.symbols('x')  # Definir la variable simbólica 'x'

    for i in range(n):
        L = 1
        for j in range(n):
            if j != i:
                L *= (x - xd[j]) / (xd[i] - xd[j])
        P += L * yd[i]
    Poly = sp.expand(P)
    return Poly

def min_c(x, y):
    Sx = sum(x)
    Sf = sum(y)
    Sx2 = sum((np.array(x) ** 2))
    Sxy = sum((np.array(x) * np.array(y)))
    n = len(x)
    a1 = (n * Sxy - Sx * Sf) / (n * Sx2 - Sx ** 2)
    a0 = (Sf * Sx2 - Sx * Sxy) / (n * Sx2 - Sx ** 2)
    return a1, a0

class InterpolationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Interpolación y Ajuste de Curvas")
        self.root.geometry("1400x600")

        # Crear un marco para los controles
        control_frame = ttk.LabelFrame(root, text="Datos de Entrada")
        control_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=10, pady=10)

        # Crear un marco para la gráfica
        self.graph_frame = ttk.LabelFrame(root, text="Gráfica")
        self.graph_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Crear etiquetas y campos de entrada para los parámetros
        ttk.Label(control_frame, text="Datos x (separados por comas):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.x_entry = ttk.Entry(control_frame)
        self.x_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Datos y (separados por comas):").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.y_entry = ttk.Entry(control_frame)
        self.y_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Dato a aproximar (opcional):").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.x_approx_entry = ttk.Entry(control_frame)
        self.x_approx_entry.grid(row=2, column=1, padx=5, pady=5)

        # Crear botones para ejecutar las funciones
        self.pol_simple_button = ttk.Button(control_frame, text="Polinomio Simple", command=self.run_pol_simple)
        self.pol_simple_button.grid(row=3, column=0, padx=5, pady=10, sticky="ew")

        self.lagrange_button = ttk.Button(control_frame, text="Lagrange", command=self.run_lagrange)
        self.lagrange_button.grid(row=3, column=1, padx=5, pady=10, sticky="ew")

        self.min_cuad_button = ttk.Button(control_frame, text="Mínimos Cuadrados", command=self.run_min_cuad)
        self.min_cuad_button.grid(row=4, column=0, columnspan=2, padx=5, pady=10, sticky="ew")

        instructions = (
            "Instrucciones de uso:\n"
            "1. Ingrese los datos x e y separados por comas en los campos correspondientes.\n"
            "2. Si desea aproximar un valor, ingrese el valor en el campo 'Dato a aproximar'.\n"
            "3. Seleccione el método de interpolación o ajuste de curvas que desea utilizar.\n"
            "4. Haga clic en el botón correspondiente al método seleccionado.\n"
            "5. Los resultados se mostrarán en el área de texto a la derecha y la gráfica se actualizará automáticamente."
        )
        self.instructions_label = ttk.Label(control_frame, text=instructions)
        self.instructions_label.grid(row=5, column=0, columnspan=2, padx=5, pady=10)

        # Crear un widget Text para mostrar los resultados
        self.result_text = tk.Text(self.graph_frame, height=10, width=60)
        self.result_text.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Inicializar la variable del lienzo de la gráfica
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def plot_solution(self, x_vals, y_vals, poly, method_name, x_approx=None):
        # Limpiar la gráfica
        self.ax.clear()

        # Graficar los datos originales
        self.ax.plot(x_vals, y_vals, 'bo', label='Datos')

        # Graficar el polinomio de ajuste o interpolación
        x_min = min(x_vals)
        x_max = max(x_vals)
        x_plot = np.linspace(x_min, x_max, 100)

        if isinstance(poly, np.poly1d):
            y_plot = poly(x_plot)
            label = f'Polinomio: {poly}'
        else:
            y_plot = [poly.evalf(subs={'x': xi}) for xi in x_plot]
            label = str(poly)

        self.ax.plot(x_plot, y_plot, label=label)

        if x_approx:
            x_approx = float(x_approx)
            y_approx = poly(x_approx) if isinstance(poly, np.poly1d) else poly.evalf(subs={'x': x_approx})
            self.ax.plot(x_approx, y_approx, 'ro', label=f'Aproximado: ({x_approx}, {y_approx:.2f})')

        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_title(f'Interpolación y Ajuste de Curvas: {method_name}')
        self.ax.legend()

        self.canvas.draw()

        self.result_text.insert(tk.END, f"{label}\n")
        if x_approx:
            self.result_text.insert(tk.END, f"Resultado aproximado en x={x_approx}: {y_approx}\n")

    def run_pol_simple(self):
        try:
            x_vals = list(map(float, self.x_entry.get().split(',')))
            y_vals = list(map(float, self.y_entry.get().split(',')))
            poly_coeffs = Pol_simple(x_vals, y_vals)
            poly = np.poly1d(poly_coeffs)
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Polinomio Simple: {poly}\n")
            self.plot_solution(x_vals, y_vals, poly, "Polinomio Simple", self.x_approx_entry.get() if self.x_approx_entry.get() else None)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def run_lagrange(self):
        try:
            x_vals = list(map(float, self.x_entry.get().split(',')))
            y_vals = list(map(float, self.y_entry.get().split(',')))
            poly = Lagrange(x_vals, y_vals)
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Polinomio de Lagrange: {poly}\n")
            self.plot_solution(x_vals, y_vals, poly, "Lagrange", self.x_approx_entry.get() if self.x_approx_entry.get() else None)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def run_min_cuad(self):
        try:
            x_vals = list(map(float, self.x_entry.get().split(',')))
            y_vals = list(map(float, self.y_entry.get().split(',')))
            a1, a0 = min_c(x_vals, y_vals)
            poly = np.poly1d([a1, a0])
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Polinomio de Mínimos Cuadrados: {poly}\n")
            self.plot_solution(x_vals, y_vals, poly, "Mínimos Cuadrados", self.x_approx_entry.get() if self.x_approx_entry.get() else None)
        except Exception as e:
            messagebox.showerror("Error", str(e))