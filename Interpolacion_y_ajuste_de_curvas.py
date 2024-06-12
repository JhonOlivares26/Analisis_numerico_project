import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
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
    a0 = (Sf * Sx2 - Sx * Sxy) / (n * Sx2 - Sx ** 2)
    a1 = (n * Sxy - Sx * Sf) / (n * Sx2 - Sx ** 2)
    return a0, a1

class InterpolationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Interpolación y Ajuste de Curvas")
        root.geometry("1200x600")

        # Crear un marco para los controles
        control_frame = tk.Frame(root, width=300)
        control_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False)

        # Crear un marco para la gráfica
        self.graph_frame = tk.Frame(root, width=700)
        self.graph_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Crear etiquetas y campos de entrada para los parámetros
        tk.Label(control_frame, text="Datos x (separados por comas):").pack()
        self.x_entry = tk.Entry(control_frame)
        self.x_entry.pack()

        tk.Label(control_frame, text="Datos y (separados por comas):").pack()
        self.y_entry = tk.Entry(control_frame)
        self.y_entry.pack()

        tk.Label(control_frame, text="Dato a aproximar:").pack()
        self.x_approx_entry = tk.Entry(control_frame)
        self.x_approx_entry.pack()

        # Crear botones para ejecutar las funciones
        self.pol_simple_button = tk.Button(control_frame, text="Polinomio Simple", command=self.run_pol_simple)
        self.pol_simple_button.pack()

        self.lagrange_button = tk.Button(control_frame, text="Lagrange", command=self.run_lagrange)
        self.lagrange_button.pack()

        self.min_cuad_button = tk.Button(control_frame, text="Mínimos Cuadrados", command=self.run_min_cuad)
        self.min_cuad_button.pack()

        # Crear un widget Text para mostrar los resultados
        self.result_text = tk.Text(control_frame, height=10)
        self.result_text.pack()

        # Inicializar la variable del lienzo de la gráfica
        self.canvas = None

    def plot_solution(self, x_vals, y_vals, poly, x_approx):
        # Limpiar el marco de la gráfica
        for widget in self.graph_frame.winfo_children():
            widget.destroy()

        # Crear una nueva figura de matplotlib
        fig, ax = plt.subplots()
        ax.plot(x_vals, y_vals, 'bo', label='Datos')

        # Plotear el polinomio de ajuste o interpolación
        x_min = min(x_vals)
        x_max = max(x_vals)
        x_plot = np.linspace(x_min, x_max, 100)

        # Manejar diferentes tipos de polinomios
        if isinstance(poly, np.ndarray):  # Polinomio simple o mínimos cuadrados
            y_plot = np.polyval(poly, x_plot)
            label = f'Polinomio: {poly}'
        else:  # Polinomio de Lagrange (SymPy)
            y_plot = [poly.evalf(subs={'x': xi}) for xi in x_plot]
            label = str(poly)

        ax.plot(x_plot, y_plot, label=label)

        # Marcar el dato a aproximar
        x_approx = float(x_approx)
        if isinstance(poly, np.ndarray):
            y_approx = np.polyval(poly, x_approx)
        else:
            y_approx = poly.evalf(subs={'x': x_approx})
        ax.plot(x_approx, y_approx, 'ro', label=f'Dato aproximado: ({x_approx}, {y_approx:.2f})')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Interpolación y Ajuste de Curvas')
        ax.legend()

        # Crear un lienzo de tkinter para la figura de matplotlib
        self.canvas = FigureCanvasTkAgg(fig, master=self.graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def run_pol_simple(self):
        # Obtener los valores de los campos de entrada
        x_vals = list(map(float, self.x_entry.get().split(',')))
        y_vals = list(map(float, self.y_entry.get().split(',')))

        # Calcular el polinomio simple
        poly = Pol_simple(x_vals, y_vals)

        # Mostrar el resultado en el widget de texto
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"Polinomio Simple: {poly}\n")

        # Plotear la solución
        self.plot_solution(x_vals, y_vals, poly, self.x_approx_entry.get())

    def run_lagrange(self):
        # Obtener los valores de los campos de entrada
        x_vals = list(map(float, self.x_entry.get().split(',')))
        y_vals = list(map(float, self.y_entry.get().split(',')))

        # Calcular el polinomio de Lagrange
        poly = Lagrange(x_vals, y_vals)

        # Mostrar el resultado en el widget de texto
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"Polinomio de Lagrange: {poly}\n")

        # Plotear la solución
        self.plot_solution(x_vals, y_vals, poly, self.x_approx_entry.get())

    def run_min_cuad(self):
        # Obtener los valores de los campos de entrada
        x_vals = list(map(float, self.x_entry.get().split(',')))
        y_vals = list(map(float, self.y_entry.get().split(',')))

        # Calcular los coeficientes de mínimos cuadrados
        a0, a1 = min_c(x_vals, y_vals)

        # Crear el polinomio de mínimos cuadrados
        poly = np.array([a0, a1])

        # Mostrar el resultado en el widget de texto
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"Polinomio de Mínimos Cuadrados: y = {a0} + {a1}*x\n")

        # Plotear la solución
        self.plot_solution(x_vals, y_vals, poly, self.x_approx_entry.get())

