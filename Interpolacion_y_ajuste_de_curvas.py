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
        while max(np.abs((x1 - xo))) > tol:
            xo = x1
            x1 = np.dot(Tg, xo) + Cg
            iteraciones += 1
        return x1
    else:
        print("El sistema iterativo no converge a la solucion unica del sistema")
        return None

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
    x = sp.Symbol('x')
    n = len(xd)
    P = 0

    for i in range(n):
        L = 1  # Polinomio base L_i(x)
        for j in range(n):
            if j != i:
                L *= (x - xd[j]) / (xd[i] - xd[j])  # Construcción de L_i(x)
        P += L * yd[i]  # Sumar el término L_i(x) * y_i al polinomio P
    Poly = sp.expand(P)
    return Poly

def min_c(x, y):
    x = np.array(x)
    y = np.array(y)
    Sx = sum(x)
    Sf = sum(y)
    Sx2 = sum((x ** 2))
    Sxy = sum((x * y))
    n = len(x)
    a0 = (Sf * Sx2 - Sx * Sxy) / (n * Sx2 - Sx ** 2)
    a1 = (n * Sxy - Sx * Sf) / (n * Sx2 - Sx ** 2)
    return a0, a1

class InterpolationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Interpolación y Ajuste de Curvas")
        root.geometry("1000x600")  # Ajustar el tamaño de la ventana

        # Crear un marco para los controles
        control_frame = tk.Frame(root, width=300)
        control_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False)

        # Crear un marco para la gráfica
        self.graph_frame = tk.Frame(root, width=700)
        self.graph_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Crear etiquetas y campos de entrada para los parámetros
        tk.Label(control_frame, text="Puntos x (separados por comas):").pack()
        self.x_entry = tk.Entry(control_frame)
        self.x_entry.pack()

        tk.Label(control_frame, text="Puntos y (separados por comas):").pack()
        self.y_entry = tk.Entry(control_frame)
        self.y_entry.pack()

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

    def plot_solution(self, x_vals, y_vals, x_plot, y_plot, method_name):
        # Limpiar el marco de la gráfica
        for widget in self.graph_frame.winfo_children():
            widget.destroy()

        # Crear una nueva figura de matplotlib
        fig, ax = plt.subplots()
        ax.scatter(x_vals, y_vals, color='red', label='Datos')
        ax.plot(x_plot, y_plot, label=method_name)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Solución usando {method_name}')
        ax.legend()

        # Crear un lienzo de tkinter para la figura de matplotlib
        self.canvas = FigureCanvasTkAgg(fig, master=self.graph_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def run_pol_simple(self):
        # Limpiar el contenido del widget de texto
        self.result_text.delete(1.0, tk.END)

        # Obtener los valores de los campos de entrada
        x_vals = list(map(float, self.x_entry.get().split(',')))
        y_vals = list(map(float, self.y_entry.get().split(',')))

        # Calcular el polinomio simple
        coefs = Pol_simple(x_vals, y_vals)
        if coefs is None:
            self.result_text.insert(tk.END, "El sistema iterativo no converge a la solución única del sistema.\n")
            return

        # Mostrar los coeficientes del polinomio
        self.result_text.insert(tk.END, f"Coeficientes del polinomio simple: {coefs}\n")

        # Evaluar el polinomio en un rango de x para la gráfica
        x_plot = np.linspace(min(x_vals), max(x_vals), 100)
        y_plot = np.polyval(coefs[::-1], x_plot)
        self.plot_solution(x_vals, y_vals, x_plot, y_plot, "Polinomio Simple")

    def run_lagrange(self):
        # Limpiar el contenido del widget de texto
        self.result_text.delete(1.0, tk.END)

        # Obtener los valores de los campos de entrada
        x_vals = list(map(float, self.x_entry.get().split(',')))
        y_vals = list(map(float, self.y_entry.get().split(',')))

        # Calcular el polinomio de Lagrange
        Poly = Lagrange(x_vals, y_vals)

        # Mostrar el polinomio de Lagrange
        self.result_text.insert(tk.END, f"Polinomio de Lagrange: {Poly}\n")

        # Evaluar el polinomio en un rango de x para la gráfica
        x_plot = np.linspace(min(x_vals), max(x_vals), 100)
        y_plot = [Poly.evalf(subs={sp.Symbol('x'): val}) for val in x_plot]
        self.plot_solution(x_vals, y_vals, x_plot, y_plot, "Lagrange")

    def run_min_cuad(self):
        # Limpiar el contenido del widget de texto
        self.result_text.delete(1.0, tk.END)

        # Obtener los valores de los campos de entrada
        x_vals = list(map(float, self.x_entry.get().split(',')))
        y_vals = list(map(float, self.y_entry.get().split(',')))

        # Calcular los coeficientes de mínimos cuadrados
        a0, a1 = min_c(x_vals, y_vals)

        # Mostrar los coeficientes
        self.result_text.insert(tk.END, f"Coeficientes de mínimos cuadrados: a0 = {a0}, a1 = {a1}\n")

        # Evaluar la recta en un rango de x para la gráfica
        x_plot = np.linspace(min(x_vals), max(x_vals), 100)
        y_plot = a0 + a1 * x_plot
        self.plot_solution(x_vals, y_vals, x_plot, y_plot, "Mínimos Cuadrados")
