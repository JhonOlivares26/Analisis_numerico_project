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
        self.root.geometry("1000x600")

        control_frame = ttk.LabelFrame(self.root, text="Datos de Entrada")
        control_frame.pack(padx=10, pady=10, fill="both", expand=True)

        self.graph_frame = ttk.LabelFrame(self.root, text="Gráfica")
        self.graph_frame.pack(padx=10, pady=10, fill="both", expand=True)

        ttk.Label(control_frame, text="Función f(t, y):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.f_entry1 = ttk.Entry(control_frame, width=50)
        self.f_entry1.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Función g(t, y):").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.f_entry2 = ttk.Entry(control_frame, width=50)
        self.f_entry2.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Función h(t, y):").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.f_entry3 = ttk.Entry(control_frame, width=50)
        self.f_entry3.grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Valor inicial t (a):").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.t0_entry = ttk.Entry(control_frame)
        self.t0_entry.grid(row=3, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Valor final t (b):").grid(row=4, column=0, padx=5, pady=5, sticky="e")
        self.tf_entry = ttk.Entry(control_frame)
        self.tf_entry.grid(row=4, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Tamaño de paso h:").grid(row=5, column=0, padx=5, pady=5, sticky="e")
        self.h_entry = ttk.Entry(control_frame)
        self.h_entry.grid(row=5, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="Valores iniciales y0 (separados por comas):").grid(row=6, column=0, padx=5, pady=5, sticky="e")
        self.y0_entry = ttk.Entry(control_frame)
        self.y0_entry.grid(row=6, column=1, padx=5, pady=5)

        self.runge_kutta_button = ttk.Button(control_frame, text="Runge Kutta", command=self.run_runge_kutta)
        self.runge_kutta_button.grid(row=7, column=0, padx=5, pady=10, sticky="ew")

        self.euler_button = ttk.Button(control_frame, text="Euler", command=self.run_euler)
        self.euler_button.grid(row=7, column=1, padx=5, pady=10, sticky="ew")

        self.instrucciones_button = ttk.Button(control_frame, text="Instrucciones", command=self.mostrar_instrucciones)
        self.instrucciones_button.grid(row=7, column=2, padx=5, pady=10, sticky="ew")

        self.result_text = tk.Text(self.graph_frame, height=10, width=60)
        self.result_text.pack(side=tk.LEFT, fill="both", expand=True)

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill="both", expand=True)

    def plot_solution(self, t_vals, y_vals, method_names):
        self.ax.clear()
        if t_vals and y_vals:
            for i, method in enumerate(method_names):
                self.ax.plot(t_vals[i], y_vals[i], label=f"Función {method}")
            self.ax.set_xlabel('t')
            self.ax.set_ylabel('y')
            self.ax.set_title('Solución de Ecuaciones Diferenciales')
            self.ax.legend()
            self.canvas.draw()
        else:
            messagebox.showwarning("Advertencia", "No hay datos para graficar.")

    def parse_function(self, f_str):
        try:
            t, y = sp.symbols('t y')
            locals_dict = {"sp": sp}
            f_expr = sp.sympify(f_str, locals=locals_dict)
            return sp.lambdify((t, y), f_expr, 'numpy')
        except Exception as e:
            messagebox.showerror("Error", f"Error en la función: {e}")
            return None

    def run_runge_kutta(self):
        f_str1 = self.f_entry1.get().strip()
        f_str2 = self.f_entry2.get().strip()
        f_str3 = self.f_entry3.get().strip()
        try:
            t0 = float(self.t0_entry.get())
            tf = float(self.tf_entry.get())
            h = float(self.h_entry.get())
            y0 = list(map(float, self.y0_entry.get().split(',')))
        except ValueError as e:
            messagebox.showerror("Error", f"Error en los valores numéricos: {e}")
            return

        functions = []
        method_names = []
        initial_values = []

        if f_str1:
            f_lambda1 = self.parse_function(f_str1)
            if f_lambda1:
                functions.append(f_lambda1)
                method_names.append("f1")
                initial_values.append(y0[0])
        if f_str2:
            f_lambda2 = self.parse_function(f_str2)
            if f_lambda2:
                functions.append(f_lambda2)
                method_names.append("f2")
                initial_values.append(y0[1])
        if f_str3:
            f_lambda3 = self.parse_function(f_str3)
            if f_lambda3:
                functions.append(f_lambda3)
                method_names.append("f3")
                initial_values.append(y0[2])

        if not functions:
            messagebox.showerror("Error", "No se proporcionaron funciones válidas.")
            return

        t_vals, y_vals = [], []
        for func, y0_val in zip(functions, initial_values):
            t_val, y_val = RungeKutta(func, t0, tf, h, [y0_val])
            t_vals.append(t_val)
            y_vals.append(y_val)

        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"t: {t_vals}\n\ny: {y_vals}\n")
        self.plot_solution(t_vals, y_vals, method_names)

    def run_euler(self):
        f_str1 = self.f_entry1.get().strip()
        f_str2 = self.f_entry2.get().strip()
        f_str3 = self.f_entry3.get().strip()
        try:
            t0 = float(self.t0_entry.get())
            tf = float(self.tf_entry.get())
            h = float(self.h_entry.get())
            y0 = list(map(float, self.y0_entry.get().split(',')))
        except ValueError as e:
            messagebox.showerror("Error", f"Error en los valores numéricos: {e}")
            return

        functions = []
        method_names = []
        initial_values = []

        if f_str1:
            f_lambda1 = self.parse_function(f_str1)
            if f_lambda1:
                functions.append(f_lambda1)
                method_names.append("f1")
                initial_values.append(y0[0])
        if f_str2:
            f_lambda2 = self.parse_function(f_str2)
            if f_lambda2:
                functions.append(f_lambda2)
                method_names.append("f2")
                initial_values.append(y0[1])
        if f_str3:
            f_lambda3 = self.parse_function(f_str3)
            if f_lambda3:
                functions.append(f_lambda3)
                method_names.append("f3")
                initial_values.append(y0[2])

        if not functions:
            messagebox.showerror("Error", "No se proporcionaron funciones válidas.")
            return

        t_vals, y_vals = [], []
        for func, y0_val in zip(functions, initial_values):
            t_val, y_val = Euler(func, t0, tf, h, [y0_val])
            t_vals.append(t_val)
            y_vals.append(y_val)

        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"t: {t_vals}\n\ny: {y_vals}\n")
        self.plot_solution(t_vals, y_vals, method_names)

    def mostrar_instrucciones(self):
        instrucciones = """
        Instrucciones de Uso:
    
        1. Función f(t, y) - obligatorio:
            - Ingrese la función f(t, y) en este campo. Por ejemplo:
              - t + y
              - t**2 + y
              - y * sp.exp(t)
            - Utilice 'sp.' antes de funciones especiales de sympy, por ejemplo:
              - sp.exp(y)
              - sp.sin(t)
    
        2. Función g(t, y) - opcional:
            - Ingrese la función g(t, y) en este campo siguiendo las mismas reglas que para f(t, y).
    
        3. Función h(t, y) - opcional:
            - Ingrese la función h(t, y) en este campo siguiendo las mismas reglas que para f(t, y).
    
        4. Valor inicial t (a) - obligatorio:
            - Ingrese el valor inicial t en este campo. Debe ser un número real.
    
        5. Valor final t (b) - obligatorio:
            - Ingrese el valor final t en este campo. Debe ser un número real mayor que el valor inicial t.
    
        6. Tamaño de paso h - obligatorio:
            - Ingrese el tamaño de paso h en este campo. Debe ser un número real positivo.
    
        7. Valores iniciales y0 (separados por comas) - obligatorio:
            - Ingrese los valores iniciales y0 separados por comas. Cada valor se asignará a la función correspondiente en el orden en que se ingresaron.
    
        8. Botones:
            - Runge Kutta: Ejecuta el método de Runge Kutta con las funciones y parámetros ingresados.
            - Euler: Ejecuta el método de Euler con las funciones y parámetros ingresados.
            - Instrucciones: Muestra esta ventana de instrucciones.
    
        Nota: Asegúrese de ingresar funciones válidas y valores numéricos correctos antes de ejecutar los métodos.
        """

        messagebox.showinfo("Instrucciones de Uso", instrucciones)




