import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk


def RungeKutta(f, a, b, h, y0):
    n = int((b - a) / h)
    t = np.linspace(a, b, n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        k1 = h * np.array(f(t[i], y[i]))
        k2 = h * np.array(f(t[i] + h / 2, y[i] + k1 / 2))
        k3 = h * np.array(f(t[i] + h / 2, y[i] + k2 / 2))
        k4 = h * np.array(f(t[i] + h, y[i] + k3))
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return t, y


def Euler(f, a, b, h, y0):
    n = int((b - a) / h)
    t = np.linspace(a, b, n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        y[i + 1] = y[i] + h * np.array(f(t[i], y[i]))
    return t, y


class DifferentialEquationsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Differential Equations")
        root.geometry("600x300")

        # Crear campos de entrada para los parámetros
        self.f_entry = tk.Entry(root)
        self.a_entry = tk.Entry(root)
        self.b_entry = tk.Entry(root)
        self.h_entry = tk.Entry(root)
        self.y0_entry = tk.Entry(root)

        # Crear botones para ejecutar las funciones
        self.runge_kutta_button = tk.Button(root, text="Runge Kutta", command=self.run_runge_kutta)
        self.euler_button = tk.Button(root, text="Euler", command=self.run_euler)

        # Posicionar los campos de entrada y los botones
        self.f_entry.pack()
        self.a_entry.pack()
        self.b_entry.pack()
        self.h_entry.pack()
        self.y0_entry.pack()
        self.runge_kutta_button.pack()
        self.euler_button.pack()

    def run_runge_kutta(self):
        # Obtener los valores de los campos de entrada
        f = self.f_entry.get()
        a = float(self.a_entry.get())
        b = float(self.b_entry.get())
        h = float(self.h_entry.get())
        y0 = list(map(float, self.y0_entry.get().split(',')))

        # Ejecutar la función RungeKutta y mostrar los resultados
        t, y = RungeKutta(f, a, b, h, y0)
        print(t, y)

    def run_euler(self):
        # Obtener los valores de los campos de entrada
        f = self.f_entry.get()
        a = float(self.a_entry.get())
        b = float(self.b_entry.get())
        h = float(self.h_entry.get())
        y0 = list(map(float, self.y0_entry.get().split(',')))

        # Ejecutar la función Euler y mostrar los resultados
        t, y = Euler(f, a, b, h, y0)
        print(t, y)


