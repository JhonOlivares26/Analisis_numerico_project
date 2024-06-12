import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np

class MetodosNumericosApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Métodos Numéricos")

        self.create_widgets()

    def create_widgets(self):
        info_text = [
            "Ingrese la matriz A con los valores separados por comas y las filas separadas por llaves (ejemplo para matriz 3x3: [1,2,3],[2,3,1],[3,2,1])",
            "Ingrese los valores del vector B separados por coma (ejemplo: 3,2,1)",
            "Ingrese el vector inicial x0 con valores de 0 por cada fila del vector a (ejemplo para vector 3x3: 0,0,0)"
        ]
        self.input_frame = ttk.LabelFrame(self.root, text="Informacion de uso")
        self.input_frame.pack(padx=10, pady=10, fill="x")

        for i, text in enumerate(info_text):
            ttk.Label(self.input_frame, text=text).grid(row=i, column=0, padx=5, pady=5, sticky="w")


        self.input_frame = ttk.LabelFrame(self.root, text="Datos de Entrada")
        self.input_frame.pack(padx=10, pady=10, fill="x")

        ttk.Label(self.input_frame, text="Matriz A:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.A_entry = ttk.Entry(self.input_frame, width=50)
        self.A_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Vector B:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.b_entry = ttk.Entry(self.input_frame)
        self.b_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Vector inicial x0:").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.x0_entry = ttk.Entry(self.input_frame)
        self.x0_entry.grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Tolerancia:").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.tol_entry = ttk.Entry(self.input_frame)
        self.tol_entry.grid(row=3, column=1, padx=5, pady=5)

        self.button_frame = ttk.Frame(self.root)
        self.button_frame.pack(padx=10, pady=10, fill="x")

        self.eliminacion_gaussiana_button = ttk.Button(self.button_frame, text="Eliminación Gaussiana", command=self.run_eliminacion_gaussiana)
        self.eliminacion_gaussiana_button.grid(row=0, column=0, padx=5, pady=5)

        self.gauss_seidel_button = ttk.Button(self.button_frame, text="Gauss-Seidel", command=self.run_gauss_seidel)
        self.gauss_seidel_button.grid(row=0, column=1, padx=5, pady=5)

        self.output_frame = ttk.LabelFrame(self.root, text="Resultados")
        self.output_frame.pack(padx=10, pady=10, fill="both", expand=True)

        self.result_text = tk.Text(self.output_frame, height=10, width=70)
        self.result_text.pack(padx=5, pady=5, fill="both", expand=True)

    def run_eliminacion_gaussiana(self):
        self.result_text.delete(1.0, tk.END)
        try:
            A = np.array(eval("[" + self.A_entry.get() + "]"))
            b = np.array(eval("[" + self.b_entry.get() + "]"))
            result = Eliminacion_Gaussiana(A, b)
        except Exception as e:
            result = f"Error: {e}"
        self.result_text.insert(tk.END, result)

    def run_gauss_seidel(self):
        self.result_text.delete(1.0, tk.END)
        try:
            A = np.array(eval("[" + self.A_entry.get() + "]"))
            b = np.array(eval("[" + self.b_entry.get() + "]"))
            x0 = np.array(eval("[" + self.x0_entry.get() + "]"))
            tol = float(self.tol_entry.get())
            result = Gauss_seidel(A, b, x0, tol)
        except Exception as e:
            result = f"Error: {e}"
        self.result_text.insert(tk.END, result)

def Eliminacion_Gaussiana(A, b):
    n = len(b)
    x = np.zeros(n)

    for k in range(n-1):
        max_index = np.argmax(abs(A[k:, k])) + k  
        if max_index != k:
            A[[k, max_index]] = A[[max_index, k]] 
            b[[k, max_index]] = b[[max_index, k]]  
        
        for i in range(k+1, n):
            lam = A[i, k] / A[k, k]
            A[i, k:n] = A[i, k:n] - lam * A[k, k:n]
            b[i] = b[i] - lam * b[k]
    
    for k in range(n-1, -1, -1):
        x[k] = (b[k] - np.dot(A[k, k+1:n], x[k+1:n])) / A[k, k]

    return x

def Gauss_seidel(A, b, x0, tol):
    n = len(b)
    x1 = np.zeros(n)
    norm = tol + 1  
    iteraciones = 0

    while norm > tol:
        x0 = x1.copy() 
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i, j] * x1[j]
            x1[i] = (b[i] - sigma) / A[i, i]
        norm = np.linalg.norm(x1 - x0, np.inf)
        iteraciones += 1
    return x1
