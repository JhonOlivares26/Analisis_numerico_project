import tkinter as tk
import numpy as np

class MetodosNumericosApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Métodos Numéricos")
        root.geometry("600x400")

        self.A_entry = tk.Entry(root, width=50)
        self.b_entry = tk.Entry(root)
        self.x0_entry = tk.Entry(root)
        self.tol_entry = tk.Entry(root)

        tk.Label(root, text="Matriz A (separada por comas):").pack()
        self.A_entry.pack()
        tk.Label(root, text="Vector b (separado por comas):").pack()
        self.b_entry.pack()
        tk.Label(root, text="Vector inicial x0 (separado por comas):").pack()
        self.x0_entry.pack()
        tk.Label(root, text="Tolerancia:").pack()
        self.tol_entry.pack()

        self.eliminacion_gaussiana_button = tk.Button(root, text="Eliminación Gaussiana", command=self.run_eliminacion_gaussiana)
        self.gauss_seidel_button = tk.Button(root, text="Gauss-Seidel", command=self.run_gauss_seidel)

        self.eliminacion_gaussiana_button.pack(pady=5)
        self.gauss_seidel_button.pack(pady=5)

        self.result_text = tk.Text(root, height=10, width=70)
        self.result_text.pack(pady=10)

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

def Gauss_seidel(A, b, xo, tol):
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

