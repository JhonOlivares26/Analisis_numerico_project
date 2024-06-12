import tkinter as tk
from Ecuaciones_Diferenciales import DifferentialEquationsApp
from Ceros_de_funciones import CerosDeFuncionesApp
import Ecuaciones_Diferenciales
import Interpolacion_y_ajuste_de_curvas
import Sistemas_de_Ecuaciones_Lineales


class MainApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Menú principal")
        root.geometry("600x300")

        # Definir funciones para los métodos numéricos
        def method1():
            print("Series de taylor")

        def open_Ceros_de_funcione_app():
            new_window = tk.Toplevel(self.root)
            CerosDeFuncionesApp(new_window)

        def method3():
            print("Sistemas de ecuaciones lineales")

        def method4():
            print("Interpolación y ajuste de curvas")

        # Crear botones para los métodos numéricos
        button1 = tk.Button(root, text="Series de taylor", command=method1)
        button2 = tk.Button(root, text="Ceros de funciones", command=open_Ceros_de_funcione_app)
        button3 = tk.Button(root, text="Sistemas de ecuaciones lineales", command=method3)
        button4 = tk.Button(root, text="Interpolación y ajuste de curvas", command=method4)
        button5 = tk.Button(root, text="Ecuaciones Diferenciales", command=self.open_differential_equations_app)

        # Posicionar los botones
        button1.pack(pady=10)
        button2.pack(pady=10)
        button3.pack(pady=10)
        button4.pack(pady=10)
        button5.pack(pady=10)

    def open_differential_equations_app(self):
        # Crear una nueva ventana
        new_window = tk.Toplevel(self.root)
        # Crear una nueva instancia de DifferentialEquationsApp
        DifferentialEquationsApp(new_window)


# Crear la ventana principal
root = tk.Tk()
app = MainApp(root)
root.mainloop()
