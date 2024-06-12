import tkinter as tk
from Ecuaciones_Diferenciales import DifferentialEquationsApp
from Ceros_de_funciones import CerosDeFuncionesApp
from Sistemas_de_Ecuaciones_Lineales import MetodosNumericosApp
from Series_de_taylor import SeriesTaylorApp
from Interpolacion_y_ajuste_de_curvas import InterpolationApp
import Ecuaciones_Diferenciales
import Interpolacion_y_ajuste_de_curvas
import Sistemas_de_Ecuaciones_Lineales


class MainApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Menú principal")
        root.geometry("600x300")

        def open_Series_de_taylor_app():
            new_window = tk.Toplevel(self.root)
            SeriesTaylorApp(new_window)

        def open_Ceros_de_funcione_app():
            new_window = tk.Toplevel(self.root)
            CerosDeFuncionesApp(new_window)

        def open_sistemas_de_ecuaciones_lineales_app():
            new_window = tk.Toplevel(self.root)
            MetodosNumericosApp(new_window)

        def open_Interpolacion_app():
            new_window = tk.Toplevel(self.root)
            InterpolationApp(new_window)

        button1 = tk.Button(root, text="Series de taylor", command=open_Series_de_taylor_app)
        button2 = tk.Button(root, text="Ceros de funciones", command=open_Ceros_de_funcione_app)
        button3 = tk.Button(root, text="Sistemas de ecuaciones lineales", command=open_sistemas_de_ecuaciones_lineales_app)
        button4 = tk.Button(root, text="Interpolación y ajuste de curvas", command=open_Interpolacion_app)
        button5 = tk.Button(root, text="Ecuaciones Diferenciales", command=self.open_differential_equations_app)

        button1.pack(pady=10)
        button2.pack(pady=10)
        button3.pack(pady=10)
        button4.pack(pady=10)
        button5.pack(pady=10)

    def open_differential_equations_app(self):
        new_window = tk.Toplevel(self.root)
        DifferentialEquationsApp(new_window)


root = tk.Tk()
app = MainApp(root)
root.mainloop()
