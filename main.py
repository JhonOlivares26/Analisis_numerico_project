import tkinter as tk


class MainApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Menú principal")
        root.geometry("600x300")

        # Definir funciones para los métodos numéricos
        def method1():
            print("Método numérico 1")

        def method2():
            print("Método numérico 2")

        # Crear botones para los métodos numéricos
        button1 = tk.Button(root, text="Método 1", command=method1)
        button2 = tk.Button(root, text="Método 2", command=method2)

        # Posicionar los botones
        button1.pack()
        button2.pack()


# Crear la ventana principal
root = tk.Tk()
app = MainApp(root)
root.mainloop()
