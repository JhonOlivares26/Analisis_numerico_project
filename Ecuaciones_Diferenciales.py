import numpy as np
import matplotlib.pyplot as plt


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
