import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from func_bounds import *


def test_A(A, x, n):
    u = np.zeros(n + 1)
    for i in range(n + 1):
        u[i] = exact_soln(x[i])

    v = np.matmul(A, u)
    s = np.sum(v)

    if s == 0:
        print("\n Perfect. SUM A*u: \n", s)
    else:
        print("\n SUM A*u: \n", s)


def test_R(R, x, n):
    u = np.zeros(n + 1)
    for i in range(n + 1):
        u[i] = exact_soln(x[i])

    v = np.matmul(R, u)
    s = np.sum(v)

    f = lambda x : (3 - 5 * np.pi + np.pi ** 2) * x + (x ** 2 - 4 * x) * np.sin(x) - 1
    I = integrate.quad(f, 0, np.pi)

    # Check if s == the integral of u from a to b
    if s == I[0]:
        print("\n Perfect. SUM R*u: \n", s)
    else:
        print("\n SUM R*u: \n", s, "\n Integral of u: \n", I[0])

    rel_error = abs(I[0] - s) / abs(I[0])
    print("\n Relative Error in R: \n", rel_error)


def l2_norm(v, dx):
    val = 0
    for i in range(len(v)):
        val += (v[i]) ** 2 * dx
    return np.sqrt(val)


def L2_norm(f, R, x):
    val = 0
    for i in range(len(R)):
        for j in range(len(R[0])):
            val += ((f[i]) * (f[j]) * R[i, j])
    return np.sqrt(val)


def relative_error(u_comp, u_exact, R, x):
    return L2_norm(u_comp - u_exact, R, x) / L2_norm(u_exact, R, x)

    # RECALL:
    # L2_norm = SQRT( SUM (R * (u_comp - u_exact) ^ 2 ) ) / SQRT( SUM (R * u_exact ^ 2) )
    # Here, square refers to element wise squaring, not dot product squaring
    

def log_10_vals(arr):
    new_arr = np.zeros(len(arr))
    for i in range(len(arr)):
        new_arr[i] = np.log10(arr[i])
    return new_arr


def log_slope(error, n_val):
    dy = np.log10(error[-1]) - np.log10(error[-2])
    dx = np.log10(n_val[-1]) - np.log10(n_val[-2])
    return dy/dx


def regular_slope(f, n_val):
    dy = f[-1] - f[-2]
    dx = n_val[-1] - n_val[-2]
    return dy/dx