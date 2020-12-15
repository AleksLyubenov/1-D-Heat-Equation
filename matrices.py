import numpy as np


def matrix(n, dx):
    D = np.zeros((n + 1, n + 1))
    R = np.zeros((n + 1, n + 1))

    for i in range(1, n):
            D[i, i] = 2 / dx
            R[i, i] = (2 / 3) * dx
            D[i, i + 1] = -1 / dx
            R[i, i + 1] = dx / 6
            D[i, i - 1] = -1 / dx
            R[i, i - 1] = dx / 6

    D[0, 1] = -1/dx
    D[n, n-1] = -1/dx
    R[0, 1] = dx / 6
    R[n, n - 1] = dx / 6

    D[0, 0] = 1 / dx
    D[n, n] = 1 / dx
    R[0, 0] = dx / 3
    R[n, n] = dx / 3

    return D, R


def modify_system(D, R, u_prev, f, dt):
    A = R + dt * D
    rhs = np.matmul(R, u_prev) + dt * np.matmul(R, f)

    return A, rhs