import numpy as np


def dual_dirichlet(A, rhs, n, c1, c2):
    # Modify the linear system to enforce the left Dirichlet boundary condition: u(a) = c1
    A[0, 0] = 1.0
    A[0, 1:n] = 0.0
    rhs[0] = c1
    # Modify the linear system to enforce the right Dirichlet boundary condition: u(b) = c2
    A[n, n] = 1.0
    A[n, 0:n] = 0.0
    rhs[n] = c2
    return A, rhs


def dirichlet_neumann(A, rhs, n, c1, c2):
    # Modify the linear system to enforce the left Dirichlet boundary condition: u(a) = c1
    A[0, 0] = 1.0
    A[0, 1:n] = 0.0
    rhs[0] = c1
    # Modify the linear system to enforce the right Neumann condition: u'(b) = c2
    rhs[n] += c2
    return A, rhs


def neumann_dirichlet(A, rhs, n, c1, c2):
    # Modify the linear system to enforce the left Neumann condition: u'(a) = c1
    rhs[0] -= c1
    # Modify the linear system to enforce the right Dirichlet boundary condition: u(b) = c2
    A[n, n] = 1.0
    A[n, 0:n] = 0.0
    rhs[n] = c2
    return A, rhs


def robin_dirichlet(A, rhs, n, c1, c2):
    # Modify the linear system to enforce the left Robin condition: u(a) - u'(a) = c1
    A[0, 0] += 1.0
    rhs[0] += c1
    # Modify the linear system to enforce the right Dirichlet condition: u(b) = c2
    A[n, n] = 1.0
    A[n, 0:n] = 0.0
    rhs[n] = c2
    return A, rhs


# THIS IS THE EXACT SOLUTION u(x)
def exact_soln(x):
    val = (3 - 5 * np.pi + np.pi ** 2) * x + (x ** 2 - 4 * x) * np.sin(x) - 1
    return val


# THIS IS THE FUNCTION f(x)
def func(x):
    val = -4 * x * np.cos(x) + 8 * np.cos(x) - 2 * np.sin(x) + x ** 2 * np.sin(x) - 4 * x * np.sin(x)
    return val
