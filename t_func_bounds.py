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


def dirichlet_neumann(A, rhs, n, c1, c2, dt):
    # Modify the linear system to enforce the left Dirichlet boundary condition: u(a) = c1
    A[0, 0] = 1.0
    A[0, 1:n] = 0.0
    rhs[0] = c1
    # Modify the linear system to enforce the right Neumann condition: u'(b) = c2
    rhs[n] += c2 * dt
    return A, rhs


def neumann_dirichlet(A, rhs, n, c1, c2, dt):
    # Modify the linear system to enforce the left Neumann condition: u'(a) = c1
    rhs[0] -= c1 * dt
    # Modify the linear system to enforce the right Dirichlet boundary condition: u(b) = c2
    A[n, n] = 1.0
    A[n, 0:n] = 0.0
    rhs[n] = c2
    return A, rhs


def neumann_neumann(A, rhs, n, c1, c2, dt):
    # Modify the linear system to enforce the left Neumann condition: u'(a) = c1
    rhs[0] -= c1 * dt
    # Modify the linear system to enforce the right Neumann condition: u'(b) = c2
    rhs[n] += c2 * dt
    return A, rhs


def robin_dirichlet(A, rhs, n, c1, c2, dt):
    # Modify the linear system to enforce the left Robin condition: u(a) - u'(a) = c1
    A[0, 0] += 1.0
    rhs[0] += c1 * dt
    # Modify the linear system to enforce the right Dirichlet condition: u(b) = c2
    A[n, n] = 1.0
    A[n, 0:n] = 0.0
    rhs[n] = c2
    return A, rhs


# THIS IS THE EXACT SOLUTION u(x)
def exact_soln(x, t):
    val = np.cos(x) * (t**2) + np.sin(x**2)
    return val


# THIS IS THE FUNCTION f(x)
def func(x, t):
    val = 2*t*np.cos(x)-((t**2)*(-np.cos(x))-4*(x**2)*np.sin(x**2)+2*np.cos(x**2))
    return val


# THIS IS THE INITIAL CONDITION u(x,0)
def initial_condition(a, b, n):
    x = np.linspace(a, b, n + 1)

    vec = np.zeros(n + 1)
    for i in range(len(vec)):
        vec[i] = exact_soln(x=x[i], t=0)
        #vec[i] = INITIAL CONDITION
        
    return vec
    
