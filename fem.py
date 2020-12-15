import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_la
from matrices import matrix
from func_bounds import *
from test import *


def fem1d(a, b, n, bound, c1, c2, toll=1e-16, max_it=10000, plot_me = True):

    #  Define the mesh with n+1 points between a and b.
    #  These will be x[0] through x[n].

    # This divides the interval [a,b] into n+1 equally spaced elements.
    x = np.linspace(a, b, n + 1)
    dx = x[1] - x[0]

    # Construct matrices A, R and the RHS vector.
    A, R = matrix(n, dx)
    f = np.zeros(n + 1)
    for i in range(n + 1):
        f[i] = func(x[i])

    # Test correctness of A and R.
    if True:
        test_A(A, x, n)
        test_R(R, x, n)

    # Print A, R and b after the tests (but before boundary modifications).
    rhs = np.matmul(R, f)
    if False:
        print('\n\n', 'A := \n\n', A)
        print('\n\n', 'R := \n\n', R)
        print('\n\n', 'b := \n\n', rhs)


    # Modify the linear system to enforce boundary conditions.
    if bound == 'dd':
        A, rhs = dual_dirichlet(A, rhs, n, c1, c2)
    elif bound == 'dn':
        A, rhs = dirichlet_neumann(A, rhs, n, c1, c2)
    elif bound == 'nd':
        A, rhs = neumann_dirichlet(A, rhs, n, c1, c2)
    elif bound == 'rd':
        A, rhs = robin_dirichlet(A, rhs, n, c1, c2)


    # Print out the matrix A and the rhs vector b after boundary modifications.
    if False:
        print('\n\n', 'A := \n\n', A)
        print('\n\n', 'R := \n\n', R)
        print('\n\n', 'b := \n\n', rhs)


    # Solve the linear system
    A_sparse = sp.csr_matrix(A)
    
    #Conjugate gradient returns an answer and an integer to tell you if the algorithm has successfully executed.
    result = sp_la.bicg(A_sparse, rhs, tol=toll, maxiter=max_it)
    u = result[0]


    # Compute the exact solution u.
    u_ex = np.zeros(n + 1)
    for i in range(n + 1):
        u_ex[i] = exact_soln(x[i])

    # Compute the Relative Error of (u - u_ex).
    rel_error = relative_error(u, u_ex, R, x)

    if False:
        print('\n Node:     Computed u:        Exact u:          Error:      Rel_Error:\n')
        for i in range(0, n + 1):
            err = abs(u_ex[i] - u[i])
            print('  %4d  %14.6g  %14.6g  %14.6g %14.6g'
                  % (i, u[i], u_ex[i], err, rel_error))

    x_ex = np.linspace(a, b, 100)
    u_ex_plot = np.zeros(100)
    for i in range(0, 100):
        u_ex_plot[i] = exact_soln(x_ex[i])

    print('\nRelative Error: ', rel_error)

    if plot_me == True:
        plt.plot(x, u, 'b', x_ex, u_ex_plot, 'r.')
        plt.savefig('u.png')
        plt.show()

    return rel_error
