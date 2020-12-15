import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as np_la
import scipy.sparse as sp
import scipy.sparse.linalg as sp_la
from matrices import *
from t_func_bounds import *
from t_test import *
from t_solver import *

def fem1d(a, b, n, bound, c1, c2,
          u_prev, 
          t_end, dt=1e-07,
          toll=1e-20, max_it=10000, 
          plot_me = True):
    
    #  Define the mesh with n+1 points between a and b.
    #  These will be x[0] through x[n].

    # This divides the interval [a,b] into n+1 equally spaced elements.
    x = np.linspace(a, b, n + 1)
    dx = x[1] - x[0]

    t = dt
    # Construct matrices A, R and the RHS vector.
    D, R = matrix(n, dx)
    
    u_all = []
    count=1
    while t <= t_end + dt/2:
        # Construct the vector representing f(x).
        f = np.zeros(n + 1)
        for i in range(n + 1):
            f[i] = func(x=x[i], t=t)

        # Test correctness of D and R.
        # test_D(D, x, n, t)
        # test_R(R, x, n, t)

        # Print D and R after the tests (but before boundary modifications).
        if False:
            print('\n\n', 'D := \n\n', D)
            print('\n\n', 'R := \n\n', R)


        A_initial, rhs_initial = modify_system(D, R, u_prev, f, dt)

        # Modify the linear system to enforce boundary conditions.
        if bound == 'dd':
            A, rhs = dual_dirichlet(A_initial, rhs_initial, n, c1(t), c2(t))
        elif bound == 'dn':
            A, rhs = dirichlet_neumann(A_initial, rhs_initial, n, c1(t), c2(t), dt)
        elif bound == 'nd':
            A, rhs = neumann_dirichlet(A_initial, rhs_initial, n, c1(t), c2(t), dt)
        elif bound == 'nn':
            A, rhs = neumann_neumann(A_initial, rhs_initial, n, c1(t), c2(t), dt)
        elif bound == 'rd':
            A, rhs = robin_dirichlet(A_initial, rhs_initial, n, c1(t), c2(t))

        # Print out the matrix A and the rhs vector b after boundary modifications.
        if False:
            print('\n\n', 'A := \n\n', A)
            print('\n\n', 'R := \n\n', R)
            print('\n\n', 'b := \n\n', rhs)

    
        # Solve the linear system
        u = modified_solver(A=A, b=rhs.transpose(), n=n, u_prev=u_prev, toll=toll, max_it=max_it)

        # Compute the exact solution u.
        u_ex = np.zeros(n + 1)
        for i in range(n + 1):
            u_ex[i] = exact_soln(x=x[i], t=t)

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
            u_ex_plot[i] = exact_soln(x=x_ex[i], t=t)
        
        # Plot the solution at every 1000th timestep
        if count%1000 == 0:
            u_all.append(u)
            print('Time Step Number ', count,' at t =', t, '\n')
            plt.plot(u)
        
        u_prev = u.copy()
        # Ensure that solution is computed at t_end regardless of dt chosen
        if t == t_end:
             break
        t += dt
        count+=1
        if t > t_end:
            t = t_end


    # At the end of the while loop t should be equal to t_end
    print("Relative Error: ", rel_error, "\nCurrent t: ", t, "     dt: ", dt, "   t_end: ", t_end)

    if plot_me == True:
        plt.show()
        plt.clf()
        
        plt.plot(x, u, 'b', x_ex, u_ex_plot, 'r')
        plt.savefig('u_final.png')
        plt.show()
    
    return rel_error
