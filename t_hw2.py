if __name__ == '__main__':
    import numpy as np
    from t_fem import *
    from t_func_bounds import *
    from t_test import *
    import matplotlib.pyplot as plt

    def INSTRUCTIONS():
        '''

        INSTRUCTIONS:

        call fem1d(a, b, n, bound, c1, c2, u_prev, t_end, dt, toll, max_it, plot_me)
        
        a       ==> LHS of interval
        b       ==> RHS of interval
        c1      ==> value of exact solution u(x) at x=a
        c2      ==> value of exact solution u(x) at x=b
        
        bound argument:
        'dd'    ==> Dirichlet Dirichlet
        'dn'    ==> Dirichlet Neumann
        'nd'    ==> Neumann Dirichlet
        'rd'    ==> Robin Dirichlet
        
        u_prev  ==> the initial condition: u(x,0) 
        dt      ==> the size of the timestep 
        t_end   ==> the final time at which the solution is computed
        toll    ==> the tollerance of the solver
        max_it  ==> the maximum number of iterations of the solver
        plot_me ==> boolean value: whether or not to plot the solution each time it is computed
        
        to specify a function f(x):
        open t_func_bounds.py and type the function into the 'func' method as shown below.

            # THIS IS THE FUNCTION f(x,t)
            def func(x,t):
                val = FUNCTION
            return val

        to specify an exact solution function u(x,t) (for graphing or test purposes):
        open t_func_bounds.py and type the function into the 'exact_soln' method as shown below.

            # THIS IS THE EXACT SOLUTION u(x,t)
            def exact_soln(x,t):
                val = FUNCTION
            return val

        :return:
        ''' 

'''
        TEST CASE:
            u(x,t) = np.cos(x) * (t**2) + np.sin(x**2)
            f(x,t) = 2*t*np.cos(x)-((t**2)*(-np.cos(x))-4*(x**2)*np.sin(x**2)+2*np.cos(x**2))
            u'(x,t) = 2*x*np.cos(x**2) - (t**2)*(sin(x))
'''

fem1d(a=0, b=np.pi, n=100, bound='dd', c1=lambda t: t**2,
                  c2=lambda t: np.sin(np.pi**2) - (t**2),
                  u_prev=initial_condition(0, np.pi, 100),
                  t_end=1e-07*10000, toll=1e-16,
                  dt=1e-07, max_it=10000,
                  plot_me = True)        
