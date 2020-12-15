if __name__ == '__main__':
    import numpy as np
    from fem import fem1d

    def INSTRUCTIONS():
        '''

        INSTRUCTIONS:

        call fem1d(a, b, n, bound, c1, c2, toll, max_it, plot_me)
        
        a       ==> LHS of interval
        b       ==> RHS of interval
        c1      ==> value of exact solution u(x) at x=a
        c2      ==> value of exact solution u(x) at x=b
        
        bound argument:
        'dd'    ==> Dirichlet Dirichlet
        'dn'    ==> Dirichlet Neumann
        'nd'    ==> Neumann Dirichlet
        'rd'    ==> Robin Dirichlet
        
        toll    ==> the tollerance of the solver
        max_it  ==> the maximum number of iterations of the solver
        plot_me ==> boolean value: whether or not to plot the solution each time it is computed
        
        to specify a function f(x):
        open func_bounds.py and type the function into the 'func' method as shown below.

            # THIS IS THE FUNCTION f(x)
            def func(x):
                val = FUNCTION
            return val

        to specify an exact solution function u(x,t) (for graphing or test purposes):
        open func_bounds.py and type the function into the 'exact_soln' method as shown below.

            # THIS IS THE EXACT SOLUTION u(x,t)
            def exact_soln(x):
                val = FUNCTION
            return val

        :return:
        '''
'''
        TEST CASE:
            u(x) = (3 - 5*np.pi + np.pi**2)*x + (x**2 - 4*x)*np.sin(x) - 1
            f(x) = -4*x*np.cos(x) + 8*np.cos(x) - 2*np.sin(x) + x**2*np.sin(x) - 4*x*np.sin(x)
            u'(x) = 2*(x - 2)*np.sin(x) + (x - 4)*x*np.cos(x) + np.pi**2 - 5*np.pi + 3
'''

fem1d(a=0, b=np.pi, n=300, bound='nd', 
        c1=3 - 5 * np.pi + np.pi ** 2,
        c2=-1 + 3 * np.pi - 5 * np.pi ** 2 + np.pi ** 3,
        toll=1e-16, max_it=10000, plot_me=True)
