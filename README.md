NOTE: In all cases, if an expression for the exact solution is not specified: THE RELATIVE ERROR RETURNED WILL BE MEANINGLESS.

# STEADY STATE: -u_xx = f
        
        To run a single test case (that is, one computation with a chosen boundary condition), run hw2.py
        Simply change the parameters in the "fem1d" function as shown in the COMPUTATION section below.

        To run a complete test case (that is, every possible boundary condition at every n in {100, 200, ..., 1500},
        run error_plot_test_case.py
        Simply change the paramenters at the top of the script as below: 
                
                #VARIABLES FOR TEST CASE EXECUTION:
                a_input = 0
                b_input = np.pi
                toll_input=1e-16
                max_it_input = 10000
                plot_me_input = False
        
## THE COMPUTATION:
        
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
        plot_me ==> boolean value: whether or not to plot the solution when it is computed
        
        to specify a function f(x):
        open func_bounds.py and type the function into the 'func' method as shown below.

            # THIS IS THE FUNCTION f(x)
            def func(x):
                val = FUNCTION
            return val

        to specify an exact solution function u(x,t) (for graphing or test purposes):
        open func_bounds.py and type the function into the 'exact_soln' method as shown below.

            # THIS IS THE EXACT SOLUTION u(x)
            def exact_soln(x):
                val = FUNCTION
            return val
       
       The DEFAULT function choices are: 
                u(x) = (3 - 5*np.pi + np.pi**2)*x + (x**2 - 4*x)*np.sin(x) - 1
                u'(x) = 2*(x - 2)*np.sin(x) + (x - 4)*x*np.cos(x) + np.pi**2 - 5*np.pi + 3
                f(x) = -4*x*np.cos(x) + 8*np.cos(x) - 2*np.sin(x) + x**2*np.sin(x) - 4*x*np.sin(x)
       
       To simulate a test on the matrices A and R for an arbitrarty choice of functions:
                - The matrix A test requires no modification
                - The matrix R test requires the user to input their exact solution u(x) as a lambda expression in the 
                  test_R function of the test.py script
       
# TIME DEPENDENT EQUATION: u_t - u_xx = f

        To run a single test case (that is, one computation with a chosen boundary condition), run t_hw2.py
        Simply change the parameters in the "fem1d" function as shown in the COMPUTATION section below.
        
        To run a complete test case (that is, every possible boundary condition at every n in {100, 200, ..., 1500}, 
        run t_error_plot_test_case.py
        Simply change the paramenters at the top of the script as below: 
        
                #VARIABLES FOR TEST CASE EXECUTION:
                dt_input = 1e-07
                a_input = 0
                b_input = np.pi
                t_end_input = dt_input * 100
                toll_input=1e-16
                max_it_input = 10000
                plot_me_input = False
        
## THE COMPUTATION:

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
        
        to specify a function f(x):
        open t_func_bounds.py and type the function into the 'initial_condition' method as as a constant,
        or enter the exact solution u(x,t) as above and allow the initial condition function to evaluate
        at t=0.
                
            # THIS IS THE INITIAL CONDITION u(x,0)
            def initial_condition(a, b, n):
                x = np.linspace(a, b, n + 1)

                vec = np.zeros(n + 1)
                for i in range(len(vec)):
                        vec[i] = exact_soln(x=x[i], t=0)
                        #vec[i] = INITIAL CONDITION

                return vec
        
        The DEFAULT function choices are: 
                u(x,t) = np.cos(x) * (t**2) + np.sin(x**2)
                u'(x,t) = 2*x*np.cos(x**2) - (t**2)*(sin(x))
                f(x,t) = 2*t*np.cos(x)-((t**2)*(-np.cos(x))-4*(x**2)*np.sin(x**2)+2*np.cos(x**2))
                
        To simulate a test on the matrices A and R for an arbitrarty choice of functions:
        - The matrix A test requires no modification
        - The matrix R test requires the user to input their exact solution u(x,t) as a lambda expression in the 
          test_R function of the t_test.py script
                
                
