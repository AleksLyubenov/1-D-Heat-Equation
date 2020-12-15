import numpy as np
import matplotlib.pyplot as plt
from t_fem import *
from t_func_bounds import *
from t_test import *


def DESCRIPTION():
    '''

        u(x,t) = np.cos(x) * (t**2) + np.sin(x**2)
        f(x,t) = 2*t*np.cos(x)-((t**2)*(-np.cos(x))-4*(x**2)*np.sin(x**2)+2*np.cos(x**2))
        u'(x,t) = 2*x*np.cos(x**2) - (t**2)*(sin(x))

        This test case will solve u_t - u_xx = f using the above functions
        at n = {100, 200, ..., 600} with all possible boundary conditions. 
        
        The relative error will be plotted for each boundary condition at each value of n.
        
    :return:
    '''


x = np.linspace(0, np.pi, 100 + 1)

rel_1_err = []
rel_2_err = []
rel_3_err = []
rel_4_err = []

#VARIABLES FOR TEST CASE EXECUTION:
dt_input = 1e-07
a_input = 0
b_input = np.pi
t_end_input = dt_input * 100
toll_input=1e-16
max_it_input = 10000
plot_me_input = False

# Run relative error analysis for n = {100, 200, ..., 600}
for i in range(100, 700, 100):
    rel_1_err.append(fem1d(a=a_input, b=b_input, n=i, bound='dd', 
                            c1=lambda t: t**2, 
                            c2=lambda t: np.sin(np.pi**2) - (t**2),
                            u_prev=initial_condition(0, np.pi, i),
                            t_end=t_end_input, toll=toll_input,
                            dt=dt_input, max_it=max_it_input,
                            plot_me=plot_me_input))

    rel_2_err.append(fem1d(a=a_input, b=b_input, n=i, bound='dn', 
                            c1=lambda t: t**2, 
                            c2=lambda t: 2*np.pi*np.cos(np.pi**2),
                            u_prev=initial_condition(0, np.pi, i),
                            t_end=t_end_input, toll=toll_input,
                            dt=dt_input, max_it=max_it_input,
                            plot_me=plot_me_input))
     
    rel_3_err.append(fem1d(a=a_input, b=b_input, n=i, bound='nd', 
                            c1=lambda t: 0, 
                            c2=lambda t: np.sin(np.pi**2) - (t**2),
                            u_prev=initial_condition(0, np.pi, i),
                            t_end=t_end_input, toll=toll_input,
                            dt=dt_input, max_it=max_it_input, 
                            plot_me=plot_me_input))

    rel_4_err.append(fem1d(a=a_input, b=b_input, n=i, bound='nn', 
                            c1=lambda t: 0,
                            c2=lambda t: 2*np.pi*np.cos(np.pi**2),
                            u_prev=initial_condition(0, np.pi, i),
                            t_end=t_end_input, toll=toll_input,
                            dt=dt_input, max_it=max_it_input, 
                            plot_me=plot_me_input))

log_rel_1_err = log_10_vals(rel_1_err)
log_rel_2_err = log_10_vals(rel_2_err)
log_rel_3_err = log_10_vals(rel_3_err)
log_rel_4_err = log_10_vals(rel_4_err)


n_val = []
for i in range(len(rel_1_err)):
    n_val.append(100 * (i + 1))

line = n_val[::-1]
log_line = log_10_vals(line) + np.log10(rel_1_err[0])

print("\nLog Slope: Dirichlet Dirichlet = ", log_slope(rel_1_err, n_val))
print("Log Slope: Dirichlet Neumann = ", log_slope(rel_2_err, n_val))
print("Log Slope: Neumann Dirichlet = ", log_slope(rel_3_err, n_val))
print("Log Slope: Neumann Neumann = ", log_slope(rel_4_err, n_val))
print("\nLog Slope: Line = ", log_slope(line, n_val))
print("Reg Slope: Line = ", regular_slope(line, n_val))
print("\nReg Slope: Dirichlet Dirichlet = ", regular_slope(log_rel_1_err, n_val))
print("Reg Slope: Dirichlet Neumann = ", regular_slope(log_rel_2_err, n_val))
print("Reg Slope: Neumann Dirichlet = ", regular_slope(log_rel_3_err, n_val))
print("Reg Slope: Neumann Neumann = ", regular_slope(log_rel_4_err, n_val))


plt.plot(n_val, log_rel_1_err, 'b', label="Dirichlet Dirichlet")
plt.plot(n_val, log_rel_2_err, 'g', label="Dirichlet Neumann")
plt.plot(n_val, log_rel_3_err, 'y', label="Neumann Dirichlet")
plt.plot(n_val, log_rel_4_err, 'k', label="Neumann Neumann")
plt.plot(n_val, log_line, 'r', label="Line")
plt.legend(loc="center right")

plt.show()
plt.savefig('timedep_error_plot.png')
