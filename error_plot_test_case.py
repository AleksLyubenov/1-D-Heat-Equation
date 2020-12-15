import numpy as np
import matplotlib.pyplot as plt
from test import *
from fem import fem1d


def DESCRIPTION():
    '''

        u(x) = (3 - 5*np.pi + np.pi**2)*x + (x**2 - 4*x)*np.sin(x) - 1
        f(x) = -4*x*np.cos(x) + 8*np.cos(x) - 2*np.sin(x) + x**2*np.sin(x) - 4*x*np.sin(x)
        u'(x) = 2*(x - 2)*np.sin(x) + (x - 4)*x*np.cos(x) + np.pi**2 - 5*np.pi + 3
        
        This test case will solve u_t - u_xx = f using the above functions
        at n = {100, 200, ..., 1500} with all possible boundary conditions. 
        
        The relative error will be plotted for each boundary condition at each value of n.
        
    :return:
    '''
    
x = np.linspace(0, np.pi, 100 + 1)

rel_1_err = []
rel_2_err = []
rel_3_err = []

#VARIABLES FOR TEST CASE EXECUTION:
a_input = 0
b_input = np.pi
toll_input=1e-16
max_it_input = 10000
plot_me_input = False

# Run relative error analysis for n = {100, 200, ..., 1500}
for i in range(100, 1600, 100):
    rel_1_err.append(fem1d(a=a_input, b=b_input, n=i, bound='dd', 
                            c1=-1, 
                            c2=-1 + 3 * np.pi - 5 * np.pi ** 2 + np.pi ** 3, 
                            toll=toll_input, max_it=max_it_input, 
                            plot_me=plot_me_input))
    rel_2_err.append(fem1d(a=a_input, b=b_input, n=i, bound='dn',
                            c1=-1, 
                            c2=3 - np.pi,
                            toll=toll_input, max_it=max_it_input, 
                            plot_me=plot_me_input))
    rel_3_err.append(fem1d(a=a_input, b=b_input, n=i, bound='nd', 
                            c1=3 - 5 * np.pi + np.pi ** 2, 
                            c2=-1 + 3 * np.pi - 5 * np.pi ** 2 + np.pi ** 3,
                            toll=toll_input, max_it=max_it_input, 
                            plot_me=plot_me_input))
log_rel_1_err = log_10_vals(rel_1_err)
log_rel_2_err = log_10_vals(rel_2_err)
log_rel_3_err = log_10_vals(rel_3_err)

n_val = []
for i in range(len(rel_1_err)):
    n_val.append(100*(i + 1))

line = n_val[::-1]
log_line = log_10_vals(line) + np.log10(rel_1_err[0])

print("\nLog Slope: Dirichlet Dirichlet = ", log_slope(rel_1_err, n_val))
print("Log Slope: Dirichlet Neumann = ", log_slope(rel_2_err, n_val))
print("Log Slope: Neumann Dirichlet = ", log_slope(rel_3_err, n_val))
print("\nLog Slope: Line = ", log_slope(line, n_val))
print("Reg Slope: Line = ", regular_slope(line, n_val))
print("\nReg Slope: Dirichlet Dirichlet = ", regular_slope(log_rel_1_err, n_val))
print("Reg Slope: Dirichlet Neumann = ", regular_slope(log_rel_2_err, n_val))
print("Reg Slope: Neumann Dirichlet = ", regular_slope(log_rel_3_err, n_val))


plt.plot(n_val, log_rel_1_err, 'b', label = "Dirichlet Dirichlet")
plt.plot(n_val, log_rel_2_err, 'g', label = "Dirichlet Neumann")
plt.plot(n_val, log_rel_3_err, 'y', label = "Neumann Dirichlet")
plt.plot(n_val, log_line, 'r', label = "Line")
plt.legend(loc="center right")

plt.show()
plt.savefig('steady_error_plot.png')
