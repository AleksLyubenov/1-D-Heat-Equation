import scipy.sparse.linalg as spla
import scipy.sparse as sp
import numpy as np
import numpy.linalg as npl

def modified_solver(A, b, n, u_prev, toll, max_it):

    if False:
        # checking A is invertible:
        print("det(A) = ", npl.det(A))

        # checking A is definite positive
        print("Eigenvalues of A = ", npl.eig(A)[0])


    A = sp.csc_matrix(A)

    # First we create the ILU preconditioner
    
    M2 = spla.spilu(A) 
    # This computes the Incomplete LU factorization 
    # Instead of inverting ILU to get our preconditioner, we create a linear solver that solves ILUx=b and returns x
    # this is equivalent to x=(ILU)^-1*b
    
    M_x = lambda x: M2.solve(x) 
    M = spla.LinearOperator((n+1,n+1), M_x) # This creates the Incomplete LU factorization 

    # careful in your time dependent equation, x0 should be your previous iteration
    x, info = spla.bicgstab(A,b,x0=u_prev,M=M,tol=toll,maxiter=max_it) 
    # bicg stab is the BIConjugate Gradient STABilized
    # you can also use gmres instead but if bicgstab works keep it.

    if False:
        print("Solution")
        print("x, info = ",x, info)
        print("||Ax-b|| = ",npl.norm(A.dot(x)-b.flatten()))
    
    return x