First import PyAMG (the others are used later):
```
import scipy
import numpy
import pyamg
import pylab
```

Let's look at the 2D Q1 discretization of the Poisson problem.  Create a stencil and the matrix
```
stencil = [[-1,-1,-1],[-1,8,-1],[-1,-1,-1]]
A = pyamg.gallery.stencil_grid(stencil, (1000,1000), dtype=float, format='csr')
```

Now check out the matrix
```
A.nnz       # Number of nonzeros
A.shape     # Matrix dimensions
A.tocsr()   # Return A in Compressed Sparse Row (CSR)
A.tocsc()   # Return A in Compressed Sparse Column (CSC)
A.tocoo()   # Return A in Coordinate (COO) format
A.T         # Return the transpose of A
```

Construct standard aggregation for A.  Aggregation is one phase of AMG based on Smoothed Aggregation (SA).
```
Agg = pyamg.aggregation.standard_aggregation(A)
```

One can do a lot with this Agg (e.g. visualize, construct full solver).  Or we could build a multigrid hierarchy directly
```
B = numpy.ones((A.shape[0],1))
ml = pyamg.smoothed_aggregation_solver(A,B,max_coarse=10)
```
This tells PyAMG to create an SA hierarchy using A, a near-null space vector B of ones, and to coarsen until the system has 10 degrees-of-freedom.

Now print the details of the multilevel object:
```
print ml
```
which gives
```
multilevel_solver
Number of Levels:     7
Operator Complexity:  1.125
Grid Complexity:      1.126
Coarse Solver:        'pinv2'
  level   unknowns     nonzeros
    0      1000000      8988004 [88.87%]
    1       111556      1000000 [ 9.89%]
    2        12544       111556 [ 1.10%]
    3         1444        12544 [ 0.12%]
    4          169         1369 [ 0.01%]
    5           25          169 [ 0.00%]
    6            4           16 [ 0.00%]
```
Here we can see that 7 levels are constructed until a coarse level of 4 dof (<max\_coarse=10>) is achieved. The operator complexity (i.e., the sum of the nnz in the operator over all levels divided by the nnz in the fine level operator) is moderately low, and the grid complexity (i.e., the sum of the number of dof on all levels divided by the number of dof on the fine level) is low.  Also, the coarsest level solver is set to be application of the pseudo-inverse (pinv2).

Now, we can initialize a place-holder for the residual history, a right-hand-side for Ax=b, and an initial guess x0:
```
residuals=[]
b = scipy.rand(A.shape[0],1)
x0 = scipy.rand(A.shape[0],1)
```

Now, let's finally solve to a tolerance of 1e-10
```
x = ml.solve(b=b,x0=x0,tol=1e-10,residuals=residuals)
```

How well did we do?  Let's look at the geometric convergence factor:
```
(residuals[-1]/residuals[0])**(1.0/len(residuals))
```
which gives (approximately),
```
.16
```
in this case.  Not bad.  Now with CG acceleration we do
```
x = ml.solve(b=b,x0=x0,tol=1e-10,residuals=residuals,accel='cg')
```
which gives a geometric convergence factor of
```
(residuals[-1]/residuals[0])**(1.0/len(residuals))
```
which gives (approximately)
```
0.072
```

Plotting the residuals
```
pylab.semilogy(residuals/residuals[0],'o-')
pylab.xlabel('iterations')
pylab.ylabel('normalized residual')
pylab.show()
```

# Summary #

```
import scipy
import numpy
import pyamg
import pylab
stencil = [[-1,-1,-1],[-1,8,-1],[-1,-1,-1]]
A = pyamg.gallery.stencil_grid(stencil, (1000,1000), dtype=float, format='csr')
B = numpy.ones((A.shape[0],1))
ml = pyamg.smoothed_aggregation_solver(A,B,max_coarse=10)
print ml
residuals=[]
b = scipy.rand(A.shape[0],1)
x0 = scipy.rand(A.shape[0],1)
x = ml.solve(b=b,x0=x0,tol=1e-10,residuals=residuals)
(residuals[-1]/residuals[0])**(1.0/len(residuals))
x = ml.solve(b=b,x0=x0,tol=1e-10,residuals=residuals,accel='cg')
(residuals[-1]/residuals[0])**(1.0/len(residuals))
pylab.semilogy(residuals/residuals[0],'o-')
pylab.xlabel('iterations')
pylab.ylabel('normalized residual')
pylab.show()
```

## Result ##
![http://pyamg.googlecode.com/svn/wiki/Images/example_convergence.png](http://pyamg.googlecode.com/svn/wiki/Images/example_convergence.png)