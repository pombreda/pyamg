The source code for these examples is available in the [Downloads](http://code.google.com/p/pyamg/downloads/list) section.

Collection of demos highlighting features of PyAMG and examples of use.



# Blackbox Solver #

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Blackbox/demo.py)

This demo highlights using PyAMG's <tt>blackbox</tt> module, which is designed to solve an arbitrary system Ax=b.  The matrix A can be non-Hermitian, indefinite, Hermitian positive-definite, etc...  Generic and robust settings for smoothed\_aggregation\_solver(..) are used to invert A.  If <tt>blackbox.solve</tt> fails to solve your system well, look at the [Detecting Best Solver](#Detecting_Best_Solver.md) example below for guidance on automatically finding better parameter settings.

To use <tt>blackbox.solve</tt>, all you need is a matrix and a right-hand-side.
```
from numpy import arange, array           
from pyamg import solve               
from pyamg.gallery import poisson 
n = 100
A = poisson((n,n),format='csr')         
b = array(arange(A.shape[0]), dtype=float)
x = solve(A,b,verb=True,tol=1e-8)
```
produces the output
```
Detected a hermitian matrix
  maxiter = 400
  iteration 1
  iteration 2
  ...
  iteration 7
  iteration 8
Residual reduction ||r_k||/||r_0|| = 1.12e-07
```


---


# Smoothed Aggregation AMG #

## Aggregation ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Aggregation/demo.py)

The first level aggregates in AMG based on smoothed aggregation are depicted in this example.  An example mesh and adjacency mesh is loaded from square.mat, the smoothed\_aggregation\_solver is called, and the first level aggregates are plotted with draw.py.  Notice from the result, that most groupings encompass entire groups of elements in the underlying mesh.  Still, there are a few aggregates that yields "strings" in the aggregation, often impacting performance.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Aggregation.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Aggregation.png)

## One Dimensional Problem ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/OneDimension/demo.py)

This example is to illustrate the effect, in 1D, of smoothed aggregation on tentative prolongation operators.  Each of the aggregates (groups of three in this case) are plotted with their corresponding (smoothed) basis functions.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/OneDimension.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/OneDimension.png)

## Visualizing Aggregation ##

[demo1.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/VisualizingAggregation/demo1.py)

[demo2.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/VisualizingAggregation/demo1.py)

In these two example we show how to use the pyamg.vis module to display aggregation in both two and three dimensions.  demo1.py considers the Poisson problem on an unstructured triangulation of the unit square (from the pyamg gallery).  In demo2.py, the same Poisson problem is considered on an unstructured tetrahedral mesh on the unit cube.  Two VTK compliant output files are generated in each case: output\_mesh.vtu and output\_aggs.vtu.  output\_mesh.vtu provides information on the underlying mesh (straight from the unit square mesh), while output\_aggs.vtu holds information on the aggregates generated from first level of Smoothed Aggregation.  The process of visulization in paraview is straightforward:
  * start paraview:
  * open file: output\_mesh.vtu
  * apply
  * under display in the object inspector: select wireframe representation
  * under display in the object inspector: select a better solid color
  * open file: output\_aggs.vtu
  * apply
  * under display in the object inspector: select surface with edges representation
  * under display in the object inspector: select a better solid color
  * under display in the object inspector: increase line width to see line aggregates (if present)
  * under display in the object inspector: increase point size to see point aggregates (if present)

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/VisualizingAggregation1.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/VisualizingAggregation1.png)
![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/VisualizingAggregation2.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/VisualizingAggregation2.png)

## Detecting Best Solver ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/SolverDiagnostics/demo.py)

The range of parameter choices can often be overwhelming, especially to those new to AMG.  The solver diagnostics function makes finding good parameter choices easier.  A brute force search is applied, and depending on the matrix characteristics (e.g., symmetry and definiteness), 60-120 different solvers are constructed and then tested.  So, this is a script to be run overnight, or to be run on a small matrix.

The function solver\_diagnostics (solver\_diagnostics.py) has many parameters, but setting those is mostly for experts.  As a typical user, all you need is a matrix, A, which can be nonsymmetric, indefinite, or symmetric-positive-definite.  The function detects symmetry and definiteness, but it is safest to specify those.

The function output is two separate files and we briefly examine this output for the second example of rotated anisotropic diffusion from [demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/SolverDiagnostics/demo.py).

  1. The first output file is <tt>rot_ani_diff_diagnostic.txt</tt>, which is a sorted table of solver statistics for all the solvers tried.  This file has detailed output for the performance of each solver, and what parameter choices were used for each solver.
```
****************************************************************
*                Begin Solver Diagnostic Results               *
*                                                              *
*        ''solver #'' refers to below solver descriptors       *
*                                                              *
*        ''iters'' refers to iterations taken                  *
*                                                              *
*        ''op complexity'' refers to operator complexity       *
*                                                              *
*        ''work per DOA'' refers to work per digit of          *
*          accuracy to solve the algebraic system, i.e. it     *
*          measures the overall efficiency of the solver       *
****************************************************************


 solver # | iters | op complexity | work per DOA 
-------------------------------------------------
    1     |   11  |      1.5      |     3.1      
    2     |   10  |      1.5      |     3.3      
    3     |   12  |      1.5      |     3.3      
    4     |   10  |      1.7      |     3.7      
    ...
    ...
    ...
****************************************************************
*                 Begin Solver Descriptors                     *
****************************************************************

Solver Descriptor 1
  Solve phase arguments:
    cycle = V
    krylov accel = cg
    tol = 1e-08
    maxiter = 300
  Setup phase arguments:
    max_levels = 25
    max_coarse = 300
    coarse_solver = pinv
    presmoother = ('block_gauss_seidel', {'sweep': 'symmetric', 'iterations': 1})
    postsmoother = ('block_gauss_seidel', {'sweep': 'symmetric', 'iterations': 1})
    B = ones((A.shape[0],1), dtype=A.dtype); BH = B.copy()
    strength = ('symmetric', {'theta': 0.0})
    aggregate = standard
    smooth = ('energy', {'weighting': 'local', 'krylov': 'cg', 'degree': 3, 'maxiter': 4})
    Bimprove = default 
 
Solver Descriptor 2
  Solve phase arguments:
    cycle = V
    krylov accel = cg
    ...
    ...
    ...
```
  1. The second file defines a function, that when given a matrix, automatically generates and uses the best solver found.  For example,
```
>>> from scipy import pi
>>> from pyamg import gallery 
>>> from pyamg.gallery import diffusion 
>>> from rot_ani_diff_diagnostic import rot_ani_diff_diagnostic
>>> stencil = diffusion.diffusion_stencil_2d(type='FE', epsilon=0.001, >>> theta=2*pi/16.0)
>>> A = gallery.stencil_grid(stencil, (50,50), format='csr')
>>> rot_ani_diff_diagnostic(A)
```
> yields the output
```
multilevel_solver
Number of Levels:     2
Operator Complexity:  1.523
Grid Complexity:      1.116
Coarse Solver:        'pinv'
  level   unknowns     nonzeros
    0         2500        21904 [65.67%]
    1          289        11449 [34.33%]

System size:                (2500, 2500)
Avg. Resid Reduction:       0.16
Iterations:                 11
Operator Complexity:        1.52
Work per DOA:               3.15
Relative residual norm:     6.12e-09
```
> plus the convergence plot

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/solver_diagnostics.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/solver_diagnostics.png)

## Complex Arithmetic ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Complex/demo.py)

The smoothed aggregation solver natively supports complex arithmetc, i.e., there is no conversion to an equivalent real system.  For example, the highlighted demo here generates a basic gauge Laplacian from quantum chromodynamics and solves the system for a random right-hand-side and random initial guess.  The resulting on-screen output is
```
---Convergence Summary---------------------------------------

             Levels: 3
   Cycle Complexity:  2.644
Operator Complexity:  1.341
    Grid Complexity:  1.183
avg geo conv factor:  0.316
               work:  5.290

level   unknowns     nnz
   0    10000        50000      [74.57%]
   1    1658         15132      [22.57%]
   2    170          1920       [ 2.86%]

---Convergence Summary (verbose)-----------------------------
            Factors:
iter       Factor     A-Mean     G-Mean     Work      
0         
1          0.600      0.600      0.775      23.829    
2          0.185      0.392      0.480      8.301     
3          0.202      0.329      0.387      6.411     
4          0.240      0.307      0.352      5.825     
5          0.267      0.299      0.336      5.579     
6          0.283      0.296      0.328      5.457 
  ...
  ...
  ...    
```

## Nonsymmetric Example ##

[demo\_recirc\_flow.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Nonsymmetric/demo_recirc_flow.py)

The smoothed aggregation solver supports nonsymmetric (i.e., non-Hermitian) and indefinite matrices, through some recent advances in multigrid research.  The demo highlighted here constructs a solver for a small nonsymmetric recirculating flow problem.  Below, is an excerpt from the on-screen output
```
  ...
  ...
  ...
Now, we use the nonsymmetric solver to accelerate GMRES. 

---Convergence Summary---------------------------------------

             Levels: 3
   Cycle Complexity:  2.820
Operator Complexity:  1.466
    Grid Complexity:  1.422
avg geo conv factor:  0.315
               work:  5.629

level   unknowns     nnz
   0    225          1849       [68.23%]
   1    74           656        [24.21%]
   2    21           205        [ 7.56%]

---Convergence Summary (verbose)-----------------------------
iter       Factor     A-Mean     G-Mean     Work      
0         
1          0.168      0.168      0.410      7.274     
2          0.367      0.267      0.395      6.987     
3          0.218      0.251      0.340      6.023     
4          0.425      0.294      0.356      6.281     
5          0.544      0.344      0.382      6.743     
6          0.310      0.339      0.371      6.542     
  ...
  ...
  ...
```


---


# Classical AMG #
## Coarse Fine Splitting ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/CoarseFineSplitting/demo.py)

The C/F splitting---i.e. the splitting of indices into strictly coarse nodes and strictly fine nodes---using Ruge-Stuben coarsening is illustrated in this example.  An example mesh and adjacency graph is loaded from square.mat, ruge\_stuben\_solver is initiated, and the first level of splittings is plotted.  Coarse nodes are highlighted red, while fine nodes are highlighted blue.  Printing the multilevel object in this case shows that the coarsening is typical: around 25% reduction in unknowns:

```
>>> print mls
multilevel_solver
Number of Levels:     2
Operator Complexity:  1.331
Grid Complexity:      1.267
Coarse Solver:        'pinv2'
  level   unknowns     nonzeros
    0          191         1243 [75.15%]
    1           51          411 [24.85%]
```

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/CoarseFineSplitting.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/CoarseFineSplitting.png)

## Compatible Relaxation ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/CompatibleRelaxation/demo.py)

The C/F splitting---i.e. the splitting of indices into strictly coarse nodes and strictly fine nodes---using Compatible Relaxation is illustrated in this example.  A 2d finite-difference matrix of the Poisson problem is used and the coarse and fine splitting is plotted using matplotlib and the vis module in pyamg.  Coarse nodes are highlighted red, while fine nodes are highlighted blue.  In this case, the coarsening is not aggressive.

```
>>> import numpy
>>> len(numpy.where(splitting==0)[0])   #coarse
1186
>>> len(numpy.where(splitting==1)[0])   #fine
1314
>>> len(splitting)
2500
```

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/CompatibleRelaxation.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/CompatibleRelaxation.png)


---


# Rootnode AMG #

The rootnode\_solver is a mixture of both classical and aggregation-based approaches to AMG, with the intent to combine their strengths, while minimizing their respective drawbacks.  As a result, this solver is more robust for some problem types, especially anisotropic diffusion.

In terms of use, the interface to pyamg.aggregation.rootnode\_solver(...) is identical to pyamg.aggregation.smoothed\_aggregation\_solver(...), meaning that the above aggregation examples can be easily changed by simply replacing calls to "smoothed\_aggregation\_solver" with "rootnode\_solver."  However, we do provide two basic rootnode examples here.

## Visualization ##

[demo.py](http://code.google.com/p/pyamg/source/browse/trunk/Examples/Rootnode/demo.py)

This example compares the rootnode coarsening to classical AMG's coarsening ([Coarse Fine Splitting Example](#Coarse_Fine_Splitting.md)) and to smoothed aggregation's coarsening ([Aggregation Example](#Aggregation.md)).  The rootnode approach mixes classical AMG and smoothed aggregation, and hence has an associated C/F splitting that splits the indices into strictly coarse (C) nodes and strictly fine (F) nodes, and also has an associated aggregation that disjointly splits the nodes into strongly connected neighborhoods.  Essentially, each aggregate has one "root" C-node associated with it, that is injected between the fine and coarse grids.

An example mesh and adjacency graph is loaded from square.mat, and the rootnode\_solver is initiated.  Then, the first-level C/F splitting and the first-level aggregation are plotted. Coarse nodes are highlighted red, while fine nodes are highlighted blue. We can also print the multilevel object,
```
>>> print mls
multilevel_solver
Number of Levels:     2
Operator Complexity:  1.187
Grid Complexity:      1.131
Coarse Solver:        'pinv2'
  level   unknowns     nonzeros
    0          191         1243 [84.21%]
    1           25          233 [15.79%]
```

When we examine the C/F splitting, we find it to contain far fewer coarse nodes than for classical AMG ([Coarse Fine Splitting Example](#Coarse_Fine_Splitting.md)).  In general, this fewer number of coarse nodes is compensated by having a somewhat denser interpolation operator than for classical AMG.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/rootnode_CF.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/rootnode_CF.png)

Last, we examine the associated aggregation, finding each aggregate has one associated "root" C-node in red.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/rootnode_aggregation.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/rootnode_aggregation.png)

## Simple Example ##

The interface to rootnode\_solver is identical to that of smoothed\_aggregation\_solver, so we provide one simple example, noting that the above aggregation examples can be easily extended to the rootnode case.
```
>>> from pyamg import rootnode_solver, util
>>> from pyamg.gallery import poisson
>>> from scipy.sparse.linalg import cg
>>> import numpy
>>> A = poisson((100,100), format='csr')
>>> b = numpy.ones((A.shape[0]))
>>> ml = rootnode_solver(A, smooth=('energy', {'degree':2}), strength='evolution' )
>>> M = ml.aspreconditioner(cycle='V')
>>> x,info = cg(A, b, tol=1e-8, maxiter=30, M=M)  
>>> print "Final relative residual is  %1.2e"%( util.linalg.norm(b - A*x) / util.linalg.norm(b) )
```
The output should be similar to
```
Final relative residual is  8.38e-10
```


---


# Finite Elements #

## Diffusion ##

[demo\_anisotropic\_convergence.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Diffusion/demo_anisotropic_convergence.py)

[demo\_anisotropic\_scalability.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Diffusion/demo_anisotropic_scalability.py)

[demo\_strength\_measures.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Diffusion/demo_strength_measures.py)

There are three examples in this directory that illustrate typical measures for multigrid convergence.  First, demo\_anisotropic\_convergence.py, considers a rotated, anisotropic diffusion problem from the pyamg gallery.  The multigrid hierarchy is printed

```
             Levels: 7
   Cycle Complexity:  4.447
Operator Complexity:  2.224
    Grid Complexity:  1.707
avg geo conv factor:  0.442
               work: 12.535

level   unknowns     nnz
   0    10000        88804      [44.97%]
   1    5000         73014      [36.97%]
   2    1249         17929      [ 9.08%]
   3    624          14692      [ 7.44%]
   4    148          2520       [ 1.28%]
   5    38           468        [ 0.24%]
   6    9            45         [ 0.02%]
```

and shows high cycle complexities and high work-per-digit-accuracy for this problem (due to the rotated anisotropy).  The residual history that is plotted, indicates that convergence is eventually stable, with an average geometric convergence factor around 0.45.  The instantaneous and running average convergence factors are also printed in order to show that not all measures reflect performance.

Second, demo\_anisotropic\_scalability.py, considers again a rotated anisotropic diffusion problem from the pyamg gallery, but focus of this test is on scalability as the problem size increases.  The problem increases from 100<sup>2 to 600</sup>2 on a regular grid.  The results display the computed convergence factor and computed total work-per-digit-accuracy (and several other metrics):

```
    n     |    nnz    |   rho   |   OpCx  |   Work  
----------------------------------------------------
  10000   |   88804   |  0.370  |  2.224  |  5.157  
  40000   |   357604  |  0.342  |  2.245  |  4.820  
  90000   |   806404  |  0.343  |  2.263  |  4.871  
  160000  |  1435204  |  0.331  |  2.265  |  4.716  
  250000  |  2244004  |  0.327  |  2.267  |  4.668  
  360000  |  3232804  |  0.332  |  2.270  |  4.744
```

We see that scalability is acheived in this situation.  Tuning the parameters of the rotation and anisotropy impact scalability.

Third, demo\_strength\_measures.py, considers different strength measures in the SA-AMG setup phase.  In particular, the Classic Strength Measure is compared to the Evolution Measure.  For this example we see that total work is reduced by using the later measure as the basis for forming aggregates:

```
       Classic Strength Measure DropTol = 0.00       
    n     |    nnz    |   rho   |   OpCx  |   Work   
-----------------------------------------------------
  10000   |   88804   |  0.752  |  1.127  |  9.125   
  40000   |   357604  |  0.833  |  1.125  |  14.188  
  90000   |   806404  |  0.836  |  1.125  |  14.462  
  160000  |  1435204  |  0.842  |  1.125  |  15.043  
```

```
         ODE Strength Measure DropTol = 4.00         
    n     |    nnz    |   rho   |   OpCx  |   Work   
-----------------------------------------------------
  10000   |   88804   |  0.599  |  1.300  |  5.847   
  40000   |   357604  |  0.691  |  1.302  |  8.098   
  90000   |   806404  |  0.724  |  1.304  |  9.281   
  160000  |  1435204  |  0.746  |  1.306  |  10.262  
```

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Diffusion1.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Diffusion1.png)
![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Diffusion2.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Diffusion2.png)


## Linear Elasticity ##

[driver.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/LinearElasticity/driver.py)

We consider the 2D linear elasticity problem from the pyamg gallery in this example.  Three near null space modes are fed to the smoothed\_aggregation\_solver (relating to rotation and two types of translation).  SA-AMG is ideal for this problem and the results are apparent.  Very low operator complexities and the convergence is quick:

```
Number of Levels:     4
Operator Complexity:  1.280
Grid Complexity:      1.191
Coarse Solver:        'pinv2'
  level   unknowns     nonzeros
    0        80000      1430416 [78.10%]
    1        13467       356409 [19.46%]
    2         1587        40401 [ 2.21%]
    3          192         4356 [ 0.24%]
```

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/LinearElasticity.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/LinearElasticity.png)


## Dolfin Formulation ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/DolfinFormulation/demo.py)

This example considers a problem formed by the Fenics Dolfin finite element package and highlights the ease of using pyamg as a solver.  No memory copies are required to fetch the data from the Dolfin sparse matrix backend and subsequently to construct the corresponding scipy.sparse csr matrix.

## FiPy Formulation ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/FiPyFormulation/demo.py)

[PyAMGSolver.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/FiPyFormulation/PyAMGSolver.py)

We use the package FiPy in this example to construct a finite volume discretization of the diffusion problem with Dirichlet and Neumann boundary conditions.  A FiPy (inherited) class is provided for the scipy.sparse matrix used in PyAMG.  This example highlights the seamless integration possible with this package.  For this problem we achieve very good results:

```
MG iterations: 10
MG convergence factor: 0.195535
```

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/FiPyFormulation.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/FiPyFormulation.png)


---


# Preconditioning #

## Krylov Methods ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Preconditioning/demo.py)

This example shows how to effectively use a constructed multilevel solvers to precondtion a Krylov method.  The first example considers the Poisson problem from the pyamg gallery and uses a constant near-nullspace vector for SA-AMG.  The second example is 2D linear elasticity also from the pyamg gallery and use the typical three ridgid body modes (rotation and translation in x and y) to coach SA-AMG.  Since both problems are symmetric and positive definite, CG acceleration is used.  The residual histories show a clear improvement in using the SA-AMG preconditioners in both cases.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Preconditioning1.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Preconditioning1.png)
![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Preconditioning2.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Preconditioning2.png)


## Eigenvalue Solvers ##

[driver.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Eigensolver/driver.py)

We use in this example a Smoothed Aggregation AMG method to precondition the LOBPCG eigensolver to findthe lowest nine eigenmodes of a Poisson problem.  With preconditioning (M=M in the loppcg call), the computation of the eigensubspace is extremely fast.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Eigensolver.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Eigensolver.png)


---


# Other Applications #

## Graph Partitioning ##

[demo.py](http://code.google.com/p/pyamg/source/browse/branches/2.0.x/Examples/Partition/demo.py)

In this example we computing a partition of a basic crack mesh (crack\_mesh.mat) using the Fiedler vector (the second lowest eigenmode of the graph laplacian).  We construct a SA-AMG preconditioner to assist LOBPCG in finding the Fiedler vector.  Positive values of the Fiedler vector are plotted in red and negative values are plotted in blue, thus illustrating the natural splitting this mesh.

![http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Partition.png](http://pyamg.googlecode.com/svn/wiki/ExamplesResults/Partition.png)