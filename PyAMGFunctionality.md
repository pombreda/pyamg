# Introduction #

Full documentation can be found [here](http://pyamg.googlecode.com/svn/branches/1.0.x/Docs/html/index.html).

Docstrings in the pyamg module offer a lot of direction.  For example TAB completing pyamg. gives a list of available modules.  The question mark, ?, gives help:
```
>>> pyamg.<TAB>
pyamg.Tester                       pyamg.aggregation
pyamg.__all__                      pyamg.amg_core
pyamg.__builtins__                 pyamg.bench
pyamg.__class__                    pyamg.classical
pyamg.__config__                   pyamg.coarse_grid_solver
pyamg.__delattr__                  pyamg.demo
pyamg.__dict__                     pyamg.gallery
pyamg.__doc__                      pyamg.graph
pyamg.__file__                     pyamg.info
pyamg.__getattribute__             pyamg.krylov
pyamg.__hash__                     pyamg.multilevel
pyamg.__init__                     pyamg.multilevel_solver
pyamg.__name__                     pyamg.relaxation
pyamg.__new__                      pyamg.ruge_stuben_solver
pyamg.__path__                     pyamg.show_config
pyamg.__reduce__                   pyamg.smoothed_aggregation_solver
pyamg.__reduce_ex__                pyamg.strength
pyamg.__repr__                     pyamg.test
pyamg.__setattr__                  pyamg.testing
pyamg.__str__                      pyamg.util
pyamg.__version__                  pyamg.version

>>> pyamg.smoothed_aggregation_solver?

Definition:       
pyamg.smoothed_aggregation_solver(A, B=None, BH=None, 
     symmetry='hermitian', strength='symmetric', aggregate='standard',
     smooth=('jacobi', {'omega': 1.3333333333333333}), 
     presmoother=('block_gauss_seidel', {'sweep': 'symmetric'}), 
     postsmoother=('block_gauss_seidel', {'sweep': 'symmetric'}),
     Bimprove='default', max_levels=10, max_coarse=500, keep=False,**kwargs)

Docstring:
    Create a multilevel solver using Smoothed Aggregation (SA)
.
.
.
```


# Details #

Summary of the most often used functions:
| **module**    | **submodule** | **funcion**                                 | **description** |
|:--------------|:--------------|:--------------------------------------------|:----------------|
|                              aggregation |                                 adaptive |                       adaptive\_sa\_solver | Create a multilevel solver using Adaptive Smoothed Aggregation (aSA) |
|                              aggregation |                                aggregate |                     standard\_aggregation | Compute the sparsity pattern of the tentative prolongator |
|                              aggregation |                                aggregate |                        lloyd\_aggregation | Aggregated nodes using Lloyd Clustering |
|                              aggregation |                                   smooth |             jacobi\_prolongation\_smoother | Jacobi prolongation smoother |
|                              aggregation |                                   smooth |         richardson\_prolongation\_smoother | Richardson prolongation smoother |
|                              aggregation |                                   smooth |             energy\_prolongation\_smoother | Minimize the energy of the coarse basis functions (columns of T) |
|                              aggregation |                                tentative |                           fit\_candidates | Fit near-nullspace candidates to form the tentative prolongator |
|                                classical |                                classical |                       ruge\_stuben\_solver | Create a multilevel solver using Classical AMG (Ruge-Stuben AMG) |
|                                classical |                                       cr |                                       CR | Use Compatible Relaxation to compute a C/F splitting  |
|                                classical |                                       cr |                              binormalize | Binormalize matrix A.  Attempt to create unit l\_1 norm rows. |
|                                classical |                              interpolate |                     direct\_interpolation | Create prolongator using direct interpolation |
|                                classical |                                    split |                                       RS | Compute a C/F splitting using Ruge-Stuben coarsening |
|                                classical |                                    split |                                     PMIS | C/F splitting using the Parallel Modified Independent Set method |
|                                classical |                                    split |                                    PMISc | C/F splitting using Parallel Modified Independent Set (in color) |
|                                classical |                                    split |                                      MIS | Compute a maximal independent set of a graph in parallel |
|                                  gallery |                               elasticity |                        linear\_elasticity | Linear elasticity problem discretizes with Q1 finite elements |
|                                  gallery |                               elasticity |                     linear\_elasticity\_p1 | P1 elements in 2 or 3 dimensions |
|                                  gallery |                                  example |                             load\_example | Load an example problem by name |
|                                  gallery |                                laplacian |                                  poisson | Returns a sparse matrix for the N-dimensional Poisson problem |
|                                  gallery |                                laplacian |                          gauge\_laplacian | Construct a Gauge Laplacian from Quantum Chromodynamics for regular 2D grids |
|                                  gallery |                            random\_sparse |                                   sprand | Returns a random sparse matrix. |
|                                  gallery |                                  stencil |                             stencil\_grid | Construct a sparse matrix form a local matrix stencil  |
|                                    graph |                                          |                  maximal\_independent\_set | Compute a maximal independent vertex set for a graph |
|                                    graph |                                          |                          vertex\_coloring | Compute a vertex coloring of a graph  |
|                                    graph |                                          |                             bellman\_ford | Bellman-Ford iteration |
|                                    graph |                                          |                            lloyd\_cluster | Perform Lloyd clustering on graph with weighted edges |
|                                    graph |                                          |                     connected\_components | Compute the connected components of a graph |
|                                   krylov |                                _bicgstab_|                                 bicgstab | Biconjugate Gradient Algorithm with Stabilization |
|                                   krylov |                                      _cg_|                                       cg | Conjugate Gradient algorithm |
|                                   krylov |                                    _cgne_|                                     cgne | Conjugate Gradient, Normal Error algorithm |
|                                   krylov |                                    _cgnr_|                                     cgnr | Conjugate Gradient, Normal Residual algorithm |
|                                   krylov |                                  _fgmres_|                                   fgmres | Flexible Generalized Minimum Residual Method (fGMRES) |
|                                   krylov |                                   _gmres_|                                    gmres | Generalized Minimum Residual Method (GMRES) |
|                               relaxation |                                chebyshev |        chebyshev\_polynomial\_coefficients | Chebyshev polynomial coefficients for the interval [a,b] |
|                               relaxation |                                chebyshev |              mls\_polynomial\_coefficients | Determine the coefficients for a MLS polynomial smoother |
|                               relaxation |                               relaxation |                                      sor | Perform SOR iteration on the linear system Ax=b |
|                               relaxation |                               relaxation |                             gauss\_seidel | Perform Gauss-Seidel iteration on the linear system Ax=b |
|                               relaxation |                               relaxation |                                   jacobi | Perform Jacobi iteration on the linear system Ax=b |
|                               relaxation |                               relaxation |                     gauss\_seidel\_indexed | Perform indexed Gauss-Seidel iteration on the linear system Ax=b |
|                               relaxation |                               relaxation |                               polynomial | Apply a polynomial smoother to the system Ax=b |
|                               relaxation |                                smoothing |                         change\_smoothers | Initialize pre- and post- smoothers throughout a multilevel\_solver, with |
|                                 strength |                                          |         classical\_strength\_of\_connection | Return a strength of connection matrix using the classical AMG measure |
|                                 strength |                                          |         symmetric\_strength\_of\_connection | Compute a strength of connection matrix using the standard symmetric measure |
|                                 strength |                                          |      energy\_based\_strength\_of\_connection | Compute a strength of connection matrix using an energy-based measure. |
|                                 strength |                                          |               ode\_strength\_of\_connection | `ode_strength_of_connection` is deprecated! |
|                                     util |                                   linalg |                                     norm | 2-norm of a vector |
|                                     util |                                   linalg |                            infinity\_norm | Infinity norm of a matrix (maximum absolute row sum).   |
|                                     util |                                   linalg |                            residual\_norm | Compute b - A\*x |
|                                     util |                                   linalg |              approximate\_spectral\_radius | Approximate the spectral radius of a matrix |
|                                     util |                                   linalg |                                  condest | Estimates the condition number of A |
|                                     util |                                   linalg |                                     cond | Returns condition number of A |
|                                     util |                                   linalg |                              ishermitian | Returns True if A is Hermitian to within tol |
|                                     util |                                    utils |                           profile\_solver | A quick solver to profile a particular multilevel object |
|                                     util |                                    utils |                              diag\_sparse | If A is a sparse matrix (e.g. csr\_matrix or csc\_matrix) |
|                                     util |                                    utils |                               scale\_rows | Scale the sparse rows of a matrix |
|                                     util |                                    utils |                            scale\_columns | Scale the sparse columns of a matrix |
|                                     util |                                    utils |                      symmetric\_rescaling | Scale the matrix symmetrically:: |
|                                     util |                                    utils |                             get\_diagonal | Return the diagonal or inverse of diagonal for  |
|                                     util |                                    utils |                              print\_table | Print a table from a list of lists representing the rows of a table |
|                                     util |                                    utils |                       hierarchy\_spectrum | Examine a multilevel hierarchy's spectrum |

| **class**     | **members**           |
|:--------------|:----------------------|
| multilevel  | levels              |
| multilevel  | coarse\_solver       |
| multilevel  | aspreconditioner    |
| multilevel  | cycle\_complexity    |
| multilevel  | grid\_complexity     |
| multilevel  | operator\_complexity |
| multilevel  | solve               |
| multilevel  | coarse\_grid\_solver  |
| level       | A                   |
| level       | R                   |
| level       | P                   |