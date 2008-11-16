"""Method to create pre- and post-smoothers on the levels of a multilevel_solver
"""

__all__ = ['setup_smoothers']

def setup_smoothers(ml, presmoother, postsmoother):
    """Initialize pre- and post- smoothers throughout a multilevel_solver

    For each level of the multilevel_solver 'ml' (except the coarsest level),
    initialize the .presmoother() and .postsmoother() methods used in the 
    multigrid cycle.  Specifically, 'ml.

    Parameters
    ----------
    ml : multilevel_solver
        Data structure that stores the multigrid hierarchy.
    pre, post : smoother configuration
        See below of for available options
        
    Returns
    -------
    Nothing, ml will be changed in place.

    Smoother Configuration
    ----------------------
    Arguments 'pre' and 'post' can be the name of a supported smoother, 
    e.g. "gauss_seidel" or a tuple of the form 'method','opts') where 
    'method' is the name of a supported smoother and 'opts' a dict of
    keyword arguments to the smoother.  See the Examples section for
    illustrations of the format.

    Smoother Methods
    ----------------
    gauss_seidel
    jacobi
    richardson
    sor
    chebyshev
    kaczmarz_gauss_seidel
    kaczmarz_jacobi
    kaczmarz_richardson

    Notes
    -----
    Parameter 'omega' of the Jacobi and Richardson methods is scaled by the 
    spectral radius of the matrix on each level.  Therefore 'omega' should 
    be in the interval (0,2).

    Examples
    --------
    >>> from pyamg import poisson, smoothed_aggregation_solver
    >>> A = poisson((100,100), format='csr')
    >>> ml = smoothed_aggregation_solver(A)
    >>> #Gauss-Seidel with default options
    >>> setup_smoother(ml, pre='gauss_seidel', post='gauss_seidel')  
    >>> #Two iterations of symmetric Gauss-Seidel
    >>> pre  = ('gauss_seidel', {'iterations': 2, 'sweep':'symmetric'})
    >>> post = ('gauss_seidel', {'iterations': 2, 'sweep':'symmetric'})
    >>> setup_smoother(ml, pre=pre, post=post)

    """

    # interpret arguments
    if isinstance(presmoother, str):
        presmoother = (presmoother,{})
    if isinstance(postsmoother, str):
        postsmoother = (postsmoother,{})

    # get function handles
    try:
        setup_presmoother = eval('setup_' + presmoother[0])
    except NameError, ne:
        raise NameError("invalid presmoother method: ", presmoother[0])
    try:
        setup_postsmoother = eval('setup_' + postsmoother[0])
    except NameError, ne:
        raise NameError("invalid postsmoother method: ", postsmoother[0])

    for lvl in ml.levels[:-1]:
        lvl.presmoother  = setup_presmoother(lvl, **presmoother[1])
        lvl.postsmoother = setup_postsmoother(lvl, **postsmoother[1])


import relaxation
from chebyshev import chebyshev_polynomial_coefficients
from pyamg.util.utils import scale_rows
from pyamg.util.linalg import approximate_spectral_radius


def rho_D_inv_A(A):
    """Return the (approx.) spectral radius of D^-1 * A 
    """
    D = A.diagonal()
    D_inv = 1.0 / D
    D_inv[D == 0] = 0
    D_inv_A = scale_rows(A, D_inv, copy=True)
    return approximate_spectral_radius(D_inv_A)

def setup_gauss_seidel(lvl, iterations=1, sweep='forward'):
    def smoother(A,x,b):
        relaxation.gauss_seidel(A, x, b, iterations=iterations, sweep=sweep)
    return smoother

def setup_jacobi(lvl, iterations=1, omega=1.0):
    omega = omega/rho_D_inv_A(lvl.A)
    def smoother(A,x,b):
        relaxation.jacobi(A, x, b, iterations=iterations, omega=omega)
    return smoother

def setup_richardson(lvl, iterations=1, omega=1.0):
    omega = omega/approximate_spectral_radius(lvl.A)
    def smoother(A,x,b):
        relaxation.polynomial(A, x, b, coeffients=[omega], iterations=iterations)
    return smoother

def setup_sor(lvl, omega=0.5, iterations=1, sweep='forward'):
    def smoother(A,x,b):
        relaxation.sor(A, x, b, omega=omega, iterations=iterations, sweep=sweep)
    return smoother

def setup_chebyshev(lvl, lower_bound=1.0/30.0, upper_bound=1.1, degree=3, iterations=1):
    rho = approximate_spectral_radius(lvl.A)
    a = rho * lower_bound
    b = rho * upper_bound
    coeffients = -chebyshev_polynomial_coefficients(a, b, degree)[:-1] # drop the constant coefficient
    def smoother(A,x,b):
        relaxation.polynomial(A, x, b, coeffients=coeffients, iterations=iterations)
    return smoother

def setup_kaczmarz_jacobi(lvl, iterations=1, omega=1.0):
    omega = omega/rho_D_inv_A(lvl.A)**2
    def smoother(A,x,b):
        relaxation.kaczmarz_jacobi(A, x, b, iterations=iterations, omega=omega)
    return smoother

def setup_kaczmarz_gauss_seidel(lvl, iterations=1, sweep='forward'):
    def smoother(A,x,b):
        relaxation.kaczmarz_gauss_seidel(A, x, b, iterations=iterations, sweep=sweep)
    return smoother

def setup_kaczmarz_richardson(lvl, iterations=1, omega=1.0):
    omega = omega/approximate_spectral_radius(lvl.A)**2
    def smoother(A,x,b):
        relaxation.kaczmarz_richardson(A, x, b, iterations=iterations, omega=omega)
    return smoother


