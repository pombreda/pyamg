"""Method to create pre- and post-smoothers on the levels of a multilevel_solver
"""

__all__ = ['setup_smoothers', 'change_smoothers']#, 'smootherlist_tostring', 'smootherlists_tostring']


def setup_smoothers(ml, presmoother, postsmoother):
    """Initialize pre- and post- smoothers throughout a multilevel_solver

    For each level of the multilevel_solver 'ml' (except the coarsest level),
    initialize the .presmoother() and .postsmoother() methods used in the 
    multigrid cycle.

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
    e.g. "gauss_seidel" or a tuple of the form ('method','opts') where 
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
    cg
    gmres
    cgne
    cgnr
    None

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


# Helper fcn for unpacking smoother descriptions
def unpack_arg(v):
    if isinstance(v,tuple):
        return v[0],v[1]
    else:
        return v,{}

def change_smoothers(ml, presmoother, postsmoother):
    '''
    Initialize pre- and post- smoothers throughout a multilevel_solver, with
    the option of having different smoothers at different levels

    For each level of the multilevel_solver 'ml' (except the coarsest level),
    initialize the .presmoother() and .postsmoother() methods used in the 
    multigrid cycle.

    Parameters
    ----------
    ml : {pyamg multilevel hierarchy}
        Data structure that stores the multigrid hierarchy.
    presmoother : {None, string, tuple, list}
        presmoother can be (1) the name of a supported smoother, 
        e.g. "gauss_seidel", (2) a tuple of the form ('method','opts') where 
        'method' is the name of a supported smoother and 'opts' a dict of
        keyword arguments to the smoother, or (3) a list of instances of options 1 or 2.
        See the Examples section for illustrations of the format. 

        If presmoother is a list, presmoother[i] determines the smoothing
        strategy for level i.  Else, presmoother defines the same strategy 
        for all levels.
        
        If len(presmoother) < len(ml.levels), then 
        presmoother[-1] is used for all remaining levels
        
        If len(presmoother) > len(ml.levels), then
        the remaining smoothing strategies are ignored
    
    postsmoother : {string, tuple, list}
        Defines postsmoother in identical fashion to presmoother

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
    cg
    gmres
    cgne
    cgnr
    None

    Returns
    -------
    ml changed in place
    ml.level[i].presmoother  <===  presmoother[i]
    ml.level[i].postsmoother  <===  postsmoother[i]

    Notes
    -----
    Parameter 'omega' of the Jacobi and Richardson methods is scaled by the 
    spectral radius of the matrix on each level.  Therefore 'omega' should 
    be in the interval (0,2).

    Examples
    --------
    >>> from pyamg import poisson, smoothed_aggregation_solver
    >>> from pyamg.relaxation.smoothing import change_smoothers
    >>> from pyamg.util.linalg import norm
    >>> from scipy import rand
    >>> A = poisson((50,50), format='csr')
    >>> b = rand(A.shape[0],)
    >>> ml = smoothed_aggregation_solver(A, max_coarse=10)
    >>> smoothers = 'gauss_seidel'
    >>> change_smoothers(ml, presmoother=smoothers, postsmoother=smoothers)
    >>> x = ml.solve(b)
    >>> print "residaul 1 = %e" % norm(b - A*x) 
    >>>
    >>> smoothers = ('gauss_seidel', {'iterations' : 3})
    >>> change_smoothers(ml, presmoother=smoothers, postsmoother=None)
    >>> x = ml.solve(b)
    >>> print "residaul 2 = %e" % norm(b - A*x) 
    >>>
    >>> smoothers = ['gauss_seidel', ('cgnr', {'maxiter' : 5})]
    >>> change_smoothers(ml, presmoother=smoothers, postsmoother=smoothers)
    >>> x = ml.solve(b)
    >>> print "residaul 3 = %e" % norm(b - A*x) 
    '''

    # interpret arguments into list
    if isinstance(presmoother, str) or isinstance(presmoother, tuple) or (presmoother == None):
        presmoother = [presmoother]
    elif not isinstance(presmoother, list):
        raise ValueError,'Unrecognized presmoother'

    if isinstance(postsmoother, str) or isinstance(postsmoother, tuple) or (postsmoother == None):
        postsmoother = [postsmoother]
    elif not isinstance(postsmoother, list):
        raise ValueError,'Unrecognized postsmoother'

    # set ml.levels[i].presmoother = presmoother[i]
    for i in range( min(len(presmoother), len(ml.levels[:-1])) ):
        # unpack presmoother[i]
        fn,kwargs = unpack_arg(presmoother[i])
        # get function handle
        try:
            setup_presmoother = eval('setup_' + str(fn))
        except NameError, ne:
            raise NameError("invalid presmoother method: ", fn)
        ml.levels[i].presmoother = setup_presmoother(ml.levels[i], **kwargs)
    
    # Fill in remaining levels
    for j in range(i+1, len(ml.levels[:-1])):
        ml.levels[j].presmoother = setup_presmoother(ml.levels[j], **kwargs)


    # set ml.levels[i].postsmoother = postsmoother[i]
    for i in range( min(len(postsmoother), len(ml.levels[:-1])) ):
        # unpack postsmoother[i]
        fn,kwargs = unpack_arg(postsmoother[i])
        # get function handle
        try:
            setup_postsmoother = eval('setup_' + str(fn))
        except NameError, ne:
            raise NameError("invalid postsmoother method: ", fn)
        ml.levels[i].postsmoother = setup_postsmoother(ml.levels[i], **kwargs)
    
    # Fill in remaining levels
    for j in range(i+1, len(ml.levels[:-1])):
        ml.levels[j].postsmoother = setup_postsmoother(ml.levels[j], **kwargs)


from numpy import array, zeros, ones, ravel
import relaxation
from chebyshev import chebyshev_polynomial_coefficients
from pyamg.util.utils import scale_rows
from pyamg.util.linalg import approximate_spectral_radius
from pyamg.krylov import gmres, cgne, cgnr, cg

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

def setup_gmres(lvl, tol=1e-12, maxiter=1, restrt=None, M=None, callback=None, residuals=None):
    def smoother(A,x,b):
            (x[:],flag) = gmres(A, b, x0=x, tol=tol, maxiter=maxiter, restrt=restrt, M=M, callback=callback, residuals=residuals)
    return smoother
            
def setup_cg(lvl, tol=1e-12, maxiter=1, M=None, callback=None, residuals=None):
    def smoother(A,x,b):
            (x[:],flag) = cg(A, b, x0=x, tol=tol, maxiter=maxiter, M=M, callback=callback, residuals=residuals)
    return smoother

def setup_cgne(lvl, tol=1e-12, maxiter=1, M=None, callback=None, residuals=None):
    def smoother(A,x,b):
            (x[:],flag) = cgne(A, b, x0=x, tol=tol, maxiter=maxiter, M=M, callback=callback, residuals=residuals)
    return smoother

def setup_cgnr(lvl, tol=1e-12, maxiter=1, M=None, callback=None, residuals=None):
    def smoother(A,x,b):
            (x[:],flag) = cgnr(A, b, x0=x, tol=tol, maxiter=maxiter, M=M, callback=callback, residuals=residuals)
    return smoother

def setup_None(lvl):
    def smoother(A,x,b):
        pass
    return smoother

def smootherlist_tostring(sm):
    '''
    Parameters
    ----------
    sm : {string, tuple, list}
        A valid smoothing strategy input for change_smoothers.
        See change_smoothers documentation and Examples section for format examples.

    Returns
    -------
    string that when printed describes the smoothing strategy, sm, in readable format

    Notes
    -----
    Designed for running numerical tests on various pre/post smoothing 
    strategies so that the test scripts can generate readable output on 
    the multilevel method.

    Examples
    --------
    >>> from pyamg.relaxation.smoothing import smootherlist_tostring
    >>> sm = 'gauss_seidel'
    >>> print smootherlist_tostring(sm)
    >>> sm = ('gauss_seidel', {'iterations' : 3})
    >>> print smootherlist_tostring(sm)
    >>> sm = [('gauss_seidel', {'iterations' : 3}), ('gmres', {'maxiter' : 5})]
    >>> print smootherlist_tostring(sm)
    '''

    # interpret argument into list
    if isinstance(sm, str) or isinstance(sm, tuple) or (sm == None):
        sm = [sm]
    elif not isinstance(sm, list):
        return "Smoothing Strategy Unrecognized\n\n"
        
    out = ""    
    for i in range(len(sm)):
        fn, kwargs= unpack_arg(sm[i])
        
        if i == (len(sm)-1):
            out += "    Level %d to N Smoother = %s\n" % (i,fn)
        else:    
            out += "    Level %d Smoother = %s\n" % (i,fn)
        
        out += "      User Paramters:\n"
        for param in kwargs:
            out += '        ' + param + " = " + str(kwargs[param]) + '\n'
        
        if kwargs == {}:
            out += '        None\n' 

    return out


def smootherlists_tostring(smlist):
    '''
    Parameters
    ----------
    sm : {list}
        List of valid smoothing strategy inputs for change_smoothers
        See change_smoothers documentation and Examples section for format examples.
        
    Returns
    -------
    string that when printed describes the smoothing strategies 
    in smlist in readable format

    Notes
    -----
    Designed for running numerical tests on various pre/post smoothing 
    strategies so that the test scripts can generate readable output on 
    the multilevel method.

    Examples
    --------
    >>> from pyamg.relaxation.smoothing import smootherlists_tostring
    >>> sm = [('gauss_seidel', {'iterations' : 3}), ('gmres', {'maxiter' : 5})]
    >>> print smootherlists_tostring(sm)
    >>> sm = [ [('gauss_seidel', {'iterations' : 3}), 'gmres'], ('cgnr', {'maxiter' : 5})]
    >>> print smootherlists_tostring(sm)

    '''
    
    # interpret argument into list
    if isinstance(smlist, str) or isinstance(smlist, tuple) or (smlist == None):
        smlist = [smlist]
    elif not isinstance(smlist, list):
        return "Smoothing Strategy Unrecognized\n\n"

    outstr = ""
    for i in range(len(smlist)):
        outstr += "Smoothing Strategy %d\n" % (i+1)
        outstr += smootherlist_tostring(smlist[i])
        outstr += '\n'

    return outstr

