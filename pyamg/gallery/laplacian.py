"""Discretizations of the Poisson problem"""

__docformat__ = "restructuredtext en"

__all__ = ['poisson', "gauge_laplacian"]

import numpy
from numpy import array, abs
from scipy import arange, empty, intc, ravel, prod, pi, exp, hstack
from scipy.sparse import coo_matrix

def poisson( grid, spacing=None, dtype=float, format=None):
    """Finite Difference approximation to the Poisson problem on a 
    regular n-dimensional grid with Dirichlet boundary conditions.
   
    
    Parameters
    ----------
    grid : tuple of integers
        grid dimensions e.g. (100,100)


    Examples
    --------

    >>> # 4 nodes in one dimension
    >>> poisson( (4,) ).todense()
    matrix([[ 2., -1.,  0.,  0.],
            [-1.,  2., -1.,  0.],
            [ 0., -1.,  2., -1.],
            [ 0.,  0., -1.,  2.]])

    >>> # rectangular two dimensional grid 
    >>> poisson( (2,3) ).todense()
    matrix([[ 4., -1.,  0., -1.,  0.,  0.],
            [-1.,  4., -1.,  0., -1.,  0.],
            [ 0., -1.,  4.,  0.,  0., -1.],
            [-1.,  0.,  0.,  4., -1.,  0.],
            [ 0., -1.,  0., -1.,  4., -1.],
            [ 0.,  0., -1.,  0., -1.,  4.]])

    """
    grid = tuple(grid)

    D = len(grid) # grid dimension

    if D < 1 or min(grid) < 1:
        raise ValueError,'invalid grid shape: %s' % str(grid)

    nodes = arange(prod(grid)).reshape(*grid)

    nnz = nodes.size 
    for i in range(D):
        nnz += 2 * prod( grid[:i] + grid[i+1:] ) * (grid[i] - 1)
    
    row  = empty(nnz, dtype=intc)
    col  = empty(nnz, dtype=intc)
    data = empty(nnz, dtype=dtype)
    
    row[:nodes.size]  = ravel(nodes)
    col[:nodes.size]  = ravel(nodes)
    data[:nodes.size] = 2*D
    data[nodes.size:] = -1
    
    ptr = nodes.size
    
    for i in range(D):
        s0 = [slice(None)] * i + [slice(0,-1)  ] + [slice(None)] * (D - i - 1)
        s1 = [slice(None)] * i + [slice(1,None)] + [slice(None)] * (D - i - 1)
    
        n0 = nodes[s0]
        n1 = nodes[s1]
    
        row0 = row[ ptr:ptr + n0.size].reshape(n0.shape)
        col0 = col[ ptr:ptr + n0.size].reshape(n0.shape)
        ptr += n0.size
    
        row1 = row[ ptr:ptr + n0.size].reshape(n0.shape)
        col1 = col[ ptr:ptr + n0.size].reshape(n0.shape)
        ptr += n0.size
    
        row0[:] = n0
        col0[:] = n1
    
        row1[:] = n1
        col1[:] = n0
    
    return coo_matrix((data,(row,col)),shape=(nodes.size,nodes.size)).asformat(format)



def gauge_laplacian( npts, spacing=1.0, beta=0.1):
    ''' Construct a Gauge Laplacian from Quantum Chromodynamics for regualar 2D grids
        Note that this function is not written efficiently, but should be fine for N x N
        grids where N is in the low hundreds.

    Parameters
    ----------
    npts : {int}
        number of pts in x and y directions

    spacing : {float}
        grid spacing between points

    beta : {float}
        temperature
        Note that if beta=0, then we get the typical 5pt Laplacian stencil


    Examples
    --------
    $ A = gauge_laplacian(10)


    Output
    ------
    A : {csr matrix}
        A is Hermitian positive definite for beta > 0.0
        A is Symmetric semi-definite for beta = 0.0


    References
    ----------
    'Algebraic Multigrid Solvers for Complex-Valued Matrices", Maclachlan, Oosterlee, 
    Vol. 30, SIAM J. Sci. Comp, 2008
    '''

    # The gauge laplacian has the same sparsity structure as a normal
    # Laplacian, so we start out with a Poisson Operator
    N = npts
    A = poisson( (N,N),  format='coo', dtype=complex)

    # alpha is a random function of a point's integer position
    # on a 1-D grid along the x or y direction.  e.g. the first
    # point at (0,0) would be evaluate at alpha_*[0], while the
    # last point at (N*spacing, N*spacing) would evaluate at alpha_*[-1]
    alpha_x = 1.0j*2.0*pi*beta*numpy.random.randn(N*N)
    alpha_y = 1.0j*2.0*pi*beta*numpy.random.randn(N*N)
    
    # Replace off diagonals of A
    for i in range(A.nnz):
        r = A.row[i]; c = A.col[i]
        diff = abs(r - c)
        index = min(r, c)
        if r > c:
            s = -1.0
        else:
            s = 1.0
        if diff == 1:
            # differencing in the x-direction
            A.data[i] = -1.0*exp(s*alpha_x[index]) 
        if diff == N:
            # differencing in the y-direction
            A.data[i] = -1.0*exp(s*alpha_y[index])

    # Handle periodic BCs
    alpha_x = 1.0j*2.0*pi*beta*numpy.random.randn(N*N)
    alpha_y = 1.0j*2.0*pi*beta*numpy.random.randn(N*N)
    new_r = []; new_c = []; new_data = []; new_diff = []
    for i in range(0, N):
        new_r.append(i)
        new_c.append(i + N*N - N)
        new_diff.append(N)
    
    for i in range(N*N - N, N*N):
        new_r.append(i)
        new_c.append(i - N*N + N)
        new_diff.append(N)
    
    for i in range(0, N*N-1, N):
        new_r.append(i)
        new_c.append(i + N - 1)
        new_diff.append(1)
    
    for i in range(N-1, N*N, N):
        new_r.append(i)
        new_c.append(i - N + 1)
        new_diff.append(1)
    
    for i in range(len(new_r)):
        r = new_r[i]; c = new_c[i]; diff = new_diff[i]
        index = min(r, c)
        if r > c:
            s = -1.0
        else:
            s = 1.0
        if diff == 1:
            # differencing in the x-direction
            new_data.append(-1.0*exp(s*alpha_x[index])) 
        if diff == N:
            # differencing in the y-direction
            new_data.append(-1.0*exp(s*alpha_y[index]))
    
    # Construct Final Matrix
    A = coo_matrix( ( hstack((A.data, array(new_data))), (hstack((A.row, array(new_r))), hstack((A.col, array(new_c))) ) ), shape=(N*N,N*N) ).tocsr()

    return (1.0/spacing**2)*A

