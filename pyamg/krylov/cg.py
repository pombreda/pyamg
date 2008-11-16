from numpy import asmatrix, asarray, zeros, dot
from numpy import asanyarray, asarray, asmatrix, array, matrix, zeros
from scipy.linalg import norm
from scipy.linalg import solve as direct_solve
from scipy.sparse.linalg.interface import aslinearoperator, LinearOperator
from numpy.random import random
from warnings import warn
import time

_coerce_rules = {('f','f'):'f', ('f','d'):'d', ('f','F'):'F',
                 ('f','D'):'D', ('d','f'):'d', ('d','d'):'d',
                 ('d','F'):'D', ('d','D'):'D', ('F','f'):'F',
                 ('F','d'):'D', ('F','F'):'F', ('F','D'):'D',
                 ('D','f'):'D', ('D','d'):'D', ('D','F'):'D',
                 ('D','D'):'D'}

def coerce(x,y):
    if x not in 'fdFD':
        x = 'd'
    if y not in 'fdFD':
        y = 'd'
    return _coerce_rules[x,y]

def id(x):
    return x
def make_system(A, M, x0, b, xtype=None):
    """Make a linear system Ax=b
 
    Parameters
    ----------
    A : LinearOperator
        sparse or dense matrix (or any valid input to aslinearoperator)
    M : {LinearOperator, Nones}
        preconditioner
        sparse or dense matrix (or any valid input to aslinearoperator)
    x0 : {array_like, None}
        initial guess to iterative method
    b : array_like
        right hand side
    xtype : {'f', 'd', 'F', 'D', None}
        dtype of the x vector
 
    Returns
    -------
    (A, M, x, b, postprocess)
        A : LinearOperator
            matrix of the linear system
        M : LinearOperator
            preconditioner
        x : rank 1 ndarray
            initial guess
        b : rank 1 ndarray
            right hand side
        postprocess : function
            converts the solution vector to the appropriate
            type and dimensions (e.g. (N,1) matrix)
    """
    A_ = A
    A = aslinearoperator(A)
 
    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix (shape=%s)' % shape)
 
    N = A.shape[0]
 
    b = asanyarray(b)
 
    if not (b.shape == (N,1) or b.shape == (N,)):
        raise ValueError('A and b have incompatible dimensions')
 
    if b.dtype.char not in 'fdFD':
        b = b.astype('d') # upcast non-FP types to double
 
    def postprocess(x):
        if isinstance(b,matrix):
            x = asmatrix(x)
        return x.reshape(b.shape)
 
    if xtype is None:
        if hasattr(A,'dtype'):
            xtype = A.dtype.char
        else:
            xtype = A.matvec(b).dtype.char
        xtype = coerce(xtype, b.dtype.char)
    else:
        warn('Use of xtype argument is deprecated. '\
                'Use LinearOperator( ... , dtype=xtype) instead.',\
                DeprecationWarning)
        if xtype == 0:
            xtype = b.dtype.char
        else:
            if xtype not in 'fdFD':
                raise ValueError, "xtype must be 'f', 'd', 'F', or 'D'"
 
    b = asarray(b,dtype=xtype) #make b the same type as x
 
    if x0 is None:
        x = zeros(N, dtype=xtype)
    else:
        x = array(x0, dtype=xtype)
        if not (x.shape == (N,1) or x.shape == (N,)):
            raise ValueError('A and x have incompatible dimensions')
        x = x.ravel()
 
    # process preconditioner
    if M is None:
        if hasattr(A_,'psolve'):
            psolve = A_.psolve
        else:
            psolve = id
        if hasattr(A_,'rpsolve'):
            rpsolve = A_.rpsolve
        else:
            rpsolve = id
        M = LinearOperator(A.shape, matvec=psolve, rmatvec=rpsolve, dtype=A.dtype)
    else:
        M = aslinearoperator(M)
        if A.shape != M.shape:
            raise ValueError('matrix and preconditioner have different shapes')
 
    return A, M, x, b, postprocess


def cg(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None):
    """
    PCG from Saad
    """

    t = time.time()
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype=None)
    x = x.ravel()
    b = b.ravel()

    matvec = A.matvec
    psolve = M.matvec
    T0a = time.time()-t
    
    t = time.time()
    n = len(b)
    if maxiter is None:
        maxiter = n/10
    
    # setup
    doneiterating = False
    iter = 0

    r = b - matvec(x0)

    normr0 = norm(r)

    if normr0 < tol:
        doneiterating = True
    T0b = time.time()-t

    t = time.time()
    z = psolve(r)

    p = z.copy()
    Ap = matvec(p)

    T0c = time.time()-t
    print('T0a : %g ms',T0a*1000)
    print('T0b : %g ms',T0b*1000)
    print('T0c : %g ms',T0c*1000)

    T1=0
    T2=0
    T3=0
    T4=0
    T5=0
    T6=0
    T7=0

    while not doneiterating:
        t=time.time()
        alpha = dot(r.ravel(),z.ravel())/dot(Ap.ravel(),p.ravel())
        T1 += time.time()-t

        t=time.time()
        x = x + alpha * p
        T2 += time.time()-t

        t=time.time()
        rnew = r - alpha * Ap
        znew = psolve(rnew)
        T3 += time.time()-t

        t=time.time()
        beta = dot(rnew.ravel(),znew.ravel())/dot(r.ravel(),z.ravel())
        T4 += time.time()-t

        t=time.time()
        p = znew + beta * p
        T5 += time.time()-t

        t=time.time()
        z = znew
        r = rnew
        Ap = matvec(p)
        iter += 1
        T6 += time.time()-t

        t=time.time()
        normr = norm(r)
        T7 += time.time()-t

        if normr/normr0 < tol:
            doneiterating = True

        if iter>(maxiter-1):
            doneiterating = True

    print('T1 : %g ms',T1*1000)
    print('T2 : %g ms',T2*1000)
    print('T3 : %g ms',T3*1000)
    print('T4 : %g ms',T4*1000)
    print('T5 : %g ms',T5*1000)
    print('T6 : %g ms',T6*1000)
    print('T7 : %g ms',T7*1000)
    print('total : %g ms',(T0a+T0b+T0c+T1+T2+T3+T4+T5+T6+T7)*1000)
    return x

if __name__ == '__main__':
    # from numpy import diag
    # A = random((4,4))
    # A = A*A.transpose() + diag([10,10,10,10])
    # b = random((4,1))
    # x0 = random((4,1))

    from pyamg.gallery import stencil_grid
    A = stencil_grid([[0,-1,0],[-1,4,-1],[0,-1,0]],(3,3),dtype=float,format='csr')
    b = random((A.shape[0],))
    x0 = random((A.shape[0],))

    import time
    from scipy.sparse.linalg.isolve import cg as icg

    t1=time.time()
    x = cg(A,b,x0,tol=1e-8,maxiter=100)
    t2=time.time()
    print '%s took %0.3f ms' % ('cg', (t2-t1)*1000.0)
    print norm(b - A*x)
    print x.shape

    t1=time.time()
    y = icg(A,b,x0,tol=1e-8,maxiter=100)
    t2=time.time()
    print '%s took %0.3f ms' % ('icg', (t2-t1)*1000.0)
    print norm(b - A*y[0])
    print y[0]
