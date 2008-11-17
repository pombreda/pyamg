from numpy import asmatrix, asarray, zeros, dot
from numpy import asanyarray, asarray, asmatrix, array, matrix, zeros
from numpy.random import random
from scipy.sparse.linalg.isolve.utils import make_system
from scipy.linalg import solve as direct_solve
from scipy.sparse.linalg.interface import aslinearoperator, LinearOperator
from pyamg.util.utils import norm
from warnings import warn
import time

def cg(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None):
    """
    Parameters
    ----------
    
    Returns
    -------    
    
    Notes
    -----

    References
    ----------
    PCG from Saad
    """

    t = time.time()
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype=None)
    x = x.ravel()
    b = b.ravel()

    T0a = time.time()-t
    
    t = time.time()
    n = len(b)
    if maxiter is None:
        maxiter = n/10
    
    # setup
    doneiterating = False
    iter = 0

    r = b - A*x0

    normr0 = norm(r)

    if normr0 < tol:
        doneiterating = True
    T0b = time.time()-t

    t = time.time()
    z = M*r

    p = z.copy()
    Ap = A*p

    T0c = time.time()-t
    print 'T0a : %g ms'%(T0a*1000)
    print 'T0b : %g ms'%(T0b*1000)
    print 'T0c : %g ms'%(T0c*1000)

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
        Ap = A*p
        iter += 1
        T6 += time.time()-t

        t=time.time()
        normr = norm(r)
        T7 += time.time()-t

        if normr/normr0 < tol:
            doneiterating = True

        if iter>(maxiter-1):
            doneiterating = True

    print 'T1 : %g ms'%(T1*1000)
    print 'T2 : %g ms'%(T2*1000)
    print 'T3 : %g ms'%(T3*1000)
    print 'T4 : %g ms'%(T4*1000)
    print 'T5 : %g ms'%(T5*1000)
    print 'T6 : %g ms'%(T6*1000)
    print 'T7 : %g ms'%(T7*1000)
    print 'total : %g ms'%((T0a+T0b+T0c+T1+T2+T3+T4+T5+T6+T7)*1000)
    print 'iterations: %d\n\n'%(iter)
    return postprocess(x)

if __name__ == '__main__':
    # from numpy import diag
    # A = random((4,4))
    # A = A*A.transpose() + diag([10,10,10,10])
    # b = random((4,1))
    # x0 = random((4,1))

    from pyamg.gallery import stencil_grid
    A = stencil_grid([[0,-1,0],[-1,4,-1],[0,-1,0]],(50,50),dtype=float,format='csr')
    b = random((A.shape[0],))
    x0 = random((A.shape[0],))

    import time
    from scipy.sparse.linalg.isolve import cg as icg

    print '\n\nTesting CG with %d x %d 2D Laplace Matrix'%(A.shape[0],A.shape[0])
    t1=time.time()
    x = cg(A,b,x0,tol=1e-8,maxiter=100)
    t2=time.time()
    print '%s took %0.3f ms' % ('cg', (t2-t1)*1000.0)
    print 'norm = %g'%(norm(b - A*x))

    t1=time.time()
    y = icg(A,b,x0,tol=1e-8,maxiter=100)
    t2=time.time()
    print '\n%s took %0.3f ms' % ('linalg cg', (t2-t1)*1000.0)
    print 'norm = %g'%(norm(b - A*y[0]))
    print 'info flag = %d'%(y[1])

    
