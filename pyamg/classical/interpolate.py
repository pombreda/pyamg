"""Classical AMG Interpolation methods"""


__docformat__ = "restructuredtext en"

import numpy
from scipy.sparse import csr_matrix, isspmatrix_csr
from pyamg import amg_core

__all__ = ['direct_interpolation']


def direct_interpolation(A, C, splitting):
    """Create prolongator using direct interpolation

    Parameters
    ----------
    A : {csr_matrix}
        NxN matrix in CSR format
    C : {csr_matrix}
        Strength-of-Connection matrix
    splitting : array
        C/F splitting stored in an array of length N

    Returns
    -------
    P : {csr_matrix}
        Prolongator using direct interpolation

    Examples
    --------
    >>> from pyamg.gallery import poisson
    >>> from pyamg.classical import direct_interpolation
    >>> import numpy
    >>> A = poisson((5,),format='csr')
    >>> P = direct_interpolation(A,A,numpy.array([1,0,1,0,1]))
    >>> P.todense()
    matrix([[ 1. ,  0. ,  0. ],
            [ 0.5,  0.5,  0. ],
            [ 0. ,  1. ,  0. ],
            [ 0. ,  0.5,  0.5],
            [ 0. ,  0. ,  1. ]])


    """
    if not isspmatrix_csr(A): 
        raise TypeError('expected csr_matrix for A')

    if not isspmatrix_csr(C): 
        raise TypeError('expected csr_matrix for C')

    Pp = numpy.empty_like( A.indptr )

    amg_core.rs_direct_interpolation_pass1( A.shape[0],
            C.indptr, C.indices, splitting,  Pp)

    nnz = Pp[-1]
    Pj = numpy.empty( nnz, dtype=Pp.dtype )
    Px = numpy.empty( nnz, dtype=A.dtype )

    amg_core.rs_direct_interpolation_pass2( A.shape[0],
            A.indptr, A.indices, A.data,
            C.indptr, C.indices, C.data,
            splitting,
            Pp,       Pj,        Px)

    return csr_matrix( (Px,Pj,Pp) )
