from pyamg.testing import *

from scipy import matrix

from numpy import array, abs
from scipy import real, imag
from scipy.linalg import eigvals
from pyamg.gallery.laplacian import *

class TestPoisson(TestCase):
    def test_poisson(self):
        cases = []

        # 1D
        cases.append( ((1,),matrix([[2]])) )
        cases.append( ((2,),matrix([[ 2,-1],
                                    [-1, 2]])) )
        cases.append( ((4,),matrix([[ 2,-1, 0, 0],
                                    [-1, 2,-1, 0],
                                    [ 0,-1, 2,-1],
                                    [ 0, 0,-1, 2]])) )
        
        # 2D        
        cases.append( ((1,1), matrix([[4]])) )
        cases.append( ((2,1), matrix([[ 4,-1],
                                      [-1, 4]])) )
        cases.append( ((1,2), matrix([[ 4,-1],
                                      [-1, 4]])) )
        cases.append( ((1,3), matrix([[ 4,-1, 0],
                                      [-1, 4,-1],
                                      [ 0,-1, 4]])) )
        cases.append( ((2,2), matrix([[ 4,-1,-1, 0],
                                      [-1, 4, 0,-1],
                                      [-1, 0, 4,-1],
                                      [ 0,-1,-1, 4]])) )
        # 3D
        cases.append( ((2,2,1), matrix([[ 6,-1,-1, 0],
                                        [-1, 6, 0,-1],
                                        [-1, 0, 6,-1],
                                        [ 0,-1,-1, 6]])) )
        cases.append( ((2,2,2), matrix([[ 6,-1,-1, 0,-1, 0, 0, 0],
                                        [-1, 6, 0,-1, 0,-1, 0, 0],
                                        [-1, 0, 6,-1, 0, 0,-1, 0],
                                        [ 0,-1,-1, 6, 0, 0, 0,-1],
                                        [-1, 0, 0, 0, 6,-1,-1, 0],
                                        [ 0,-1, 0, 0,-1, 6, 0,-1],
                                        [ 0, 0,-1, 0,-1, 0, 6,-1],
                                        [ 0, 0, 0,-1, 0,-1,-1, 6]])) )

        for grid,expected in cases:
            assert_equal( poisson(grid).todense(), expected )


class TestGaugeLaplacian(TestCase):
    def test_gaugelaplacian(self):
        cases = []
        
        beta = 0.0; npts=5
        A = gauge_laplacian( npts, spacing=1.0, beta=beta)
        cases.append( (A, beta) )
        
        beta = 0.11; npts=3;
        A = gauge_laplacian( npts, spacing=1.0, beta=beta)
        cases.append( (A, beta) )
        npts=8
        A = gauge_laplacian( npts, spacing=1.0, beta=beta)
        cases.append( (A, beta) )
        A = gauge_laplacian( npts, spacing=0.1, beta=beta)
        cases.append( (A, beta) )

        beta = 0.42
        A = gauge_laplacian( npts, spacing=1.0, beta=beta)
        cases.append( (A, beta) )
        A = gauge_laplacian( npts, spacing=11.1, beta=beta)
        cases.append( (A, beta) )

        for A,beta in cases:
            # Check Hermitian
            diff = A - A.H
            assert_equal(diff.data, array([]))

            # Check for Definiteness
            e = eigvals(A.todense())
            if beta == 0.0:
                # Here, semi-definiteness
                assert_almost_equal(min(abs(e)), 0.0)
            else:
                # zero imaginary part
                assert_almost_equal( min(abs(imag(e))), 0.0)
                # positive real part
                assert( min(real(e)) > 0.0 )

