from pyamg.testing import *

from numpy import matrix, array, diag, zeros, sqrt, abs
from scipy import rand, real, imag, mat, diag
from scipy.sparse import csr_matrix
from scipy.linalg import eigvals, norm

from pyamg.utils import *
from pyamg.utils import symmetric_rescaling

class TestUtils(TestCase):
    def test_approximate_spectral_radius(self):
        cases = []

        cases.append( matrix([[-4]]) )
        cases.append( array([[-4]]) )
        
        cases.append( array([[2,0],[0,1]]) )
        cases.append( array([[-2,0],[0,1]]) )
      
        cases.append( array([[100,0,0],[0,101,0],[0,0,99]]) )
        
        for i in range(1,5):
            cases.append( rand(i,i) )
       
        # method should be almost exact for small matrices
        for A in cases:
            A = A.astype(float)
            Asp = csr_matrix(A)

            expected = max([norm(x) for x in eigvals(A)])
            assert_almost_equal( approximate_spectral_radius(A),   expected )
            assert_almost_equal( approximate_spectral_radius(Asp), expected )
        
        # try symmetric matrices
        for A in cases:
            A = A + A.transpose()
            A = A.astype(float) 
            Asp = csr_matrix(A)

            expected = max([norm(x) for x in eigvals(A)])
            assert_almost_equal( approximate_spectral_radius(A,symmetric=True),   expected )
            assert_almost_equal( approximate_spectral_radius(Asp,symmetric=True), expected )
      
        #TODO test larger matrices
    
    def test_infinity_norm(self):
        A = matrix([[-4]])
        assert_equal(infinity_norm(csr_matrix(A)),4)

        A = matrix([[1,0,-5],[-2,5,0]])
        assert_equal(infinity_norm(csr_matrix(A)),7)

        A = matrix([[0,1],[0,-5]])
        assert_equal(infinity_norm(csr_matrix(A)),5)

        A = matrix([[1.3,-4.7,0],[-2.23,5.5,0],[9,0,-2]])
        assert_equal(infinity_norm(csr_matrix(A)),11)

    def test_diag_sparse(self):
        #check sparse -> array
        A = matrix([[-4]])
        assert_equal(diag_sparse(csr_matrix(A)),[-4])

        A = matrix([[1,0,-5],[-2,5,0]])
        assert_equal(diag_sparse(csr_matrix(A)),[1,5])

        A = matrix([[0,1],[0,-5]])
        assert_equal(diag_sparse(csr_matrix(A)),[0,-5])

        A = matrix([[1.3,-4.7,0],[-2.23,5.5,0],[9,0,-2]])
        assert_equal(diag_sparse(csr_matrix(A)),[1.3,5.5,-2])

        #check array -> sparse
        A = matrix([[-4]])
        assert_equal(diag_sparse(array([-4])).todense(),csr_matrix(A).todense())

        A = matrix([[1,0],[0,5]])
        assert_equal(diag_sparse(array([1,5])).todense(),csr_matrix(A).todense())

        A = matrix([[0,0],[0,-5]])
        assert_equal(diag_sparse(array([0,-5])).todense(),csr_matrix(A).todense())

        A = matrix([[1.3,0,0],[0,5.5,0],[0,0,-2]])
        assert_equal(diag_sparse(array([1.3,5.5,-2])).todense(),csr_matrix(A).todense())


    def test_symmetric_rescaling(self):
        cases = []
        cases.append( diag_sparse(array([1,2,3,4])) )
        cases.append( diag_sparse(array([1,0,3,4])) )

        A = array([ [ 5.5,  3.5,  4.8],
                    [ 2. ,  9.9,  0.5],
                    [ 6.5,  2.6,  5.7]])
        A = csr_matrix( A )
        cases.append(A)
        P = diag_sparse([1,0,1])
        cases.append( P*A*P )
        P = diag_sparse([0,1,0])
        cases.append( P*A*P )
        P = diag_sparse([1,-1,1])
        cases.append( P*A*P )

        for A in cases:
            D_sqrt,D_sqrt_inv,DAD = symmetric_rescaling(A)

            assert_almost_equal( diag_sparse(A) > 0, diag_sparse(DAD) )
            assert_almost_equal( diag_sparse(DAD), D_sqrt*D_sqrt_inv )

            D_sqrt,D_sqrt_inv = diag_sparse(D_sqrt),diag_sparse(D_sqrt_inv)
            assert_almost_equal((D_sqrt_inv*A*D_sqrt_inv).todense(), DAD.todense())


    def test_profile_solver(self):
        from scipy.sparse.linalg import cg
        from pyamg.gallery import poisson
        from pyamg.aggregation import smoothed_aggregation_solver

        A = poisson((100,100), format='csr')
        ml = smoothed_aggregation_solver(A)

        opts = []
        opts.append( {} )
        opts.append( {'accel' : cg } )
        opts.append( {'accel' : cg, 'tol' : 1e-10 } )

        for kwargs in opts:
            residuals = profile_solver(ml, **kwargs)


class TestComplexUtils(TestCase):
    def test_approximate_spectral_radius(self):
        cases = []

        cases.append( matrix([[-4-4.0j]]) )
        cases.append( array([[-4+8.2j]]) )
        
        cases.append( array([[2.0-2.9j,0],[0,1.5]]) )
        cases.append( array([[-2.0-2.4j,0],[0,1.21]]) )
      
        cases.append( array([[100+1.0j,0,0],[0,101-1.0j,0],[0,0,99+9.9j]]) )
        
        for i in range(1,5):
            cases.append( rand(i,i)+1.0j*rand(i,i) )
       
        # method should be almost exact for small matrices
        for A in cases:
            Asp = csr_matrix(A)
            e = eigvals(A)
            expected = max(abs(e))
            assert_almost_equal( approximate_spectral_radius(A),   expected )
            assert_almost_equal( approximate_spectral_radius(Asp), expected )
            
            AA = mat(A).H*mat(A)
            AAsp = csr_matrix(AA)
            e = eigvals(AA)
            expected = max(abs(e))
            assert_almost_equal( approximate_spectral_radius(AA, symmetric=True),   expected )
            assert_almost_equal( approximate_spectral_radius(AAsp, symmetric=True), expected )
 
    def test_infinity_norm(self):
        A = matrix([[-4-3.0j]])
        assert_equal(infinity_norm(csr_matrix(A)),5.0)

        A = matrix([[1,0,4.0-3.0j],[-2,5,0]])
        assert_equal(infinity_norm(csr_matrix(A)),7)

        A = matrix([[0,1],[0,-4.0+3.0j]])
        assert_equal(infinity_norm(csr_matrix(A)),5.0)

    def test_diag_sparse(self):
        #check sparse -> array
        A = matrix([[-4-4.0j]])
        assert_equal(diag_sparse(csr_matrix(A)),[-4-4.0j])

        A = matrix([[1,0,-5],[-2,5-2.0j,0]])
        assert_equal(diag_sparse(csr_matrix(A)),[1,5-2.0j])

        #check array -> sparse
        A = matrix([[-4+1.0j]])
        assert_equal(diag_sparse(array([-4+1.0j])).todense(),csr_matrix(A).todense())

        A = matrix([[1,0],[0,5-2.0j]])
        assert_equal(diag_sparse(array([1,5-2.0j])).todense(),csr_matrix(A).todense())

    def test_symmetric_rescaling(self):
        cases = []
        A = array([ [ 5.5+1.0j,  3.5,    4.8   ],
                    [ 2. ,       9.9,  0.5-2.0j],
                    [ 6.5,       2.6,  5.7+1.0j]])
        A = csr_matrix( A )
        cases.append(A)
        P = diag_sparse([1,0,1.0j])
        cases.append( P*A*P )
        P = diag_sparse([0,1+1.0j,0])
        cases.append( P*A*P )

        for A in cases:
            D_sqrt,D_sqrt_inv,DAD = symmetric_rescaling(A)
            assert_almost_equal( diag_sparse(A) != 0, real(diag_sparse(DAD)) )
            assert_almost_equal( diag_sparse(DAD), D_sqrt*D_sqrt_inv )

            D_sqrt,D_sqrt_inv = diag_sparse(D_sqrt),diag_sparse(D_sqrt_inv)
            assert_almost_equal((D_sqrt_inv*A*D_sqrt_inv).todense(), DAD.todense())

    def test_get_diagonal(self):
        cases = []
        for i in range(1,6):
            A = rand(i,i)
            Ai = A + 1.0j*rand(i,i)
            cases.append(csr_matrix(A)) 
            cases.append(csr_matrix(Ai)) 


        for A in cases:
            D_A      = get_diagonal(A, norm_eq=False, inv=False)
            D_A_inv  = get_diagonal(A, norm_eq=False, inv=True)
            D_AA     = get_diagonal(A, norm_eq=True, inv=False)
            D_AA_inv = get_diagonal(A, norm_eq=True, inv=True)
            
            D = diag(A.todense())
            assert_almost_equal(D, D_A)
            D = 1.0/D
            assert_almost_equal(D, D_A_inv)
            
            D = diag((A*A.H).todense())
            assert_almost_equal(D, D_AA)
            D = 1.0/D
            assert_almost_equal(D, D_AA_inv)


    def test_profile_solver(self):
        from scipy.sparse.linalg import cg
        from pyamg.gallery import poisson
        from pyamg.aggregation import smoothed_aggregation_solver

        A = poisson((100,100), format='csr')
        A.data = A.data + 1e-5*rand(A.nnz)
        ml = smoothed_aggregation_solver(A)

        opts = []
        opts.append( {} )
        opts.append( {'accel' : cg } )
        opts.append( {'accel' : cg, 'tol' : 1e-10 } )

        for kwargs in opts:
            residuals = profile_solver(ml, **kwargs)


