from pyamg.testing import *

from numpy import matrix, array, diag, zeros, sqrt, abs, ravel, ones, arange, eye
from scipy import rand, linalg, real, imag, mat, diag, isscalar, ones, hstack
from scipy.sparse import csr_matrix, isspmatrix, bsr_matrix, isspmatrix_bsr, spdiags

import pyamg
from pyamg.util.utils import *

class TestUtils(TestCase):
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
    
    def test_symmetric_rescaling_sa(self):
        cases = []
        # case 1
        e = ones((5,1)).ravel()
        data = [ -1*e, 2*e, -1*e ]
        A = spdiags(data,[-1,0,1],5,5).tocsr()
        B = e.copy().reshape(-1,1)
        DAD_answer = array([[ 1. ,-0.5, 0. , 0. , 0. ],
                            [-0.5, 1. ,-0.5, 0. , 0. ],
                            [ 0. ,-0.5, 1. ,-0.5, 0. ],
                            [ 0. , 0. ,-0.5, 1. ,-0.5],
                            [ 0. , 0. , 0. ,-0.5, 1. ]])
        DB_answer = sqrt(2*e.reshape(-1,1))
        #            matrix   B    BH   expected matrix  expected B  expected BH 
        cases.append( (A,     B,  None,   DAD_answer,      DB_answer,   None) )

        # case 2
        A2 = A.copy()
        A2.symmetry = 'nonsymmetric'
        cases.append( (A2, B.copy(), B.copy(), DAD_answer, DB_answer, DB_answer) )
        
        # case 3
        A3 = A.copy()
        A3.symmetry = 'hermitian'
        cases.append( (A3, B.copy(), None, DAD_answer, DB_answer, None) )

        # case 4
        B4 = hstack( (B.copy(), 2*B.copy()) )
        DB4_answer = sqrt(2)*B4
        A4 = A.copy()
        A4.symmetry = 'nonsymmetric'
        cases.append( (A4, B4, B4.copy(), DAD_answer, DB4_answer, DB4_answer) )
        

        for case in cases:
            [A,B,BH,DAD_answer,DB_answer,DBH_answer] = case
            
            [DAD, DB, DBH] = symmetric_rescaling_sa(A,B,BH=BH)
            assert_array_almost_equal(DAD.todense(), DAD_answer)
            assert_array_almost_equal(DB, DB_answer)
            if DBH_answer != None:
                assert_array_almost_equal(DBH, DBH_answer)


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

    def test_get_block_diag(self):
        from scipy import arange, ravel, array
        from scipy.sparse import csr_matrix
        A = csr_matrix(arange(1,37, dtype=float).reshape(6,6))
        block_diag = get_block_diag(A, blocksize=1, inv_flag=False)
        assert_array_almost_equal(ravel(block_diag), A.diagonal())

        block_diag = get_block_diag(A, blocksize=2, inv_flag=False)
        answer = array([[[  1.,   2.],
                         [  7.,   8.]],
                        [[ 15.,  16.],
                         [ 21.,  22.]],
                        [[ 29.,  30.],
                         [ 35.,  36.]]])
        assert_array_almost_equal(ravel(block_diag), ravel(answer))

        block_diag_inv = get_block_diag(A, blocksize=2, inv_flag=True)
        answer = array([[[-1.33333333,  0.33333333],
                         [ 1.16666667, -0.16666667]],
                        [[-3.66666667,  2.66666667],
                         [ 3.5       , -2.5       ]],
                        [[-6.        ,  5.        ],
                         [ 5.83333333, -4.83333333]]])
        assert_array_almost_equal(ravel(block_diag_inv), ravel(answer), decimal=3)

        # try with singular (1,1) block, a zero (2,2) block and a zero (0,2) block
        A = bsr_matrix(array([ 1.,  2.,  3.,  4.,  0.,  0.,
                               5.,  6.,  7.,  8.,  0.,  0.,
                               9., 10., 11., 11., 12., 13.,
                              14., 15., 16., 16., 18., 19.,
                              20., 21., 22., 23.,  0.,  0.,
                              26., 27., 28., 29.,  0.,  0.,]).reshape(6,6), blocksize=(3,3))
        block_diag_inv = get_block_diag(A, blocksize=2, inv_flag=True)
        answer = array([[[-1.5       ,  0.5       ],
                         [ 1.25      , -0.25      ]],
                        [[ 0.01458886,  0.02122016],
                         [ 0.01458886,  0.02122016]],
                        [[ 0.        ,  0.        ],
                         [ 0.        ,  0.        ]]])
        assert_array_almost_equal(ravel(block_diag_inv), ravel(answer), decimal=3)

        # try with different types of zero blocks
        A = bsr_matrix(array([ 0.,  0.,  3.,  4.,  0.,  0.,
                               0.,  0.,  7.,  8.,  0.,  0.,
                               0.,  0.,  0.,  0.,  0.,  0.,
                               0.,  0.,  0.,  0.,  0.,  0.,
                               0.,  0.,  0.,  0., 22., 23.,
                               0.,  0.,  0.,  0., 28., 29.,]).reshape(6,6), blocksize=(2,2))
        block_diag_inv = get_block_diag(A, blocksize=2, inv_flag=False)
        answer = array([[[ 0.        ,  0.        ],
                         [ 0.        ,  0.        ]],
                        [[ 0.        ,  0.        ],
                         [ 0.        ,  0.        ]],
                        [[22.        , 23.        ],
                         [28.        , 29.        ]]])
        assert_array_almost_equal(ravel(block_diag_inv), ravel(answer), decimal=3)


    def test_relaxation_as_linear_operator(self):
        As = []
        bs =[]
        xs = []
        methods = ['gauss_seidel', 'jacobi', 'block_gauss_seidel', 'block_jacobi']
        params = [{}, {'iterations' : 2}]
        As.append(pyamg.gallery.poisson( (10,10), format='csr'))
        As.append(1.0j*pyamg.gallery.poisson( (10,10), format='csr'))
        As.append(1.0j*pyamg.gallery.elasticity.linear_elasticity( (20,20) )[0] )
        As.append(pyamg.gallery.elasticity.linear_elasticity( (20,20) )[0] )
        for A in As:
            if A.dtype == 'complex':
                xs.append(rand(A.shape[0],1)+1.0j*rand(A.shape[0],1))
                bs.append(rand(A.shape[0],1)+1.0j*rand(A.shape[0],1))
            else:
                bs.append(rand(A.shape[0],1))
                xs.append(rand(A.shape[0],1))

        for method in methods:
            for kwargs in params:    
                for (A,x,b) in zip(As,xs,bs):
                    kwargs_linop = dict(kwargs)
                    ##
                    # run relaxation as a linear operator
                    if kwargs_linop == dict({}):
                        relax = relaxation_as_linear_operator(method, A, b)
                    else:
                        relax = relaxation_as_linear_operator((method,kwargs_linop), A, b)
                    x_linop = relax*x
                    
                    ##
                    # manually run the relaxation routine
                    relax2 = eval('pyamg.relaxation.' + method)
                    x_gold = x.copy()
                    blockflag = False
                    kwargs_gold = dict(kwargs)
                    # deal with block matrices
                    if method.startswith('block') and isspmatrix_bsr(A):
                        blockflag = True
                        kwargs_gold['blocksize'] = A.blocksize[0]
                    # deal with omega and jacobi
                    # --> note that we assume the default setup for jacobi uses omega = 1/rho
                    if method.endswith('jacobi'):
                        if blockflag:
                            kwargs_gold['omega'] = 1.0/A.rho_block_D_inv
                        else:
                            kwargs_gold['omega'] = 1.0/A.rho_D_inv

                    relax2(A,x_gold,b,**kwargs_gold)
                    
                    assert_array_almost_equal(x_linop, x_gold)

    def test_filter_operator(self):
        ##
        # Basic tests of dimension 1 and 2 problems
        # 1x1
        A = csr_matrix(array([[1.2]]))
        C = csr_matrix(array([[1.]]))
        B = array([[0.5]])
        Bf= array([[1.5]])
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = matrix([[3.0]])
        assert_array_almost_equal(A_known, A_filter)
        # 1x1, but with no entries in C
        C = csr_matrix(array([[0.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[0.0]])
        assert_array_almost_equal(A_known, A_filter)
        # 1x1, but with no entries in A
        A = csr_matrix(array([[0.]]))
        C = csr_matrix(array([[1.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[3.0]])
        assert_array_almost_equal(A_known, A_filter)


        # 1x2
        A = csr_matrix(array([[1.2, 1.]]))
        C = csr_matrix(array([[1.,  1.]]))
        B = array([[0.5], [0.5]])
        Bf= array([[1.5]])
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = matrix([[ 1.6,  1.4]])
        assert_array_almost_equal(A_known, A_filter)
        # 1x2, but sparser
        C = csr_matrix(array([[0.,  1.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[0., 3.]])
        assert_array_almost_equal(A_known, A_filter)
        # 1x2, but with no entries
        C = csr_matrix(array([[0.,  0.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[0., 0.]])
        assert_array_almost_equal(A_known, A_filter)

        # 2x1
        A = csr_matrix(array([[1.2], [1.]]))
        C = csr_matrix(array([[1.],  [1.]]))
        B = array([[0.5]])
        Bf= array([[1.5], [0.4]])
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = matrix([[ 3.],  [0.8]])
        assert_array_almost_equal(A_known, A_filter)
        # 2x1, but sparser
        C = csr_matrix(array([[0.],  [1.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[0.], [.8]])
        assert_array_almost_equal(A_known, A_filter)
        # 2x1, but with no entries
        C = csr_matrix(array([[0.],  [0.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[0.], [0.]])
        assert_array_almost_equal(A_known, A_filter)

        # 2x2
        A = csr_matrix(array([[1.2, 1.1], [1., 0.5]]))
        C = csr_matrix(array([[1.2, 1.1], [1., 0.]]))
        B = array([[0.5, 1.0],[0.5, 1.1]])
        Bf = array([[0.5,1.0],[0.5,1.1]])
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[1., 0.], [1.08, 0.]])
        assert_array_almost_equal(A_known, A_filter)
        # 1x2, but sparser
        C = csr_matrix(array([[0., 0.], [1., 0.]]))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[0., 0.], [1.08, 0.]])
        assert_array_almost_equal(A_known, A_filter)
        # Try block structure
        A = A.tobsr((2,2))
        C = C.tobsr((2,2))
        A_filter = filter_operator(A, C, B, Bf).todense()
        A_known = array([[1., 0.], [0., 1.]])
        assert_array_almost_equal(A_known, A_filter)
        
        ##
        # Basic tests, with easy to compute answers
        # test one, the constant
        A = array([ [1.,1,1],[1,1,1],[0,1,0],[0,1,0],[0,0,1],[0,0,1]])
        C = array([ [1.,1,0],[1,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1]])
        B = ones((3,1))
        Bf = ones((6,1))
        A_filter = filter_operator(csr_matrix(A), csr_matrix(C), B, Bf).todense()
        A_known = matrix([[ 0.5,  0.5,  0. ],
                          [ 0.5,  0.5,  0. ],
                          [ 0. ,  1. ,  0. ],
                          [ 0. ,  1. ,  0. ],
                          [ 0. ,  0. ,  1. ],
                          [ 0. ,  0. ,  1. ]])
        assert_array_almost_equal(A_known, A_filter)
        # test two, the constant and linears
        B = hstack( (B, arange(B.shape[0]).reshape(-1,1)) )
        Bf = hstack( (Bf, arange(Bf.shape[0]).reshape(-1,1)) )
        A_filter = filter_operator(csr_matrix(A), csr_matrix(C), B, Bf).todense()
        A_known = matrix([[ 1. ,  0. ,  0. ],
                          [ 0. ,  1. ,  0. ],
                          [ 0. ,  1.5,  0. ],
                          [ 0. ,  2. ,  0. ],
                          [ 0. ,  0. ,  1.8],
                          [ 0. ,  0. ,  2.2]])
        assert_array_almost_equal(A_known, A_filter)
        
        ##
        # Run two tests based on the Laplacian
        # first test, constants
        from pyamg.gallery import poisson
        A = poisson((10,10), format='csr')
        C = A.copy()
        C.data[arange(0, C.nnz, 5) ] = 0.0
        C.eliminate_zeros()
        B = ones((A.shape[0],1))
        Bf = ones((A.shape[0],1))
        A_filter = filter_operator(A, C, B, Bf)
        assert_array_almost_equal(A_filter*B, Bf)
        # second test, constants and linears
        B = hstack( (B, arange(B.shape[0]).reshape(-1,1)) )
        Bf = hstack( (Bf, arange(Bf.shape[0]).reshape(-1,1)) )
        A_filter = filter_operator(A, C, B, Bf)
        assert_array_almost_equal(A_filter*B, Bf)

    def test_scale_T(self):
        from scipy.sparse import csr_matrix, bsr_matrix
        from scipy import matrix

        ##
        # Trivially sized tests
        # 1x1
        T   = matrix([[ 1.1 ]] )
        P_I = matrix([[ 1.0 ]] )
        I_F = matrix([[ 0.0 ]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 1.0 ]] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 1.1 ]] )
        P_I = matrix([[ 0.0 ]] )
        I_F = matrix([[ 1.0 ]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 1.1 ]] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 0.0 ]] )
        P_I = matrix([[ 0.0 ]] )
        I_F = matrix([[ 1.0 ]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 0. ]] )
        assert_array_almost_equal(T_answer, T_scaled)

        # 2x1
        T   = matrix([[ 1.5 ], [1.2] ] )
        P_I = matrix([[ 1.  ], [0. ] ] )
        I_F = matrix([[ 0., 0. ], [0., 1.]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 1. ], [0.8] ] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 0.  ], [1.2] ] )
        P_I = matrix([[ 1.  ], [0. ] ] )
        I_F = matrix([[ 0., 0. ], [0., 1.]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 1. ], [0. ] ] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 0.  ], [0. ] ] )
        P_I = matrix([[ 1.  ], [0. ] ] )
        I_F = matrix([[ 0., 0. ], [0., 1.]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 1. ], [0. ] ] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 0.  ], [0. ] ] )
        P_I = matrix([[ 0.  ], [0. ] ] )
        I_F = matrix([[ 1., 0. ], [0., 1.]] )
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        T_answer = matrix([[ 0. ], [0. ] ] )
        assert_array_almost_equal(T_answer, T_scaled)

        # 2x2
        T   = matrix([[ 2., 0. ], [1., 1. ] ] )
        P_I = matrix([[ 1., 0. ], [0., 1. ] ] )
        I_F = matrix([[ 0., 0. ], [0., 0. ]] )
        T_scaled = scale_T(bsr_matrix(T,  blocksize=(1,1)), 
                           bsr_matrix(P_I,blocksize=(1,1)), 
                           bsr_matrix(I_F,blocksize=(1,1))).todense()
        T_answer = matrix([[ 1., 0. ], [0., 1. ] ] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 2., 0. ], [1., 1. ] ] )
        P_I = matrix([[ 1., 0. ], [0., 1. ] ] )
        I_F = matrix([[ 0., 0. ], [0., 0. ]] )
        T_scaled = scale_T(bsr_matrix(T,  blocksize=(2,2)), 
                           bsr_matrix(P_I,blocksize=(2,2)), 
                           bsr_matrix(I_F,blocksize=(2,2))).todense()
        T_answer = matrix([[ 1., 0. ], [0., 1. ] ] )
        assert_array_almost_equal(T_answer, T_scaled)
        ##
        T   = matrix([[ 2., 0. ], [1., 1. ] ] )
        P_I = matrix([[ 0., 0. ], [0., 0. ] ] )
        I_F = matrix([[ 1., 0. ], [0., 1. ]] )
        T_scaled = scale_T(bsr_matrix(T,  blocksize=(2,2)), 
                           bsr_matrix(P_I,blocksize=(2,2)), 
                           bsr_matrix(I_F,blocksize=(2,2))).todense()
        T_answer = matrix([[ 2., 0. ], [1., 1. ] ] )
        assert_array_almost_equal(T_answer, T_scaled)

        ##
        # Test for one CSR and one BSR example
        T = matrix([[ 1.0,  0.,   0. ],
                    [ 0.5,  0.,   0. ],
                    [ 0. ,  1.,   0. ],
                    [ 0. ,  0.5,  0. ],
                    [ 0. ,  0.,   1. ],
                    [ 0. ,  0.,   0.25 ]])
        P_I = matrix([[ 0.,  0.,   0. ],
                      [ 1.,  0.,   0. ],
                      [ 0.,  1.,   0. ],
                      [ 0.,  0.,   0. ],
                      [ 0.,  0.,   0. ],
                      [ 0.,  0.,   1. ]])
        I_F = matrix([[ 1.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  1.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  1.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.]])
        T_answer = matrix([[ 2. ,  0. ,  0. ],
                           [ 1. ,  0. ,  0. ],
                           [ 0. ,  1. ,  0. ],
                           [ 0. ,  0.5,  0. ],
                           [ 0. ,  0. ,  4. ],
                           [ 0. ,  0. ,  1. ]])
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        assert_array_almost_equal(T_answer, T_scaled)

        ##
        # BSR test
        T = matrix([[ 1.0, 1., 0.,  0. ],
                    [ 0.5, 1., 0.,  0. ],
                    [ 1. , 0., 0.,  0. ],
                    [ 0. , 1., 0.,  0. ],
                    [ 0. , 0., 2.,  1. ],
                    [ 0. , 0., 3.,  1. ],
                    [ 0. , 0., 4.,  1. ],
                    [ 0. , 0., 2.,  0. ]])
        P_I = matrix([[ 0.,  0.,  0.,   0.],
                      [ 0.,  0.,  0.,   0.],
                      [ 1.,  0.,  0.,   0.],
                      [ 0.,  1.,  0.,   0.],
                      [ 0.,  0.,  1.,   0.],
                      [ 0.,  0.,  0.,   1.],
                      [ 0.,  0.,  0.,   0.],
                      [ 0.,  0.,  0.,   0.]])
        I_F = matrix([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.]])
        T_answer = matrix([[ 1. , 1.,  0.,  0. ],
                           [ 0.5, 1.,  0.,  0. ],
                           [ 1. , 0.,  0.,  0. ],
                           [ 0. , 1.,  0.,  0. ],
                           [ 0. , 0.,  1.,  0. ],
                           [ 0. , 0.,  0.,  1. ],
                           [ 0. , 0., -1.,  2. ],
                           [ 0. , 0., -2.,  2. ]]) 
        T = bsr_matrix(T, blocksize=(2,2))
        P_I = bsr_matrix(P_I, blocksize=(2,2))
        I_F = bsr_matrix(I_F, blocksize=(2,2))
        T_scaled = scale_T(T, P_I, I_F).todense()
        assert_array_almost_equal(T_answer, T_scaled)

        ##
        # BSR test
        T = matrix([[ 1.0, 1., 0.,  0. ],
                    [ 0.5, 1., 0.,  0. ],
                    [ 1. , 1., 0.,  0. ],
                    [ 1. , 1., 0.,  0. ],
                    [ 0. , 0., 2.,  1. ],
                    [ 0. , 0., 3.,  1. ],
                    [ 0. , 0., 4.,  1. ],
                    [ 0. , 0., 2.,  0. ]])
        P_I = matrix([[ 0.,  0.,  0.,   0.],
                      [ 0.,  0.,  0.,   0.],
                      [ 1.,  0.,  0.,   0.],
                      [ 0.,  1.,  0.,   0.],
                      [ 0.,  0.,  1.,   0.],
                      [ 0.,  0.,  0.,   1.],
                      [ 0.,  0.,  0.,   0.],
                      [ 0.,  0.,  0.,   0.]])
        I_F = matrix([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.]])
        T_answer = matrix([[ 0.5  , 0.5  , 0.,  0. ],
                           [ 0.375, 0.375, 0.,  0. ],
                           [ 1.   , 0.   , 0.,  0. ],
                           [ 0.   , 1.   , 0.,  0. ],
                           [ 0.   , 0.   , 1.,  0. ],
                           [ 0.   , 0.   , 0.,  1. ],
                           [ 0.   , 0.   ,-1.,  2. ],
                           [ 0.   , 0.   ,-2.,  2. ]]) 
        Cpts = array([1,2])
        T = bsr_matrix(T, blocksize=(2,2))
        P_I = bsr_matrix(P_I, blocksize=(2,2))
        I_F = bsr_matrix(I_F, blocksize=(2,2))
        T_scaled = scale_T(T, P_I, I_F).todense()
        assert_array_almost_equal(T_answer, T_scaled)
    
    def test_get_Cpt_params(self):
        from pyamg.gallery import poisson
        from scipy.sparse import csr_matrix, bsr_matrix

        ##
        # Begin with trivially sized tests
        # 1x1
        A = csr_matrix(array([[1.2]]))
        Cpts = array([0])
        AggOp = csr_matrix(array([[1. ]]))
        T = AggOp.copy().tobsr()
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[ 1.]]), blocksize=(1,1))
        I_F = bsr_matrix( array([[ 0.]]), blocksize=(1,1))
        P_I = bsr_matrix( array([[ 1.]]), blocksize=(1,1) )
        assert_equal(array([0]),  params['Cpts'])
        assert_equal(array([ ]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)
        ##
        A = csr_matrix(array([[1.2]]))
        Cpts = array([ ])
        AggOp = csr_matrix(array([[1. ]]))
        T = AggOp.copy().tobsr()
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[ 0.]]), blocksize=(1,1))
        I_F = bsr_matrix( array([[ 1.]]), blocksize=(1,1))
        P_I = bsr_matrix( array([[ 0.]]), blocksize=(1,1) )
        assert_equal(array([ ]),  params['Cpts'])
        assert_equal(array([0]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)
        ##
        # 2x2
        A = csr_matrix(array([[1., 1.],[1., 1.]]))
        Cpts = array([0])
        AggOp = csr_matrix(array([[1.], [1.]]))
        T = AggOp.copy().tobsr()
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[1., 0.],[0., 0.]]), blocksize=(1,1))
        I_F = bsr_matrix( array([[0., 0.],[0., 1.]]), blocksize=(1,1))
        P_I = bsr_matrix( array([[1.], [0.]]), blocksize=(1,1) )
        assert_equal(array([0]),  params['Cpts'])
        assert_equal(array([1]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)
        ##
        Cpts = array([0,1])
        AggOp = csr_matrix(array([[1.,0], [0.,1.]]))
        T = AggOp.copy().tobsr()
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[1., 0.],[0., 1.]]), blocksize=(1,1))
        I_F = bsr_matrix( array([[0., 0.],[0., 0.]]), blocksize=(1,1))
        P_I = bsr_matrix( array([[1., 0.], [0., 1.]]), blocksize=(1,1) )
        assert_equal(array([0,1]),  params['Cpts'])
        assert_equal(array([ ]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)
        ##
        Cpts = array([ ])
        AggOp = csr_matrix(array([[0.], [0.]]))
        T = AggOp.copy().tobsr()
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[0., 0.],[0., 0.]]), blocksize=(1,1))
        I_F = bsr_matrix( array([[1., 0.],[0., 1.]]), blocksize=(1,1))
        P_I = bsr_matrix( array([[ 0.], [0. ]]), blocksize=(1,1) )
        assert_equal(array([ ]),  params['Cpts'])
        assert_equal(array([0,1]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)
        ##
        A = A.tobsr( blocksize=(2,2) )
        Cpts = array([0])
        AggOp = csr_matrix(array([[1.]]) )
        T = bsr_matrix(array([[1., 1.], [1., 2.]]), blocksize=(2,2))
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[1., 0.],[0., 1.]]), blocksize=(2,2))
        I_F = bsr_matrix( array([[0., 0.],[0., 0.]]), blocksize=(2,2))
        P_I = bsr_matrix( array([[1., 0.],[0., 1.]]), blocksize=(2,2))
        assert_equal(array([0,1]),  params['Cpts'])
        assert_equal(array([ ]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)
        ##
        Cpts = array([ ])
        AggOp = csr_matrix(array([[1.]]) )
        T = bsr_matrix(array([[1., 1.], [1., 2.]]), blocksize=(2,2))
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix( array([[0., 0.],[0., 0.]]), blocksize=(2,2))
        I_F = bsr_matrix( array([[1., 0.],[0., 1.]]), blocksize=(2,2))
        P_I = bsr_matrix( array([[0., 0.],[0., 0.]]), blocksize=(2,2))
        assert_equal(array([ ]),  params['Cpts'])
        assert_equal(array([0,1]),  params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)

        ##
        # Begin more "realistic" tests
        A = poisson((10,), format='csr')
        Cpts = array([3, 7])
        AggOp = ([[ 1., 0.],
                  [ 1., 0.],
                  [ 1., 0.],
                  [ 1., 0.],
                  [ 1., 0.],
                  [ 0., 1.],
                  [ 0., 1.],
                  [ 0., 1.],
                  [ 0., 1.],
                  [ 0., 1.]])
        AggOp = csr_matrix(AggOp)
        T = AggOp.copy().tobsr()

        ##
        # CSR Test
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix(( array([[[ 1.]], [[ 1.]]]), 
                array([3, 7]), array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2]) ), 
                shape=(10,10) )
        I_F = bsr_matrix(( 
              array([[[ 1.]], [[ 1.]], [[ 1.]], [[ 1.]], [[ 1.]], [[ 1.]], [[ 1.]], [[ 1.]]]), 
              array([0, 1, 2, 4, 5, 6, 8, 9]), 
              array([0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8]) ), 
              shape=(10,10) )
        P_I = matrix([[ 0.,  0.],
                      [ 0.,  0.],
                      [ 0.,  0.],
                      [ 1.,  0.],
                      [ 0.,  0.],
                      [ 0.,  0.],
                      [ 0.,  0.],
                      [ 0.,  1.],
                      [ 0.,  0.],
                      [ 0.,  0.]]) 
        P_I = bsr_matrix(P_I, blocksize=(1,1))
        Fpts = array([0,1,2,4,5,6,8,9])
        assert_equal(Cpts, params['Cpts'])
        assert_equal(Fpts, params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)

        ##
        # BSR Test
        A = A.tobsr(blocksize=(2,2)) 
        Cpts = array([1, 3])
        AggOp = ([[ 1., 0.],
                  [ 1., 0.],
                  [ 1., 0.],
                  [ 0., 1.],
                  [ 0., 1.]])
        AggOp = csr_matrix(AggOp)
        T = hstack((T.todense(), T.todense()))[:,[0,2,1,3]]
        T = bsr_matrix(T, blocksize=(2,2)) 
        params = get_Cpt_params(A, Cpts, AggOp, T)
        I_C = bsr_matrix(( array([ [[ 1.,  0.],[ 0.,  1.]],
                                   [[ 1.,  0.],[ 0.,  1.]]]), 
                array([1, 3]), 
                array([0, 0, 1, 1, 2, 2]) ), 
                shape=(10,10) )
        I_F = bsr_matrix(( 
              array([[[ 1.,  0.],[ 0.,  1.]], 
                     [[ 1.,  0.],[ 0.,  1.]], 
                     [[ 1.,  0.],[ 0.,  1.]]]), 
              array([0, 2, 4]), 
              array([0, 1, 1, 2, 2, 3]) ), 
              shape=(10,10) )
        P_I = matrix([[ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.],
                      [ 1.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  1.,  0.],
                      [ 0.,  0.,  0.,  1.],
                      [ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.]]) 
        P_I = bsr_matrix(P_I, blocksize=(2,2))
        Fpts = array([0, 1, 4, 5, 8, 9])
        Cpts = array([2, 3, 6, 7])
        assert_equal(Cpts, params['Cpts'])
        assert_equal(Fpts, params['Fpts'])
        assert_equal(I_C.indptr,  params['I_C'].indptr)
        assert_equal(I_C.indices, params['I_C'].indices)
        assert_equal(I_C.data,    params['I_C'].data)
        assert_equal(I_F.indptr,  params['I_F'].indptr)
        assert_equal(I_F.indices, params['I_F'].indices)
        assert_equal(I_F.data,    params['I_F'].data)
        assert_equal(P_I.indptr,  params['P_I'].indptr)
        assert_equal(P_I.indices, params['P_I'].indices)
        assert_equal(P_I.data,    params['P_I'].data)

    def test_compute_BtBinv(self):
        ##
        # Trivially sized tests
        # 1x1x1
        T = matrix([[ 1.]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[1.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([ [[ 1. ]] ])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[0.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([ [[ 0. ]] ])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[0.5]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([ [[ 4. ]] ])
        assert_array_almost_equal(BtBinv, answer)
        ##
        # 2x1x1
        T = matrix([[ 1.,0.], [1.,1.]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[1.], [1.] ])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 1. ]], [[ 0.5]]])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.,0.], [1.,1.]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[0.], [1.] ])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 0. ]], [[ 1.]]])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.,0.], [1.,1.]])
        T = bsr_matrix(T, blocksize=(2,2))
        B = array([[0.], [2.] ])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 0.25 ]]])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.,  0.], [ 1.,  0.],
                    [ 0.,  .5], [ 0.,  .25]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[1.],[2.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 1.  ]], [[ 1.  ]],
                        [[ 0.25]], [[ 0.25]]])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.,  0.], [ 0.,  .25]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[1., 1.],[2., 1.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 0.25,  0.25], [ 0.25,  0.25]],
                        [[ 0.16,  0.08], [ 0.08,  0.04]]])
        assert_array_almost_equal(BtBinv, answer)
        ##
        T = matrix([[ 1.,  0.], [ 0.,  .25]])
        T = bsr_matrix(T, blocksize=(2,2))
        B = array([[1., 1.],[1., 1.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 0.125,  0.125],
                         [ 0.125,  0.125]]])
        assert_array_almost_equal(BtBinv, answer)

        ##
        # Simple BSR test
        T = matrix([[ 1.  ,  1.  ,  0.  ,  0.  ],
                    [ 1.  ,  1.  ,  0.  ,  0.  ],
                    [ 0.  ,  0.  ,  0.5 ,  0.5 ],
                    [ 0.  ,  0.  ,  0.25,  0.25]])
        T = bsr_matrix(T, blocksize=(2,2))
        B = array([[1., 1.],[1., 2.],[1., 1.],[1., 3.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 5. , -3. ], [-3. ,  2. ]],
                        [[ 2.5, -1. ], [-1. ,  0.5]]])
        assert_array_almost_equal(BtBinv, answer)

class TestComplexUtils(TestCase):
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

    def test_symmetric_rescaling_sa(self):
        cases = []
        # case 1
        e = 1.0j*ones((5,1)).ravel()
        data = [ -1*e, 2*e, -1*e ]
        A = 1.0j*spdiags(data,[-1,0,1],5,5).tocsr()
        B = e.copy().reshape(-1,1)
        DAD_answer = array([[ 1. ,-0.5, 0. , 0. , 0. ],
                            [-0.5, 1. ,-0.5, 0. , 0. ],
                            [ 0. ,-0.5, 1. ,-0.5, 0. ],
                            [ 0. , 0. ,-0.5, 1. ,-0.5],
                            [ 0. , 0. , 0. ,-0.5, 1. ]])
        DB_answer = sqrt(2)*1.0j*e.reshape(-1,1)
        #            matrix   B    BH   expected matrix  expected B  expected BH 
        cases.append( (A,     B,  None,   DAD_answer,      DB_answer,   None) )

        for case in cases:
            [A,B,BH,DAD_answer,DB_answer,DBH_answer] = case
            
            [DAD, DB, DBH] = symmetric_rescaling_sa(A,B,BH=BH)
            assert_array_almost_equal(DAD.todense(), DAD_answer)
            assert_array_almost_equal(DB, DB_answer)
            if DBH_answer != None:
                assert_array_almost_equal(DBH, DBH_answer)


    def test_get_diagonal(self):
        cases = []
        for i in range(1,6):
            A = rand(i,i)
            Ai = A + 1.0j*rand(i,i)
            cases.append(csr_matrix(A)) 
            cases.append(csr_matrix(Ai)) 


        for A in cases:
            D_A       = get_diagonal(A, norm_eq=False, inv=False)
            D_A_inv   = get_diagonal(A, norm_eq=False, inv=True)
            D_AA      = get_diagonal(A, norm_eq=1, inv=False)
            D_AA_inv  = get_diagonal(A, norm_eq=1, inv=True)
            D_AA2     = get_diagonal(A, norm_eq=2, inv=False)
            D_AA_inv2 = get_diagonal(A, norm_eq=2, inv=True)
            
            D = diag(A.todense())
            assert_almost_equal(D, D_A)
            D = 1.0/D
            assert_almost_equal(D, D_A_inv)
            
            D = diag((A.H*A).todense())
            assert_almost_equal(D, D_AA)
            D = 1.0/D
            assert_almost_equal(D, D_AA_inv)
            
            D = diag((A*A.H).todense())
            assert_almost_equal(D, D_AA2)
            D = 1.0/D
            assert_almost_equal(D, D_AA_inv2)


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

    def test_to_type(self):
        w = 1.2
        x = ones((5,1))
        y = rand(3,2)
        z = csr_matrix(rand(2,2))
        inlist = [w, x, y, z]

        out = to_type(complex, inlist)
        for i in range(len(out)):
            assert( out[i].dtype==complex ) 
            if isspmatrix(out[i]):
                diff = ravel(out[i].data - inlist[i].data)
            else:
                diff = out[i] - inlist[i]
            assert_equal( max(abs(ravel(diff))), 0.0)

    def test_type_prep(self):
        w = 1.2
        x = ones((5,1))
        y = rand(3,2)
        z = csr_matrix(rand(2,2))
        inlist = [w, x, y, z]

        out = type_prep(complex, inlist)
        for i in range(len(out)):
            assert( out[i].dtype==complex ) 
            assert( not isscalar(out[i]) )
            if isspmatrix(out[i]):
                diff = ravel(out[i].data - inlist[i].data)
            else:
                diff = out[i] - inlist[i]
            assert_equal( max(abs(ravel(diff))), 0.0)
    
    def test_filter_operator(self):

        ##
        # Basic tests, with easy to compute answers
        # test one, the constant
        A = array([ [1.+0.j,1,1],[1,1,1],[0,1,0],[0,1,0],[0,0,1],[0,0,1]])
        C = array([ [1.+0.j,1,0],[1,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1]])
        B = ones((3,1)) + 0.j
        Bf = ones((6,1)) + 1.0j * ones((6,1))
        A_filter = filter_operator(csr_matrix(A), csr_matrix(C), B, Bf).todense()
        A_known = matrix([[ 0.5+0.5j,  0.5+0.5j,  0.0+0.j ],
                          [ 0.5+0.5j,  0.5+0.5j,  0.0+0.j ],
                          [ 0.0+0.j ,  1.0+1.j ,  0.0+0.j ],
                          [ 0.0+0.j ,  1.0+1.j ,  0.0+0.j ],
                          [ 0.0+0.j ,  0.0+0.j ,  1.0+1.j ],
                          [ 0.0+0.j ,  0.0+0.j ,  1.0+1.j ]])
        assert_array_almost_equal(A_known, A_filter)
        
        ##
        # test two, the constant and linears
        # Note that for the rows with only one nonzero, Bf can't be
        # approximated exactly
        B = hstack( (B, arange(B.shape[0]).reshape(-1,1)) )
        Bf = hstack( (Bf, arange(Bf.shape[0]).reshape(-1,1) + 1.0j*arange(Bf.shape[0]).reshape(-1,1) ) )
        A_filter = filter_operator(csr_matrix(A), csr_matrix(C), B, Bf).todense()
        A_known = matrix([[ 1.0+1.j ,  0.0+0.j ,  0.0+0.j ],
                          [ 0.0+0.j ,  1.0+1.j ,  0.0+0.j ],
                          [ 0.0+0.j ,  1.5+1.5j,  0.0+0.j ],
                          [ 0.0+0.j ,  2.0+2.j ,  0.0+0.j ],
                          [ 0.0+0.j ,  0.0+0.j ,  1.8+1.8j],
                          [ 0.0+0.j ,  0.0+0.j ,  2.2+2.2j]])
        assert_array_almost_equal(A_known, A_filter)
        
    def test_scale_T(self):
        from scipy.sparse import csr_matrix, bsr_matrix
        from scipy import matrix

        ##
        # Test for one CSR and one BSR example
        T = matrix([[ 1.0,  0.,   0. ],
                    [ 0.5j, 0.,   0. ],
                    [ 0. ,  1.,   0. ],
                    [ 0. ,  .5j,  0. ],
                    [ 0. ,  0.,   1.j],
                    [ 0. ,  0.,   0.25 ]])
        T_answer = matrix([[ 0.-2.j ,  0.+0.j ,  0.+0.j ],
                           [ 1.+0.j ,  0.+0.j ,  0.+0.j ],
                           [ 0.+0.j ,  1.+0.j ,  0.+0.j ],
                           [ 0.+0.j ,  0.+0.5j,  0.+0.j ],
                           [ 0.+0.j ,  0.+0.j ,  0.+4.j ],
                           [ 0.+0.j ,  0.+0.j ,  1.+0.j ]])
        P_I = matrix([[ 0.,  0.,   0. ],
                      [ 1.,  0.,   0. ],
                      [ 0.,  1.,   0. ],
                      [ 0.,  0.,   0. ],
                      [ 0.,  0.,   0. ],
                      [ 0.,  0.,   1. ]])
        I_F = matrix([[ 1.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  1.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  1.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.]])
        T_scaled = scale_T(bsr_matrix(T), bsr_matrix(P_I), bsr_matrix(I_F)).todense()
        assert_array_almost_equal(T_answer, T_scaled)

        ##
        # BSR test
        T = matrix([[ 1.j, 1., 0.,  0. ],
                    [ 0.5, 1., 0.,  0. ],
                    [ 1. , 0., 0.,  0. ],
                    [ 0. , 1., 0.,  0. ],
                    [ 0. , 0., 2.j, 0. ],
                    [ 0. , 0., 0.,  1. ],
                    [ 0. , 0., 1.,  1. ],
                    [ 0. , 0., 1.,  1. ]])
        P_I = matrix([[ 0.,  0.,  0.,   0.],
                      [ 0.,  0.,  0.,   0.],
                      [ 1.,  0.,  0.,   0.],
                      [ 0.,  1.,  0.,   0.],
                      [ 0.,  0.,  1.,   0.],
                      [ 0.,  0.,  0.,   1.],
                      [ 0.,  0.,  0.,   0.],
                      [ 0.,  0.,  0.,   0.]])
        I_F = matrix([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.]])
        T_answer = matrix([[ 0.0+1.j ,  1.0+0.j ,  0.0+0.j ,  0.0+0.j ],
                           [ 0.5+0.j ,  1.0+0.j ,  0.0+0.j ,  0.0+0.j ],
                           [ 1.0+0.j ,  0.0+0.j ,  0.0+0.j ,  0.0+0.j ],
                           [ 0.0+0.j ,  1.0+0.j ,  0.0+0.j ,  0.0+0.j ],
                           [ 0.0+0.j ,  0.0+0.j ,  1.0+0.j ,  0.0+0.j ],
                           [ 0.0+0.j ,  0.0+0.j ,  0.0+0.j ,  1.0+0.j ],
                           [ 0.0+0.j ,  0.0+0.j ,  0.0-0.5j,  1.0+0.j ],
                           [ 0.0+0.j ,  0.0+0.j ,  0.0-0.5j,  1.0+0.j ]]) 
        T = bsr_matrix(T, blocksize=(2,2))
        P_I = bsr_matrix(P_I, blocksize=(2,2))
        I_F = bsr_matrix(I_F, blocksize=(2,2))
        T_scaled = scale_T(T, P_I, I_F).todense()
        assert_array_almost_equal(T_answer, T_scaled)
    
    def test_compute_BtBinv(self):
        from pyamg.gallery import poisson
        
        ##
        # Simple CSR test
        T = matrix([[ 1.j,  0.], [ 1.,  0.],
                    [ 0.,  .5], [ 0.,  .25]])
        T = bsr_matrix(T, blocksize=(1,1))
        B = array([[1.+1.j],[2.j]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 0.50+0.j]], [[ 0.50+0.j]],
                        [[ 0.25+0.j]], [[ 0.25+0.j]]])
        assert_array_almost_equal(BtBinv, answer)

        ##
        # Simple BSR test
        T = matrix([[ 1.  ,  0.  ,  0.   ,  1. ],
                    [ 1.  ,  0.  ,  0.   ,  1. ],
                    [ 0.  ,  0.  ,  0.5  ,  0. ],
                    [ 0.  ,  0.  ,  0.25 ,  0. ]])
        T = bsr_matrix(T, blocksize=(2,2))
        B = array([[1.j, 1.],[1.j, 3.],[1.j, 4.],[1.j, 2.]])
        BtBinv = compute_BtBinv(B, T)
        answer = array([[[ 1.5+0.j ,  0.0+0.5j], [ 0.0-0.5j,  0.2+0.j ]],
                        [[ 5.0+0.j ,  0.0+1.5j], [ 0.0-1.5j,  0.5+0.j ]]])
        assert_array_almost_equal(BtBinv, answer)

    ## JBS: no explicitly complex tests necessary for get_Cpt_params
    ## def test_get_Cpt_params(self):

