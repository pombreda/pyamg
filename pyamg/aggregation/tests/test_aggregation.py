from pyamg.testing import *

import numpy
import scipy.sparse
from numpy import sqrt, ones, arange, array
from scipy import rand
from scipy.sparse import csr_matrix

from pyamg.utils import diag_sparse
from pyamg.gallery import poisson, linear_elasticity

from pyamg.aggregation.aggregation import smoothed_aggregation_solver

class TestParameters(TestCase):
    def setUp(self):
        self.cases = []

        self.cases.append(( poisson( (100,),  format='csr'), None))
        self.cases.append(( poisson( (10,10), format='csr'), None))
        self.cases.append( linear_elasticity( (10,10), format='bsr') )

    def run_cases(self, opts):
        for A,B in self.cases:
            ml = smoothed_aggregation_solver(A, B, max_coarse=5, **opts)

            numpy.random.seed(0) #make tests repeatable

            x = rand(A.shape[0])
            b = A*rand(A.shape[0])

            x_sol,residuals = ml.solve(b, x0=x, maxiter=30, tol=1e-10, return_residuals=True)
            convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            assert(convergence_ratio < 0.9)


    def test_strength_of_connection(self): 
        for strength in ['symmetric','ode']:
            self.run_cases( {'strength' : strength} )
    
    def test_aggregation_method(self): 
        for aggregate in ['standard','lloyd']:
            self.run_cases( {'aggregate' : aggregate} )
    
    def test_prolongation_smoother(self): 
        for smooth in ['jacobi','richardson','energy']:
            self.run_cases( {'smooth' : smooth} )

    def test_smoothers(self): 
        smoothers = []
        smoothers.append('gauss_seidel')
        #smoothers.append( ('sor',{'omega':0.9}) )
        smoothers.append( ('gauss_seidel',{'sweep' : 'symmetric'}) )

        for pre in smoothers:
            for post in smoothers:
                self.run_cases( {'presmoother' : pre, 'postsmoother' : post} )
    
    def test_coarse_solvers(self): 
        solvers = []
        solvers.append('splu')
        solvers.append('lu')
        solvers.append('cg')

        for solver in solvers:
            self.run_cases( {'coarse_solver' : solver} )

class TestComplexParameters(TestCase):
    def setUp(self):
        self.cases = []
        
        # Consider "helmholtz" like problems with an imaginary shift so that the operator 
        #   should still be SPD in a sense and SA should perform well.
        # There are better near nullspace vectors than the default, 
        #   but a constant should give a convergent solver, nonetheless.
        A = poisson( (100,),  format='csr'); A = A + 1.0j*scipy.sparse.eye(A.shape[0], A.shape[1])
        self.cases.append((A, None))
        A = poisson( (10,10),  format='csr'); A = A + 1.0j*scipy.sparse.eye(A.shape[0], A.shape[1])
        self.cases.append((A, None))

    def run_cases(self, opts):
        for A,B in self.cases:
            ml = smoothed_aggregation_solver(A, B, max_coarse=5, **opts)

            numpy.random.seed(0) #make tests repeatable

            x = rand(A.shape[0]) + 1.0j*rand(A.shape[0])
            b = A*rand(A.shape[0])

            x_sol,residuals = ml.solve(b, x0=x, maxiter=30, tol=1e-10, return_residuals=True)
            convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            print convergence_ratio
            assert(convergence_ratio < 0.9)


    def test_strength_of_connection(self): 
        for strength in ['classical', 'symmetric']:                 #,'ode']:
            self.run_cases( {'strength' : strength} )
    
    def test_aggregation_method(self): 
        for aggregate in ['standard','lloyd']:
            self.run_cases( {'aggregate' : aggregate} )
    
    def test_prolongation_smoother(self): 
        for smooth in ['jacobi','richardson','kaczmarz_jacobi', 'kaczmarz_richardson']:    #, 'energy']:
            self.run_cases( {'smooth' : smooth} )

    def test_smoothers(self): 
        smoothers = []
        smoothers.append('gauss_seidel')
        smoothers.append( ('gauss_seidel',{'sweep' : 'symmetric'}) )
        smoothers.append( ('kaczmarz_gauss_seidel',{'sweep' : 'symmetric'}) )

        for pre in smoothers:
            for post in smoothers:
                self.run_cases( {'presmoother' : pre, 'postsmoother' : post} )
    
    def test_coarse_solvers(self): 
        solvers = []
        solvers.append('splu')
        solvers.append('lu')
        solvers.append('cg')
        solvers.append('pinv2')

        for solver in solvers:
            self.run_cases( {'coarse_solver' : solver} )


class TestSolverPerformance(TestCase):
    def setUp(self):
        self.cases = []

        self.cases.append(( poisson( (10000,),  format='csr'), None))
        self.cases.append(( poisson( (100,100), format='csr'), None))
        self.cases.append( linear_elasticity( (100,100), format='bsr') )
        # TODO add unstructured tests


    def test_basic(self):
        """check that method converges at a reasonable rate"""

        for A,B in self.cases:
            ml = smoothed_aggregation_solver(A, B, max_coarse=10)

            numpy.random.seed(0) #make tests repeatable

            x = rand(A.shape[0])
            b = A*rand(A.shape[0])

            x_sol,residuals = ml.solve(b,x0=x,maxiter=20,tol=1e-10,return_residuals=True)

            avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            
            assert(avg_convergence_ratio < 0.3)

    def test_DAD(self):
        A = poisson( (50,50), format='csr' )        

        x = rand(A.shape[0])
        b = rand(A.shape[0])
 
        D     = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
        D_inv = diag_sparse(1.0/D.data)
 
        DAD   = D*A*D
 
        B = ones((A.shape[0],1))
 
        #TODO force 2 level method and check that result is the same
        kwargs = {'max_coarse' : 1, 'max_levels' : 2, 'coarse_solver' : 'splu'}

        sa = smoothed_aggregation_solver(D*A*D, D_inv * B, **kwargs)
        x_sol,residuals = sa.solve(b,x0=x,maxiter=10,tol=1e-12,return_residuals=True)
        avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))

        assert(avg_convergence_ratio < 0.25)


