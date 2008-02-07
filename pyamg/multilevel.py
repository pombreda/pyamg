__all__ = ['multilevel_solver']

import scipy
import numpy
from numpy import ones, zeros, zeros_like, array, asarray, empty
from numpy.linalg import norm
from scipy.splinalg import spsolve

from relaxation import gauss_seidel,jacobi,sor
from utils import symmetric_rescaling, diag_sparse


class multilevel_solver:
    def __init__(self, As, Ps, Rs=None, preprocess=None, postprocess=None):
        self.As = As
        self.Ps = Ps
        self.preprocess  = preprocess
        self.postprocess = postprocess

        if Rs is None:
            self.Rs = [P.T for P in self.Ps]
        else:
            self.Rs = Rs

    def __repr__(self):
        output = 'multilevel_solver\n'
        output += 'Number of Levels:     %d\n' % len(self.As)
        output += 'Operator Complexity: %6.3f\n' % self.operator_complexity()
        output += 'Grid Complexity:     %6.3f\n' % self.grid_complexity()

        total_nnz =  sum([A.nnz for A in self.As])

        output += '  level   unknowns     nonzeros\n'
        for n,A in enumerate(self.As):
            output += '   %2d   %10d   %10d [%5.2f%%]\n' % (n,A.shape[1],A.nnz,(100*float(A.nnz)/float(total_nnz)))

        return output

    def operator_complexity(self):
        """number of nonzeros on all levels / number of nonzeros on the finest level"""
        return sum([A.nnz for A in self.As])/float(self.As[0].nnz)

    def grid_complexity(self):
        """number of unknowns on all levels / number of unknowns on the finest level"""
        return sum([A.shape[0] for A in self.As])/float(self.As[0].shape[0])


    def psolve(self, b):
        return self.solve(b,maxiter=1)

    def solve(self, b, x0=None, tol=1e-5, maxiter=100, callback=None, return_residuals=False):
        """
        TODO
        """

        if x0 is None:
            x = zeros_like(b)
        else:
            x = array(x0) #copy

        if self.preprocess is not None:
            x,b = self.preprocess(x,b)

        #TODO change use of tol (relative tolerance) to agree with other iterative solvers
        A = self.As[0]
        residuals = [ norm(b-A*x) ]

        while len(residuals) <= maxiter and residuals[-1]/residuals[0] > tol:
            self.__solve(0,x,b)

            residuals.append( norm(b-A*x) )

            if callback is not None:
                callback(x)

        if self.postprocess is not None:
            x = self.postprocess(x)

        if return_residuals:
            return x,residuals
        else:
            return x


    def __solve(self,lvl,x,b):
        A = self.As[lvl]

        if len(self.As) == 1:
            #TODO make spsolve preserve dimensions
            x[:] = spsolve(A.tocsc(),b).reshape(x.shape)
            return

        self.presmoother(A,x,b)

        residual = b - A*x

        coarse_b = self.Rs[lvl] * residual
        coarse_x = zeros_like(coarse_b)

        if lvl == len(self.As) - 2:
            #use direct solver on coarsest level
            #TODO reuse factors for efficiency?
            coarse_x[:] = spsolve(self.As[-1].tocsc(),coarse_b).reshape(coarse_x.shape)
            #coarse_x[:] = scipy.linalg.cg(self.As[-1],coarse_b,tol=1e-12)[0].reshape(coarse_x.shape)
            #A_inv = asarray(scipy.linalg.pinv2(self.As[-1].todense()))
            #coarse_x[:] = scipy.dot(A_inv,coarse_b)
            #print "coarse residual norm",scipy.linalg.norm(coarse_b - self.As[-1]*coarse_x)
        else:
            self.__solve(lvl+1,coarse_x,coarse_b)

        x += self.Ps[lvl] * coarse_x   #coarse grid correction

        self.postsmoother(A,x,b)


    def presmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="symmetric")

    def postsmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="symmetric")


