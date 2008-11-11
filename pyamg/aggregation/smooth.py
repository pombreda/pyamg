"""Methods to smooth tentative prolongation operators"""

__docformat__ = "restructuredtext en"

from numpy import ones
from pyamg.utils import approximate_spectral_radius, scale_rows, get_diagonal

__all__ = ['jacobi_prolongation_smoother', 'richardson_prolongation_smoother', 
        'energy_prolongation_smoother', 'kaczmarz_richardson_prolongation_smoother',
        'kaczmarz_jacobi_prolongation_smoother']


def jacobi_prolongation_smoother(S, T, omega=4.0/3.0, degree=1):
    """Jacobi prolongation smoother
   
    Parameters
    ----------
    S : {csr_matrix, bsr_matrix}
        Sparse NxN matrix used for smoothing.  Typically, A or the
        "filtered matrix" obtained from A by lumping weak connections
        onto the diagonal of A.
    T : {csr_matrix, bsr_matrix}
        Tentative prolongator
    omega : {scalar}
        Damping parameter

    Returns
    -------
    P : {csr_matrix, bsr_matrix}
        Smoothed (final) prolongator defined by P = (I - omega/rho(K) K) * T
        where K = diag(S)^-1 * S and rho(K) is an approximation to the 
        spectral radius of K.

    """

    D = S.diagonal()
    D_inv = 1.0 / D
    D_inv[D == 0] = 0

    D_inv_S = scale_rows(S, D_inv, copy=True)
    D_inv_S *= omega/approximate_spectral_radius(D_inv_S)

    P = T
    for i in range(degree):
        P = P - (D_inv_S*P)

    return P


def richardson_prolongation_smoother(S, T, omega=4.0/3.0, degree=1):
    """Richardson prolongation smoother
   
    Parameters
    ----------
    S : {csr_matrix, bsr_matrix}
        Sparse NxN matrix used for smoothing.  Typically, A or the
        "filtered matrix" obtained from A by lumping weak connections
        onto the diagonal of A.
    T : {csr_matrix, bsr_matrix}
        Tentative prolongator
    omega : {scalar}
        Damping parameter

    Returns
    -------
    P : {csr_matrix, bsr_matrix}
        Smoothed (final) prolongator defined by P = (I - omega/rho(S) S) * T
        where rho(S) is an approximation to the spectral radius of S.

    """

    weight = omega/approximate_spectral_radius(S)

    P = T
    for i in range(degree):
        P = P - weight*(S*P)

    return P

def kaczmarz_jacobi_prolongation_smoother(S, T, omega=4.0/3.0, degree=1):
    """Jacobi prolongation smoother for the normal equations (i.e. Kaczmarz)
   
    Parameters
    ----------
    S : {csr_matrix, bsr_matrix}
        Sparse NxN matrix used for smoothing.  Typically, A or the
        "filtered matrix" obtained from A by lumping weak connections
        onto the diagonal of A.
    T : {csr_matrix, bsr_matrix}
        Tentative prolongator
    omega : {scalar}
        Damping parameter

    Returns
    -------
    P : {csr_matrix, bsr_matrix}
        Smoothed (final) prolongator
    """

    # Form Dinv for S*S.H
    D_inv = get_diagonal(S, norm_eq=2, inv=True)
    D_inv_S = scale_rows(S, D_inv, copy=True)

    # Approximate Spectral radius by defining a matvec for S.T*D_inv_S
    ST = S.conjugate().T.asformat(D_inv_S.format)
    class matvec_mat:
    
        def __init__(self, matvec, shape, dtype):
            self.shape = shape
            self.matvec = matvec
            self.__mul__ = matvec
            self.dtype = dtype
    
    def matmul(A,B,x):
        return A*(B*x)
    
    StDS_matmul = lambda x:matmul(ST, D_inv_S, x)
    StDS = matvec_mat(StDS_matmul, S.shape, S.dtype)
    omega = omega/approximate_spectral_radius(StDS)

    P = T
    for i in range(degree):
        P = P - omega*(ST*(D_inv_S*P))

    return P

def kaczmarz_richardson_prolongation_smoother(S, T, omega=4.0/3.0, degree=1):
    """Richardson prolongation smoother for the normal equations (i.e. Kaczmarz)
   
    Parameters
    ----------
    S : {csr_matrix, bsr_matrix}
        Sparse NxN matrix used for smoothing.  Typically, A or the
        "filtered matrix" obtained from A by lumping weak connections
        onto the diagonal of A.
    T : {csr_matrix, bsr_matrix}
        Tentative prolongator
    omega : {scalar}
        Damping parameter

    Returns
    -------
    P : {csr_matrix, bsr_matrix}
        Smoothed (final) prolongator 
    """

    # Approximate Spectral radius by defining a matvec for S*S.H
    ST = S.conjugate().T.asformat(S.format)
    class matvec_mat:
    
        def __init__(self, matvec, shape, dtype):
            self.shape = shape
            self.matvec = matvec
            self.__mul__ = matvec
            self.dtype = dtype
    
    def matmul(A,B,x):
        return A*(B*x)
    
    StS_matmul = lambda x:matmul(ST, S, x)
    StS = matvec_mat(StS_matmul, S.shape, S.dtype)
    omega = omega/approximate_spectral_radius(StS)

    P = T
    for i in range(degree):
        P = P - omega*(ST*(S*P))

    return P


""" sa_energy_min + helper functions minimize the energy of a tentative prolongator for use in SA """

from numpy import ones, zeros, asarray, dot, array_split, diff, ravel, asarray, ones_like, conjugate, mat
from scipy.sparse import csr_matrix, isspmatrix_csr, bsr_matrix, isspmatrix_bsr
from scipy.linalg import pinv2
from pyamg.utils import UnAmal
import pyamg.multigridtools

########################################################################################################
#   Helper function for the energy minimization prolongator generation routine

def Satisfy_Constraints(U, B, BtBinv):
    """Update U to satisfy U*B = 0

    Input
    =====
    U                     BSR Matrix to operate on
    B                     Near nullspace vectors
    BtBinv                Local inv(B'*B) matrices for each dof, i.  
        
    Output
    ======
    Updated U, so that U*B = 0.  Update is computed by orthogonally (in 2-norm)
    projecting out the components of span(B) in U in a row-wise fashion

    """
    
    RowsPerBlock = U.blocksize[0]
    ColsPerBlock = U.blocksize[1]
    num_blocks = U.indices.shape[0]
    num_block_rows = U.shape[0]/RowsPerBlock

    UB = ravel(U*B)

    # Apply constraints, noting that we need the conjugate of B 
    # for use as Bi.H in local projection
    pyamg.multigridtools.satisfy_constraints_helper(RowsPerBlock, ColsPerBlock, 
            num_blocks, num_block_rows, 
            conjugate(ravel(B)), UB, ravel(BtBinv), 
            U.indptr, U.indices, ravel(U.data))
        
    return U

    

########################################################################################################


def energy_prolongation_smoother(A, T, Atilde, B, SPD=True, maxiter=4, tol=1e-8, degree=1):
    """Minimize the energy of the coarse basis functions (columns of T)

    Parameters
    ----------

    A : {csr_matrix, bsr_matrix}
        Sparse NxN matrix
    T : {bsr_matrix}
        Tentative prolongator, a NxM sparse matrix (M < N)
    Atilde : {csr_matrix}
        Strength of connection matrix
    B : {array}
        Near-nullspace modes for coarse grid.  Has shape (M,k) where
        k is the number of coarse candidate vectors.
    SPD : boolean
        Booolean denoting symmetric (or Hermitian) positive-definiteness of A
    maxiter : integer
        Number of energy minimization steps to apply to the prolongator
    tol : scalar
        Minimization tolerance
   
    Returns
    -------
    P : {bsr_matrix}
        Smoothed prolongator

    References
    ----------

        Jan Mandel, Marian Brezina, and Petr Vanek
        "Energy Optimization of Algebraic Multigrid Bases"
        Computing 62, 205-228, 1999
        http://dx.doi.org/10.1007/s006070050022
    
    """
    
    #====================================================================
    #Test Inputs
    if maxiter < 0:
        raise ValueError('maxiter must be > 0')
    if tol > 1:
        raise ValueError('tol must be <= 1') 
   
    if isspmatrix_csr(A):
        A = A.tobsr(blocksize=(1,1), copy=False)
    elif isspmatrix_bsr(A):
        pass
    else:
        raise TypeError("A must be csr_matrix or bsr_matrix")

    if Atilde is None:
        Atilde = csr_matrix( (ones(len(A.indices)), A.indices.copy(), A.indptr.copy()), shape=(A.shape[0]/A.blocksize[0], A.shape[1]/A.blocksize[1]))

    if not isspmatrix_csr(Atilde):
        raise TypeError("Atilde must be csr_matrix")

    if T.blocksize[0] != A.blocksize[0]:
        raise ValueError("T's row-blocksize should be the same as A's blocksize")

    if min(T.nnz, Atilde.nnz, A.nnz) == 0:
        return T

    #====================================================================
    
    
    #====================================================================
    # Retrieve problem information
    Nfine = T.shape[0]
    Ncoarse = T.shape[1]
    NullDim = B.shape[1]
    #Number of PDEs per point is defined implicitly by block size
    numPDEs = A.blocksize[0]
    #====================================================================
    
    
    #====================================================================
    # Unamalgate Atilde if (numPDEs > 1)

    # UnAmal returns a BSR matrix, so the mat-mat will be between BSR mats. 
    #TODO replace large matmat with smaller matmat, then expand
    T.sort_indices()
    Sparsity_Pattern = bsr_matrix( (ones_like(T.data), T.indices, T.indptr), shape=T.shape)
    X = UnAmal(Atilde, numPDEs, numPDEs)
    for i in range(degree):
        Sparsity_Pattern = X * Sparsity_Pattern
    del X
    Sparsity_Pattern.data[:] = 1.0
    Sparsity_Pattern.sort_indices()

    
    #====================================================================
    #Construct array of inv(Bi'Bi), where Bi is B restricted to row i's sparsity pattern in 
    #   Sparsity Pattern.  This array is used multiple times in the Satisfy_Constraints routine.

    ColsPerBlock = Sparsity_Pattern.blocksize[1]
    RowsPerBlock = Sparsity_Pattern.blocksize[0]
    Nnodes = Nfine/RowsPerBlock

    BtBinv = zeros((Nnodes,NullDim,NullDim), dtype=B.dtype) 
    BsqCols = sum(range(NullDim+1))
    Bsq = zeros((Ncoarse,BsqCols), dtype=B.dtype)
    counter = 0
    for i in range(NullDim):
        for j in range(i,NullDim):
            Bsq[:,counter] = conjugate(ravel(asarray(B[:,i])))*ravel(asarray(B[:,j]))
            counter = counter + 1
    
    pyamg.multigridtools.invert_BtB(NullDim, Nnodes, ColsPerBlock, ravel(asarray(Bsq)), 
        BsqCols, ravel(asarray(BtBinv)), Sparsity_Pattern.indptr, Sparsity_Pattern.indices)
    # TODO extend pseudoinverse in C++ to complex
    for i in range(Nnodes):
        BtBinv[i,:,:] = pinv2(BtBinv[i,:,:]) 
    #====================================================================
    
    #====================================================================
    #Iteratively minimize the energy of T subject to the constraints of Sparsity_Pattern
    #   and maintaining T's effect on B, i.e. T*B = (T+Update)*B, i.e. Update*B = 0 
    i = 0
    if SPD:
        #Apply CG with diagonal preconditioning
        Dinv = get_diagonal(A, norm_eq=False, inv=True)

        #Calculate initial residual
        R = -A*T  
        
        #Enforce constraints on R.  First the sparsity pattern, then the nullspace vectors.
        R = R.multiply(Sparsity_Pattern)
        if R.nnz < Sparsity_Pattern.nnz:
            # ugly hack to give R the same sparsity pattern as Sparsity_Pattern
            # It is dangerous to leave the 1e-100 values in there as this can give
            # coarser levels 1e-100 type entries on the diagonal and mess everything up
            R = R + 1e-100*Sparsity_Pattern 
            Rshape = R.data.shape
            R.data = R.data.reshape(-1,)
            R.data[R.data == 1e-100] = 0.0
            R.data = R.data.reshape(Rshape)

        Satisfy_Constraints(R, B, BtBinv)
    
        if R.nnz == 0:
            print "Error in sa_energy_min(..).  Initial R no nonzeros on a level.  Calling Default Prolongator Smoother\n"
            return jacobi_prolongation_smoother(Atilde, T)
        
        #Calculate max norm of the residual
        resid = abs(R.data.flatten()).max()
        #print "Energy Minimization of Prolongator --- Iteration 0 --- r = " + str(resid)

        while i < maxiter and resid > tol:
            #Apply diagonal preconditioner
            Z = scale_rows(R, Dinv)
    
            #Frobenius innerproduct of (R,Z) = sum( conjugate(rk).*zk)
            newsum = (R.conjugate().multiply(Z)).sum()
                
            #P is the search direction, not the prolongator, which is T.    
            if(i == 0):
                P = Z
            else:
                beta = newsum/oldsum
                P = Z + beta*P
            oldsum = newsum
    
            #Calculate new direction and enforce constraints
            AP = A*P
            AP = AP.multiply(Sparsity_Pattern)
            if AP.nnz < Sparsity_Pattern.nnz:
                # ugly hack to give AP the same sparsity pattern as Sparsity_Pattern
                # It is dangerous to leave the 1e-100 values in there as this can give
                # coarser levels 1e-100 type entries on the diagonal and mess everything up
                AP = AP + 1e-100*Sparsity_Pattern
                APshape = AP.data.shape
                AP.data = AP.data.reshape(-1,)
                AP.data[AP.data == 1e-100] = 0.0
                AP.data = AP.data.reshape(APshape)           
            
            Satisfy_Constraints(AP, B, BtBinv)
            
            #Frobenius innerproduct of (P, AP)
            alpha = newsum/(P.conjugate().multiply(AP)).sum()
    
            #Update the prolongator, T
            T = T + alpha*P 
    
            #Update residual
            R = R - alpha*AP
            
            i += 1
            resid = abs(R.data).max()
            #print "Energy Minimization of Prolongator --- Iteration " + str(i) + " --- r = " + str(resid)
     
    else:   
        #For non-SPD system, apply CG on Normal Equations with Diagonal Preconditioning (requires transpose)
        Ah = A.H
        
        # D for A.H*A
        Dinv = get_diagonal(A, norm_eq=1, inv=True)

        #Calculate initial residual
        R = -Ah*(A*T)  
        
        #Enforce constraints on R.  First the sparsity pattern, then the nullspace vectors.
        R = R.multiply(Sparsity_Pattern)
        if R.nnz < Sparsity_Pattern.nnz:
            # ugly hack to give R the same sparsity pattern as Sparsity_Pattern
            # It is dangerous to leave the 1e-100 values in there as this can give
            # coarser levels 1e-100 type entries on the diagonal and mess everything up
            R = R + 1e-100*Sparsity_Pattern 
            Rshape = R.data.shape
            R.data = R.data.reshape(-1,)
            R.data[R.data == 1e-100] = 0.0
            R.data = R.data.reshape(Rshape)

        Satisfy_Constraints(R, B, BtBinv)
    
        if R.nnz == 0:
            print "Error in sa_energy_min(..).  Initial R no nonzeros on a level.  Calling Default Prolongator Smoother\n"
            return jacobi_prolongation_smoother(Atilde, T)
        
        #Calculate max norm of the residual
        resid = abs(R.data.flatten()).max()
        #print "Energy Minimization of Prolongator --- Iteration 0 --- r = " + str(resid)

        while i < maxiter and resid > tol:
            #Apply diagonal preconditioner
            Z = scale_rows(R, Dinv)
    
            #Frobenius innerproduct of (R,Z) = sum(rk.*zk)
            newsum = (R.conjugate().multiply(Z)).sum()
                
            #P is the search direction, not the prolongator, which is T.    
            if(i == 0):
                P = Z
            else:
                beta = newsum/oldsum
                P = Z + beta*P
            oldsum = newsum
    
            #Calculate new direction and enforce constraints
            AP = Ah*(A*P)
            AP = AP.multiply(Sparsity_Pattern)
            if AP.nnz < Sparsity_Pattern.nnz:
                # ugly hack to give AP the same sparsity pattern as Sparsity_Pattern
                # It is dangerous to leave the 1e-100 values in there as this can give
                # coarser levels 1e-100 type entries on the diagonal and mess everything up
                AP = AP + 1e-100*Sparsity_Pattern
                APshape = AP.data.shape
                AP.data = AP.data.reshape(-1,)
                AP.data[AP.data == 1e-100] = 0.0
                AP.data = AP.data.reshape(APshape)           
            
            Satisfy_Constraints(AP, B, BtBinv)
            
            #Frobenius innerproduct of (P, AP)
            alpha = newsum/(P.conjugate().multiply(AP)).sum()
    
            #Update the prolongator, T
            T = T + alpha*P 
    
            #Update residual
            R = R - alpha*AP
            
            i += 1
            resid = abs(R.data).max()
            #print "Energy Minimization of Prolongator --- Iteration " + str(i) + " --- r = " + str(resid)
    
#        # Standard CGNE type algorithm
#        # It is not straight forward how to apply the constraints to this algorithm, i.e. 
#        # the search direction, A A.T P, is never explicitly formed, so its not clear where
#        # to apply the constraints
#
#        #Calculate initial residual
#        R = -A*T  
#        
#        #Calculate max norm of the residual
#        resid = abs(R.data.flatten()).max()
#        #print "Energy Minimization of Prolongator --- Iteration 0 --- r = " + str(resid)
#
#        # Apply diagonal preconditioner
#        Z = scale_rows(R, Dinv)
#        
#        # P is the search direction, not prolongator, which is T.
#        P = Ah*Z
#
#        #Enforce constraints on P.  First the sparsity pattern, then the nullspace vectors.
#        P = P.multiply(Sparsity_Pattern)
#        if P.nnz < Sparsity_Pattern.nnz:
#            # ugly hack to give P the same sparsity pattern as Sparsity_Pattern
#            # It is dangerous to leave the 1e-100 values in there as this can give
#            # coarser levels 1e-100 type entries on the diagonal and mess everything up
#            P = P + 1e-100*Sparsity_Pattern 
#            Pshape = P.data.shape
#            P.data = P.data.reshape(-1,)
#            P.data[P.data == 1e-100] = 0.0
#            P.data = P.data.reshape(Pshape)
#
#        Satisfy_Constraints(P, B, BtBinv)
#    
#        if P.nnz == 0:
#            print "Error in sa_energy_min(..).  Initial R no nonzeros on a level.  Calling Default Prolongator Smoother\n"
#            return jacobi_prolongation_smoother(Atilde, T)
#
#
#        # Frobenius innerproduct of (R,Z) = sum(rk.*zk)
#        old_zr = (R.multiply(Z)).sum()
#        
#        while i < maxiter and resid > tol:
#            
#            #Frobenius innerproduct of (P, P)
#            alpha = old_zr/(P.multiply(P)).sum()
#            
#            #Update the prolongator, T
#            T = T + alpha*P 
#            
#            #Calculate new direction and enforce constraints
#            AP = A*P
#            #AP = AP.multiply(Sparsity_Pattern)
#            #if AP.nnz < Sparsity_Pattern.nnz:
#            #    # ugly hack to give AP the same sparsity pattern as Sparsity_Pattern
#            #    AP = AP + 1e-100*Sparsity_Pattern
#            #    APshape = AP.data.shape
#            #    AP.data = AP.data.reshape(-1,)
#            #    AP.data[AP.data == 1e-100] = 0.0
#            #    AP.data = AP.data.reshape(APshape)           
#            #
#            #Satisfy_Constraints(AP, B, BtBinv)
#        
#            # Update residual
#            R = R - alpha*AP
#            
#            # Apply diagonal preconditioner
#            Z = scale_rows(R, Dinv)
#        
#            # Frobenius innerproduct of (R,Z) = sum(rk.*zk)
#            new_zr = (R.multiply(Z)).sum()
#            beta = new_zr/old_zr
#            old_zr = new_zr
#
#            # Update search direction
#            P = Ah*Z + beta*P
#
#            #Calculate new direction and enforce constraints
#            P = P.multiply(Sparsity_Pattern)
#            if P.nnz < Sparsity_Pattern.nnz:
#                # ugly hack to give AP the same sparsity pattern as Sparsity_Pattern
#                P = P + 1e-100*Sparsity_Pattern
#                Pshape = P.data.shape
#                P.data = P.data.reshape(-1,)
#                P.data[P.data == 1e-100] = 0.0
#                P.data = P.data.reshape(Pshape)           
#
#            Satisfy_Constraints(P, B, BtBinv)
#
#                
#            i += 1
#            resid = abs(R.data).max()
#            #print "Energy Minimization of Prolongator --- Iteration " + str(i) + " --- r = " + str(resid)

# experimentally and theoretically this does not seem like a good idea
#        ## Convert T back to original system
#        #T = Ah*T
#        #T = T.multiply(Sparsity_Pattern)
#        #if T.nnz < Sparsity_Pattern.nnz:
#        #    # ugly hack to give AP the same sparsity pattern as Sparsity_Pattern
#        #    T = T + 1e-100*Sparsity_Pattern
#        #Satisfy_Constraints(T, B, BtBinv)

# Previous non-SPD minimization strategy
#        #Apply min-res to the nonsymmetric system
#        while i < maxiter and resid > tol:
#    
#            #P is the search direction, not the prolongator
#            P = A*R
#    
#            #Enforce constraints on P
#            P = P.multiply(Sparsity_Pattern)
#            if P.nnz < Sparsity_Pattern.nnz:
#                # ugly hack to give P the same sparsity pattern as Sparsity_Pattern
#                P = P + 1e-100*Sparsity_Pattern
#                Pshape = P.data.shape
#                P.data = P.data.reshape(-1,)
#                P.data[P.data == 1e-100] = 0.0
#                P.data = P.data.reshape(Pshape)           
#
#            Satisfy_Constraints(P, B, BtBinv)
#    
#            #Frobenius innerproduct of (P, R)
#            numer = (P.multiply(R)).sum()
#            
#            #Frobenius innerproduct of (P, P)
#            denom = (P.multiply(P)).sum()
#    
#            alpha = numer/denom
#    
#            #Update prolongator
#            T = T + alpha*R
#    
#            #Update residual
#            R = R - alpha*P
#            
#            i += 1
#            resid = max(R.data.flatten().__abs__())
#            #print "Energy Minimization of Prolongator --- Iteration " + str(i) + " --- r = " + str(resid)
    #====================================================================
    
    T.eliminate_zeros()
    return T

