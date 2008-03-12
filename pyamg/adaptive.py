__all__ = ['adaptive_sa_solver']


import numpy,scipy,scipy.sparse

from numpy import sqrt, ravel, diff, zeros, zeros_like, inner, concatenate, \
                  asarray, hstack, ascontiguousarray, isinf, dot
from numpy.random import randn, rand
from scipy.linalg import norm
from scipy.sparse import csr_matrix, coo_matrix, bsr_matrix

from relaxation import gauss_seidel
from multilevel import multilevel_solver
from utils import approximate_spectral_radius,hstack_csr,vstack_csr,diag_sparse
from sa import sa_standard_aggregation, sa_smoothed_prolongator, \
        sa_fit_candidates, sa_strong_connections





def sa_hierarchy(A,B,AggOps):
    """
    Construct multilevel hierarchy using Smoothed Aggregation
        Inputs:
          A  - matrix
          B  - fine-level near nullspace candidates be approximated

        Ouputs:
          (As,Ps,Ts,Bs) - tuple of lists
                  - As -  
                  - Ps - smoothed prolongators
                  - Ts - tentative prolongators
                  - Bs - near nullspace candidates
    """
    As = [A]
    Ps = []
    Ts = []
    Bs = [B]

    for AggOp in AggOps:
        P,B = sa_fit_candidates(AggOp,B)
        I   = sa_smoothed_prolongator(A,P)
        A   = I.T.asformat(I.format) * A * I
        As.append(A)
        Ts.append(P)
        Ps.append(I)
        Bs.append(B)
    return As,Ps,Ts,Bs


def adaptive_sa_solver(A, num_candidates=1, mu=5, improvement=0, max_levels=10, max_coarse=100, \
        epsilon=0.0, omega=4.0/3.0, symmetric=True, rescale=True,
        aggregation=None):
    """Create a multilevel solver using Adaptive Smoothed Aggregation (aSA)


    Parameters
    ----------

    A : {csr_matrix, bsr_matrix}
        Square matrix in CSR or BSR format
    num_candidates : {integer} : default 1
        Number of near-nullspace candidates to generate.
    mu : {integer} : default 5
        Number of cycles used at each level of the adaptive setup phase.

    Optional Parameters
    -------------------

    improvement : {integer}
        Number of times each candidate is improved
    max_levels: {integer}
        Maximum number of levels to be used in the multilevel solver.
    max_coarse: {integer}
        Maximum number of variables permitted on the coarse grid.
    epsilon: {float}
        Strength of connection parameter used in aggregation.
    omega: {float} : default 4.0/3.0
        Damping parameter used in prolongator smoothing (0 < omega < 2)
    symmetric: {boolean} : default True
        True if A is symmetric, False otherwise
    rescale: {boolean} : default True
        If True, symmetrically rescale A by the diagonal
        i.e. A -> D * A * D,  where D is diag(A)^-0.5
    aggregation: {sequence of csr_matrix objects}
        List of csr_matrix objects that describe a user-defined
        multilevel aggregation of the degrees of freedom.
        For instance [ Agg0, Agg1 ] defines a three-level hierarchy
        where the dimensions of A, Agg0 and Agg1 are compatible, i.e.
        Agg0.shape[1] == A.shape[0] and Agg1.shape[1] == Agg0.shape[0].
  
                    

    Notes
    -----
        Unlike the standard Smoothed Aggregation (SA) method, adaptive SA
        does not require knowledge of near-nullspace candidate vectors.
        Instead, an adaptive procedure computes one or more candidates 
        'from scratch'.  This approach is useful when no candidates are known
        or the candidates have been invalidated due to changes to matrix A.
        

    Example
    -------
        TODO

    References
    ----------

        Brezina, Falgout, MacLachlan, Manteuffel, McCormick, and Ruge
        "Adaptive Smoothed Aggregation ($\alpha$SA) Multigrid"
        SIAM Review Volume 47 ,  Issue 2  (2005)
        http://www.cs.umn.edu/~maclach/research/aSA2.pdf

    """
    
    if A.shape[0] <= max_coarse:
        return multilevel_solver( [A], [] )

    #first candidate
    x,AggOps = asa_initial_setup_stage(A, max_levels = max_levels, \
            max_coarse = max_coarse, mu = mu, epsilon = epsilon, \
            aggregation = aggregation )

    #TODO make sa_fit_candidates work for small Bs
    x /= norm(x)
    
    #create SA using x here
    As,Ps,Ts,Bs = sa_hierarchy(A,x,AggOps)

    for i in range(num_candidates - 1):
        x = asa_general_setup_stage(As,Ps,Ts,Bs,AggOps,mu=mu)
        x /= norm(x)

        B = hstack((Bs[0],x))
        As,Ps,Ts,Bs = sa_hierarchy(A,B,AggOps)

    #improve candidates?
    for i in range(improvement):
        B = Bs[0]
        for i in range(B.shape[1]):
            B = B[:,1:]
            As,Ps,Ts,Bs = sa_hierarchy(A,B,AggOps)
            x = asa_general_setup_stage(As,Ps,Ts,Bs,AggOps,mu)
            B = hstack((B,x))
        As,Ps,Ts,Bs = sa_hierarchy(A,B,AggOps)

    return multilevel_solver(As,Ps)

def asa_initial_setup_stage(A, max_levels, max_coarse, mu, epsilon, aggregation):
    """Computes a complete aggregation and the first near-nullspace candidate


    """
    if aggregation is not None:
        max_coarse = 0
        max_levels = len(aggregation) + 1

    # aSA parameters
    # mu      - number of test relaxation iterations
    # epsilon - minimum acceptable relaxation convergence factor

    #step 1
    A_l = A
    x   = rand(A_l.shape[0],1) # TODO see why randn() fails here
    skip_f_to_i = False

    #step 2
    gauss_seidel(A_l, x, zeros_like(x), iterations=mu, sweep='symmetric')

    #step 3
    #TODO test convergence rate here

    As     = [A]
    Ps     = []
    AggOps = []

    while len(AggOps) + 1 < max_levels and A_l.shape[0] > max_coarse:
        if aggregation is None:
            C_l = sa_strong_connections(A_l,epsilon)
            W_l = sa_standard_aggregation(C_l) #step 4b
        else:
            W_l = aggregation[len(AggOps)]
        P_l,x = sa_fit_candidates(W_l,x)                   #step 4c
        I_l   = sa_smoothed_prolongator(A_l,P_l)           #step 4d
        A_l   = I_l.T.asformat(I_l.format) * A_l * I_l     #step 4e
        #TODO change variable names I_l -> P, P_l -> T

        AggOps.append(W_l)
        Ps.append(I_l)
        As.append(A_l)

        if A_l.shape <= max_coarse:  break

        if not skip_f_to_i:
            x_hat = x.copy()                                                   #step 4g
            gauss_seidel(A_l,x,zeros_like(x),iterations=mu,sweep='symmetric')  #step 4h
            x_A_x = dot(x.T,A_l*x)
            if (x_A_x/dot(x_hat.T,A_l*x_hat))**(1.0/mu) < epsilon:             #step 4i
                print "sufficient convergence, skipping"
                skip_f_to_i = True
                if x_A_x == 0:
                    x = x_hat  #need to restore x

    #update fine-level candidate
    for A_l,I in reversed(zip(As[1:],Ps)):
        gauss_seidel(A_l,x,zeros_like(x),iterations=mu,sweep='symmetric')         #TEST
        x = I * x
    gauss_seidel(A,x,zeros_like(x),iterations=mu,sweep='symmetric')         #TEST

    return x,AggOps  #first candidate,aggregation


def asa_general_setup_stage(As, Ps, Ts, Bs, AggOps, mu):
    A = As[0]

    x = rand(A.shape[0],1)
    b = zeros_like(x)

    x = multilevel_solver(As,Ps).solve(b, x0=x, tol=1e-10, maxiter=mu)

    #TEST FOR CONVERGENCE HERE

    temp_Ps = []
    temp_As = [A]

    def make_bridge(P):
        M,N  = P.shape
        K    = P.blocksize[0]
        bnnz = P.indptr[-1]
        data = zeros( (bnnz, K+1, K), dtype=P.dtype )
        data[:,:-1,:] = P.data
        return bsr_matrix( (data, P.indices, P.indptr), shape=( (K+1)*(M/K), N) )

    for i in range(len(As) - 2):
        B_old = Bs[i]
        B = zeros( (x.shape[0], B_old.shape[1] + 1), dtype=x.dtype)

        B[:B_old.shape[0],:B_old.shape[1]] = B_old
        B[:,-1] = x.reshape(-1)

        T,R = sa_fit_candidates(AggOps[i],B)

        P = sa_smoothed_prolongator(A,T)
        A = P.T.asformat(P.format) * A * P

        temp_Ps.append(P)
        temp_As.append(A)

        bridge = make_bridge(Ps[i+1])

        solver = multilevel_solver( [A] + As[i+2:], [bridge] + Ps[i+2:] )

        x = R[:,-1].reshape(-1,1)
        x = solver.solve(zeros_like(x), x0=x, tol=1e-8, maxiter=mu)

    for A,P in reversed(zip(temp_As,temp_Ps)):
        x = P * x
        gauss_seidel(A,x,zeros_like(x),iterations=mu,sweep='symmetric')

    return x

#    def __augment_cycle(self,As,Ps,Ts,Bs,AggOps,x):
#        A = As[0]
#
#        new_As = [A]
#        new_Ps = []
#        new_Ts = []
#        new_Bs = [ hstack((Bs[0],x)) ]
#
#        for i in range(len(As) - 1):
#            T,R = augment_candidates(AggOps[i], Ts[i], Bs[i+1], x)
#
#            P = sa_smoothed_prolongator(A,T)
#            A = P.T.asformat(P.format) * A * P
#
#            new_As.append(A)
#            new_Ps.append(P)
#            new_Ts.append(T)
#            new_Bs.append(R)
#
#            x = R[:,-1].reshape(-1,1)
#
#        return new_As,new_Ps,new_Ts,new_Bs


#def augment_candidates(AggOp, old_Q, old_R, new_candidate):
#    #TODO update P and A also
#
#    K = old_R.shape[1]
#
#    #determine blocksizes
#    if new_candidate.shape[0] == old_Q.shape[0]:
#        #then this is the first prolongator
#        old_bs = (1,K)
#        new_bs = (1,K+1)
#    else:
#        old_bs = (K,K)
#        new_bs = (K+1,K+1)
#
#    #AggOp = expand_into_blocks(AggOp,new_bs[0],1).tocsr() #TODO switch to block matrix
#
#    # tentative prolongator
#    Q_data = zeros( (AggOp.indptr[-1],) + new_bs)
#    Q_data[:,:old_bs[0],:old_bs[1]] = old_Q.data.reshape((-1,) + old_bs)
#
#    # coarse candidates
#    R = zeros((AggOp.shape[1],K+1,K+1))
#    R[:,:K,:K] = old_R.reshape(-1,K,K)
#
#    c = new_candidate.reshape(-1)[diff(AggOp.indptr) == 1]  #eliminate DOFs that aggregation misses
#    threshold = 1e-10 * abs(c).max()   # cutoff for small basis functions
#
#    X = bsr_matrix((c,AggOp.indices,AggOp.indptr),shape=AggOp.shape)
#
#    #orthogonalize X against previous
#    for i in range(K):
#        old_c = ascontiguousarray(Q_data[:,:,i].reshape(-1))
#        D_AtX = csr_matrix((old_c*X.data,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X
#        R[:,i,K] = D_AtX
#        X.data -= D_AtX[X.indices] * old_c
#
#    #normalize X
#    D_XtX = csr_matrix((X.data**2,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X
#    col_norms = sqrt(D_XtX)
#    mask = col_norms < threshold  # find small basis functions
#    col_norms[mask] = 0           # and set them to zero
#
#    R[:,K,K] = col_norms      # store diagonal entry into R
#
#    col_norms = 1.0/col_norms
#    col_norms[mask] = 0
#    X.data *= col_norms[X.indices]
#    Q_data[:,:,-1] = X.data.reshape(-1,new_bs[0])
#
#    Q_data = Q_data.reshape(-1)  #TODO BSR change
#    R = R.reshape(-1,K+1)
#
#    Q = csr_matrix((Q_data,Q_indices,Q_indptr),shape=(AggOp.shape[0],(K+1)*AggOp.shape[1]))
#
#    return Q,R


