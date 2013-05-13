"""Support for aggregation-based AMG"""

__docformat__ = "restructuredtext en"

import numpy
import scipy
import types
from warnings import warn
from scipy.sparse import csr_matrix, isspmatrix_csr, isspmatrix_bsr, eye

from pyamg import relaxation
from pyamg import amg_core
from pyamg.aggregation.aggregation import preprocess_str_or_agg, preprocess_smooth, preprocess_Bimprove
from pyamg.multilevel import multilevel_solver
from pyamg.relaxation.smoothing import change_smoothers
from pyamg.util.utils import symmetric_rescaling_sa, diag_sparse, amalgamate, \
                             relaxation_as_linear_operator, scale_rows, \
                             get_diagonal, scale_T, get_Cpt_params, \
                             eliminate_diag_dom_nodes
from pyamg.util.linalg import pinv_array, approximate_spectral_radius, \
                              _approximate_eigenvalues
from pyamg.strength import classical_strength_of_connection, \
        symmetric_strength_of_connection, evolution_strength_of_connection, \
        energy_based_strength_of_connection, distance_strength_of_connection
from aggregate import standard_aggregation, naive_aggregation, lloyd_aggregation
from tentative import fit_candidates
from smooth import jacobi_prolongation_smoother, richardson_prolongation_smoother, \
        energy_prolongation_smoother

__all__ = ['rootnode_solver']
 
def blocksize(A):
    # Helper Function: return the blocksize of a matrix 
    if isspmatrix_bsr(A):
        return A.blocksize[0]
    else:
        return 1

def nPDEs(levels):
    # Helper Function:return number of PDEs (i.e. blocksize) at coarsest level
    return blocksize(levels[-1].A) 

def rootnode_solver(A, B=None, BH=None,
        symmetry='hermitian', strength='symmetric', 
        aggregate='standard', smooth='energy',
        presmoother=('block_gauss_seidel',{'sweep':'symmetric'}),
        postsmoother=('block_gauss_seidel',{'sweep':'symmetric'}),
        Bimprove='default', max_levels = 10, max_coarse = 500, 
        coarsen_diag_dom=False, keep=False, **kwargs):
    """
    Create a multilevel solver using root-node based Smoothed Aggregation (SA).  
    See the notes below, for the major differences with the classical-style 
    smoothed aggregation solver in aggregation.smoothed_aggregation_solver.

    Parameters
    ----------
    A : {csr_matrix, bsr_matrix}
        Sparse NxN matrix in CSR or BSR format
    B : {None, array_like}
        Right near-nullspace candidates stored in the columns of an NxK array.
        K must be >= the blocksize of A (see reference [2]). The default value
        B=None is equivalent to choosing the constant over each block-variable,
        B = np.kron(np.ones((A.shape[0]/blocksize(A),1)), np.eye(blocksize(A)))
    BH : {None, array_like}
        Left near-nullspace candidates stored in the columns of an NxK array.
        BH is only used if symmetry is 'nonsymmetric'.  K must be >= the
        blocksize of A (see reference [2]). The default value B=None is
        equivalent to choosing the constant over each block-variable,
        B = np.kron(np.ones((A.shape[0]/blocksize(A),1)), np.eye(blocksize(A)))
    symmetry : {string}
        'symmetric' refers to both real and complex symmetric
        'hermitian' refers to both complex Hermitian and real Hermitian
        'nonsymmetric' i.e. nonsymmetric in a hermitian sense
        Note that for the strictly real case, symmetric and hermitian are the same
        Note that this flag does not denote definiteness of the operator.
    strength : {list} : default ['symmetric', 'classical', 'evolution', ('predefined', {'C' : csr_matrix}), None]
        Method used to determine the strength of connection between unknowns of
        the linear system.  Method-specific parameters may be passed in using a
        tuple, e.g. strength=('symmetric',{'theta' : 0.25 }). If strength=None,
        all nonzero entries of the matrix are considered strong.
        See notes below for varying this parameter on a per level basis.  Also,
        see notes below for using a predefined strength matrix on each level.
    aggregate : {list} : default ['standard', 'lloyd', 'naive', ('predefined', {'AggOp' : csr_matrix})]
        Method used to aggregate nodes.  See notes below for varying this
        parameter on a per level basis.  Also, see notes below for using a
        predefined aggregation on each level.
    smooth : {list} : default ['energy', None]
        Method used to smooth the tentative prolongator.  Method-specific
        parameters may be passed in using a tuple, e.g.  smooth=
        ('energy',{'krylov' : 'gmres'}).  Only 'energy' and None are valid
        prolongation smoothing options.  See notes below for varying this
        parameter on a per level basis.
    presmoother : {tuple, string, list} : default ('block_gauss_seidel', {'sweep':'symmetric'})
        Defines the presmoother for the multilevel cycling.  The default block
        Gauss-Seidel option defaults to point-wise Gauss-Seidel, if the matrix
        is CSR or is a BSR matrix with blocksize of 1.  See notes below for
        varying this parameter on a per level basis.
    postsmoother : {tuple, string, list}
        Same as presmoother, except defines the postsmoother.
    Bimprove : {list} : default [('block_gauss_seidel', {'sweep':'symmetric'}), None]
        The ith entry defines the method used to improve the candidates B on
        level i.  If the list is shorter than max_levels, then the last entry
        will define the method for all levels lower.
        The list elements are relaxation descriptors of the form used for
        presmoother and postsmoother.  A value of None implies no action on B.
    max_levels : {integer} : default 10
        Maximum number of levels to be used in the multilevel solver.
    max_coarse : {integer} : default 500
        Maximum number of variables permitted on the coarse grid.
    coarsen_diag_dom : {bool, tuple} : default False
        If True (or the first tuple entry is True), then avoid coarsening
        diagonally dominant rows.  The second tuple entry requires a
        dictionary, where the key value 'theta' is used to tune the diagonal
        dominance threshold.
    keep : {bool} : default False
        Flag to indicate keeping extra operators in the hierarchy for
        diagnostics.  For example, if True, then strength of connection (C),
        tentative prolongation (T), aggregation (AggOp), and arrays 
        storing the C-points (Cpts) and F-points (Fpts) are kept at
        each level.

    Other Parameters
    ----------------
    cycle_type : ['V','W','F']
        Structrure of multigrid cycle
    coarse_solver : ['splu', 'lu', 'cholesky, 'pinv', 'gauss_seidel', ... ]
        Solver used at the coarsest level of the MG hierarchy.
            Optionally, may be a tuple (fn, args), where fn is a string such as
        ['splu', 'lu', ...] or a callable function, and args is a dictionary of
        arguments to be passed to fn.

    Returns
    -------
    ml : multilevel_solver
        Multigrid hierarchy of matrices and prolongation operators

    See Also
    --------
    multilevel_solver, aggregation.smoothed_aggregation_solver, 
    classical.ruge_stuben_solver

    Notes
    -----
         - Root-node style SA differs from classical SA primarily by preserving
           and identity block in the interpolation operator, P.  Each aggregate
           has a "root-node" or "center-node" associated with it, and this 
           root-node is injected from the coarse grid to the fine grid.  The 
           injection corresponds to the identity block.

         - Only smooth={'energy', None} is supported for prolongation smoothing.
           See reference [2] below for more details on why the 'energy' prolongation
           smoother is the natural counterpart to root-node style SA.

        - The additional parameters are passed through as arguments to
          multilevel_solver.  Refer to pyamg.multilevel_solver for additional
          documentation.

        - At each level, four steps are executed in order to define the coarser
          level operator.

          1. Matrix A is given and used to derive a strength matrix, C.

          2. Based on the strength matrix, indices are grouped or aggregated.

          3. The aggregates define coarse nodes and a tentative prolongation 
             operator T is defined by injection 

          4. The tentative prolongation operator is smoothed by a relaxation
             scheme to improve the quality and extent of interpolation from the
             aggregates to fine nodes.

        - The parameters smooth, strength, aggregate, presmoother, postsmoother can
          be varied on a per level basis.  For different methods on different
          levels, use a list as input so that the i-th entry defines the method at
          the i-th level.  If there are more levels in the hierarchy than list
          entries, the last entry will define the method for all levels lower.

          Examples are:
          smooth=[('jacobi', {'omega':1.0}), None, 'jacobi']
          presmoother=[('block_gauss_seidel', {'sweep':symmetric}), 'sor']
          aggregate=['standard', 'naive']
          strength=[('symmetric', {'theta':0.25}), ('symmetric',{'theta':0.08})]

        - Predefined strength of connection and aggregation schemes can be
          specified.  These options are best used together, but aggregation can be
          predefined while strength of connection is not.

          For predefined strength of connection, use a list consisting of tuples of
          the form ('predefined', {'C' : C0}), where C0 is a csr_matrix and each
          degree-of-freedom in C0 represents a supernode.  For instance to
          predefine a three-level hierarchy, use [('predefined', {'C' : C0}),
          ('predefined', {'C' : C1}) ].

          Similarly for predefined aggregation, use a list of tuples.  For instance
          to predefine a three-level hierarchy, use [('predefined', {'AggOp' :
          Agg0}), ('predefined', {'AggOp' : Agg1}) ], where the dimensions of A,
          Agg0 and Agg1 are compatible, i.e.  Agg0.shape[1] == A.shape[0] and
          Agg1.shape[1] == Agg0.shape[0].  Each AggOp is a csr_matrix.

          Because this is a root-nodes solver, if a member of the predefined
          aggregation list is predefined, it must be of the form 
          ('predefined', {'AggOp' : Agg, 'Cnodes' : Cnodes}).

    Examples
    --------
    >>> from pyamg import rootnode_solver
    >>> from pyamg.gallery import poisson
    >>> from scipy.sparse.linalg import cg
    >>> import numpy
    >>> A = poisson((100,100), format='csr')           # matrix
    >>> b = numpy.ones((A.shape[0]))                   # RHS
    >>> ml = rootnode_solver(A)                     # AMG solver
    >>> M = ml.aspreconditioner(cycle='V')             # preconditioner
    >>> x,info = cg(A, b, tol=1e-8, maxiter=30, M=M)   # solve with CG

    References
    ----------
    .. [1] Vanek, P. and Mandel, J. and Brezina, M., 
       "Algebraic Multigrid by Smoothed Aggregation for 
       Second and Fourth Order Elliptic Problems", 
       Computing, vol. 56, no. 3, pp. 179--196, 1996.
       http://citeseer.ist.psu.edu/vanek96algebraic.html
    .. [2] Olson, L. and Schroder, J. and Tuminaro, R.,
       "A general interpolation strategy for algebraic 
       multigrid using energy minimization", SIAM Journal
       on Scientific Computing (SISC), vol. 33, pp. 
       966--991, 2011.
    """

    if not (isspmatrix_csr(A) or isspmatrix_bsr(A)):
        try:
            A = csr_matrix(A)
            warn("Implicit conversion of A to CSR", scipy.sparse.SparseEfficiencyWarning)
        except:
            raise TypeError('Argument A must have type csr_matrix or bsr_matrix,\
                             or be convertible to csr_matrix')

    A = A.asfptype()

    if (symmetry != 'symmetric') and (symmetry != 'hermitian') and (symmetry != 'nonsymmetric'):
        raise ValueError('expected \'symmetric\', \'nonsymmetric\' or \'hermitian\' for the symmetry parameter ')
    A.symmetry = symmetry

    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix')
    ##
    # Right near nullspace candidates, use constant for each variable as default
    if B is None:
        B = numpy.kron(numpy.ones((A.shape[0]/blocksize(A),1), dtype=A.dtype), 
                       numpy.eye(blocksize(A)))         
    else:
        B = numpy.asarray(B, dtype=A.dtype)
        if len(B.shape) == 1:
            B = B.reshape(-1,1)
        if B.shape[0] != A.shape[0]:
            raise ValueError('The near null-space modes B have incorrect dimensions for matrix A')
        if B.shape[1] < blocksize(A):
            raise ValueError('B.shape[1] must be >= the blocksize of A') 
    
    ##
    # Left near nullspace candidates
    if A.symmetry == 'nonsymmetric':
        if BH is None:
            BH = B.copy()
        else:
            BH = numpy.asarray(BH, dtype=A.dtype)
            if len(BH.shape) == 1:
                BH = BH.reshape(-1,1)
            if BH.shape[1] != B.shape[1]:
                raise ValueError('The number of left and right near null-space modes,' + \
                                 ' B and BH, must be equal')
            if BH.shape[0] != A.shape[0]:
                raise ValueError('The near null-space modes BH have incorrect dimensions' + \
                                 ' for matrix A')

    ##
    # Preprocess parameters
    max_levels, max_coarse, strength = preprocess_str_or_agg(strength, max_levels, max_coarse)
    max_levels, max_coarse, aggregate = preprocess_str_or_agg(aggregate, max_levels, max_coarse)
    Bimprove = preprocess_Bimprove(Bimprove, A, max_levels)
    smooth = preprocess_smooth(smooth, max_levels)
   
    ##
    # Construct multilevel structure
    levels = []
    levels.append( multilevel_solver.level() )
    levels[-1].A = A          # matrix
   
    ##
    # Append near nullspace candidates
    levels[-1].B = B          # right candidates
    if A.symmetry == 'nonsymmetric':
        levels[-1].BH = BH    # left candidates
    
    while len(levels) < max_levels and levels[-1].A.shape[0]/nPDEs(levels) > max_coarse:
        extend_hierarchy(levels, strength, aggregate, smooth, Bimprove, coarsen_diag_dom, keep)
    
    ml = multilevel_solver(levels, **kwargs)
    change_smoothers(ml, presmoother, postsmoother)
    return ml

def extend_hierarchy(levels, strength, aggregate, smooth, Bimprove, 
                     coarsen_diag_dom=False, keep=True):
    """Service routine to implement the strength of connection, aggregation,
    tentative prolongation construction, and prolongation smoothing.  Called by
    smoothed_aggregation_solver.
    """

    def unpack_arg(v):
        if isinstance(v,tuple):
            return v[0],v[1]
        else:
            return v,{}

    A = levels[-1].A
    B = levels[-1].B
    if A.symmetry == "nonsymmetric":
        AH = A.H.asformat(A.format)
        BH = levels[-1].BH
 
    ##
    # Begin constructing next level
    fn, kwargs = unpack_arg(strength[len(levels)-1])
    if fn == 'symmetric':
        C = symmetric_strength_of_connection(A, **kwargs)
        C = C + eye(C.shape[0], C.shape[1], format='csr')   # Diagonal must be nonzero
    elif fn == 'classical':
        C = classical_strength_of_connection(A, **kwargs)
        C = C + eye(C.shape[0], C.shape[1], format='csr')   # Diagonal must be nonzero
        if isspmatrix_bsr(A):
            C = amalgamate(C, A.blocksize[0])
    elif fn == 'distance':
        C = distance_strength_of_connection(A, **kwargs)
    elif (fn == 'ode') or (fn == 'evolution'):
        C = evolution_strength_of_connection(A, B, **kwargs)
    elif fn == 'energy_based':
        C = energy_based_strength_of_connection(A, **kwargs)
    elif fn == 'predefined':
        C = kwargs['C'].tocsr()
    elif fn is None:
        C = A.tocsr()
    else:
        raise ValueError('unrecognized strength of connection method: %s' % str(fn))
    
    # In SA, strength represents "distance", so we take magnitude of complex values
    if C.dtype == complex:
        C.data = numpy.abs(C.data)
    
    # Create a unified strength framework so that large values represent strong
    # connections and small values represent weak connections
    if (fn == 'ode') or (fn == 'evolution') or (fn == 'distance') or (fn == 'energy_based'):
        C.data = 1.0/C.data

    # Avoid coarsening diagonally dominant rows
    flag,kwargs = unpack_arg( coarsen_diag_dom )
    if flag:
        C = eliminate_diag_dom_nodes(A, C, **kwargs)

    ##
    # aggregation
    fn, kwargs = unpack_arg(aggregate[len(levels)-1])
    if fn == 'standard':
        AggOp,Cnodes = standard_aggregation(C, **kwargs)
    elif fn == 'naive':
        AggOp,Cnodes = naive_aggregation(C, **kwargs)
    elif fn == 'lloyd':
        AggOp,Cnodes = lloyd_aggregation(C, **kwargs)
    elif fn == 'predefined':
        AggOp = kwargs['AggOp'].tocsr()
        Cnodes = kwargs['Cnodes']
    else:
        raise ValueError('unrecognized aggregation method %s' % str(fn))

    ##
    # Improve near nullspace candidates (important to place after the call to
    # evolution_strength_of_connection)
    if Bimprove[len(levels)-1] is not None:
        b = numpy.zeros((A.shape[0],1), dtype=A.dtype)
        B = relaxation_as_linear_operator(Bimprove[len(levels)-1], A, b) * B
        levels[-1].B = B
        if A.symmetry == "nonsymmetric":
            BH = relaxation_as_linear_operator(Bimprove[len(levels)-1], AH, b) * BH 
            levels[-1].BH = BH

    ##
    # tentative prolongator
    T,dummy = fit_candidates(AggOp,B[:,0:blocksize(A)])
    del dummy
    if A.symmetry == "nonsymmetric":
        TH,dummyH = fit_candidates(AggOp,BH[:,0:blocksize(A)])
        del dummyH

    ##
    # Create necessary root node matrices
    Cpt_params = (True, get_Cpt_params(A, Cnodes, AggOp, T))
    T = scale_T(T, Cpt_params[1]['P_I'], Cpt_params[1]['I_F'])
    if A.symmetry == "nonsymmetric":
        TH = scale_T(TH, Cpt_params[1]['P_I'], Cpt_params[1]['I_F'])
        
    ##
    # Set coarse grid modes as injected fine grid modes
    B = Cpt_params[1]['P_I'].T*levels[-1].B
    if A.symmetry == "nonsymmetric":
        BH = Cpt_params[1]['P_I'].T*levels[-1].BH

    ##
    # tentative prolongator smoother
    fn, kwargs = unpack_arg(smooth[len(levels)-1])
    if fn == 'energy':
        P = energy_prolongation_smoother(A, T, C, B, levels[-1].B, Cpt_params=Cpt_params, **kwargs)
    elif fn is None:
        P = T
    else:
        raise ValueError('unrecognized prolongation smoother method %s' % str(fn))

    ##
    # Choice of R reflects A's structure
    symmetry = A.symmetry
    if symmetry == 'hermitian':
        R = P.H
    elif symmetry == 'symmetric':
        R = P.T
    elif symmetry == 'nonsymmetric':
        fn, kwargs = unpack_arg(smooth[len(levels)-1])
        if fn == 'energy':
            R = energy_prolongation_smoother(AH, TH, C, BH, levels[-1].BH, Cpt_params=Cpt_params, **kwargs)
            R = R.H
        elif fn is None:
            R = T.H
        else:
            raise ValueError('unrecognized prolongation smoother method %s' % str(fn))

    if keep:
        levels[-1].C     = C                      # strength of connection matrix
        levels[-1].AggOp = AggOp                  # aggregation operator
        levels[-1].T     = T                      # tentative prolongator
        levels[-1].Fpts  = Cpt_params[1]['Fpts']  # Fpts    
        levels[-1].P_I   = Cpt_params[1]['P_I']   # Injection operator 
        levels[-1].I_F   = Cpt_params[1]['I_F']   # Identity on F-pts
        levels[-1].I_C   = Cpt_params[1]['I_C']   # Identity on C-pts

    levels[-1].P     = P                          # smoothed prolongator
    levels[-1].R     = R                          # restriction operator 
    levels[-1].Cpts  = Cpt_params[1]['Cpts']      # Cpts (i.e., rootnodes) 

    levels.append( multilevel_solver.level() )
    A = R * A * P                                 # Galerkin operator
    A.symmetry = symmetry
    levels[-1].A = A
    levels[-1].B = B                              # right near nullspace candidates

    if A.symmetry == "nonsymmetric":
        levels[-1].BH = BH                        # left near nullspace candidates

