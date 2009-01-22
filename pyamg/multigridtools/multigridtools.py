# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.34
#
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _multigridtools
import new
new_instancemethod = new.instancemethod
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


F_NODE = _multigridtools.F_NODE
C_NODE = _multigridtools.C_NODE
U_NODE = _multigridtools.U_NODE

def standard_aggregation(*args):
  """standard_aggregation(int n_row, int Ap, int Aj, int x) -> int"""
  return _multigridtools.standard_aggregation(*args)

def rs_cf_splitting(*args):
  """rs_cf_splitting(int n_nodes, int Sp, int Sj, int Tp, int Tj, int splitting)"""
  return _multigridtools.rs_cf_splitting(*args)

def rs_direct_interpolation_pass1(*args):
  """rs_direct_interpolation_pass1(int n_nodes, int Sp, int Sj, int splitting, int Bp)"""
  return _multigridtools.rs_direct_interpolation_pass1(*args)

def maximal_independent_set_serial(*args):
  """
    maximal_independent_set_serial(int num_rows, int Ap, int Aj, int active, int C, int F, 
        int x) -> int
    """
  return _multigridtools.maximal_independent_set_serial(*args)

def vertex_coloring_mis(*args):
  """vertex_coloring_mis(int num_rows, int Ap, int Aj, int x) -> int"""
  return _multigridtools.vertex_coloring_mis(*args)

def vertex_coloring_jones_plassmann(*args):
  """vertex_coloring_jones_plassmann(int num_rows, int Ap, int Aj, int x, double y) -> int"""
  return _multigridtools.vertex_coloring_jones_plassmann(*args)

def vertex_coloring_LDF(*args):
  """vertex_coloring_LDF(int num_rows, int Ap, int Aj, int x, double y) -> int"""
  return _multigridtools.vertex_coloring_LDF(*args)

def breadth_first_search(*args):
  """breadth_first_search(int Ap, int Aj, int seed, int order, int level)"""
  return _multigridtools.breadth_first_search(*args)

def connected_components(*args):
  """connected_components(int num_nodes, int Ap, int Aj, int components)"""
  return _multigridtools.connected_components(*args)


def rs_direct_interpolation_pass2(*args):
  """
    rs_direct_interpolation_pass2(int n_nodes, int Ap, int Aj, float Ax, int Sp, int Sj, 
        float Sx, int splitting, int Bp, int Bj, 
        float Bx)
    rs_direct_interpolation_pass2(int n_nodes, int Ap, int Aj, double Ax, int Sp, int Sj, 
        double Sx, int splitting, int Bp, int Bj, 
        double Bx)
    """
  return _multigridtools.rs_direct_interpolation_pass2(*args)

def satisfy_constraints_helper(*args):
  """
    satisfy_constraints_helper(int RowsPerBlock, int ColsPerBlock, int num_blocks, 
        int num_block_rows, float x, float y, float z, 
        int Sp, int Sj, float Sx)
    satisfy_constraints_helper(int RowsPerBlock, int ColsPerBlock, int num_blocks, 
        int num_block_rows, double x, double y, double z, 
        int Sp, int Sj, double Sx)
    satisfy_constraints_helper(int RowsPerBlock, int ColsPerBlock, int num_blocks, 
        int num_block_rows, npy_cfloat_wrapper x, npy_cfloat_wrapper y, 
        npy_cfloat_wrapper z, int Sp, 
        int Sj, npy_cfloat_wrapper Sx)
    satisfy_constraints_helper(int RowsPerBlock, int ColsPerBlock, int num_blocks, 
        int num_block_rows, npy_cdouble_wrapper x, npy_cdouble_wrapper y, 
        npy_cdouble_wrapper z, 
        int Sp, int Sj, npy_cdouble_wrapper Sx)
    """
  return _multigridtools.satisfy_constraints_helper(*args)

def invert_BtB(*args):
  """
    invert_BtB(int NullDim, int Nnodes, int ColsPerBlock, float b, 
        int BsqCols, float x, int Sp, int Sj)
    invert_BtB(int NullDim, int Nnodes, int ColsPerBlock, double b, 
        int BsqCols, double x, int Sp, int Sj)
    invert_BtB(int NullDim, int Nnodes, int ColsPerBlock, npy_cfloat_wrapper b, 
        int BsqCols, npy_cfloat_wrapper x, 
        int Sp, int Sj)
    invert_BtB(int NullDim, int Nnodes, int ColsPerBlock, npy_cdouble_wrapper b, 
        int BsqCols, npy_cdouble_wrapper x, 
        int Sp, int Sj)
    """
  return _multigridtools.invert_BtB(*args)

def incomplete_BSRmatmat(*args):
  """
    incomplete_BSRmatmat(int Ap, int Aj, float Ax, int Bp, int Bj, float Bx, 
        int Sp, int Sj, float Sx, int n, int brows, 
        int bcols)
    incomplete_BSRmatmat(int Ap, int Aj, double Ax, int Bp, int Bj, double Bx, 
        int Sp, int Sj, double Sx, int n, int brows, 
        int bcols)
    incomplete_BSRmatmat(int Ap, int Aj, npy_cfloat_wrapper Ax, int Bp, int Bj, 
        npy_cfloat_wrapper Bx, int Sp, int Sj, npy_cfloat_wrapper Sx, 
        int n, int brows, int bcols)
    incomplete_BSRmatmat(int Ap, int Aj, npy_cdouble_wrapper Ax, int Bp, int Bj, 
        npy_cdouble_wrapper Bx, int Sp, int Sj, 
        npy_cdouble_wrapper Sx, int n, int brows, int bcols)
    """
  return _multigridtools.incomplete_BSRmatmat(*args)

def classical_strength_of_connection(*args):
  """
    classical_strength_of_connection(int n_row, float theta, int Ap, int Aj, float Ax, int Sp, 
        int Sj, float Sx)
    classical_strength_of_connection(int n_row, double theta, int Ap, int Aj, double Ax, 
        int Sp, int Sj, double Sx)
    classical_strength_of_connection(int n_row, float theta, int Ap, int Aj, npy_cfloat_wrapper Ax, 
        int Sp, int Sj, npy_cfloat_wrapper Sx)
    classical_strength_of_connection(int n_row, double theta, int Ap, int Aj, npy_cdouble_wrapper Ax, 
        int Sp, int Sj, npy_cdouble_wrapper Sx)
    """
  return _multigridtools.classical_strength_of_connection(*args)

def symmetric_strength_of_connection(*args):
  """
    symmetric_strength_of_connection(int n_row, float theta, int Ap, int Aj, float Ax, int Sp, 
        int Sj, float Sx)
    symmetric_strength_of_connection(int n_row, double theta, int Ap, int Aj, double Ax, 
        int Sp, int Sj, double Sx)
    symmetric_strength_of_connection(int n_row, float theta, int Ap, int Aj, npy_cfloat_wrapper Ax, 
        int Sp, int Sj, npy_cfloat_wrapper Sx)
    symmetric_strength_of_connection(int n_row, double theta, int Ap, int Aj, npy_cdouble_wrapper Ax, 
        int Sp, int Sj, npy_cdouble_wrapper Sx)
    """
  return _multigridtools.symmetric_strength_of_connection(*args)

def apply_distance_filter(*args):
  """
    apply_distance_filter(int n_row, float epsilon, int Sp, int Sj, float Sx)
    apply_distance_filter(int n_row, double epsilon, int Sp, int Sj, double Sx)
    """
  return _multigridtools.apply_distance_filter(*args)

def ode_strength_helper(*args):
  """
    ode_strength_helper(float Sx, int Sp, int Sj, int nrows, float x, float y, 
        float b, int BDBCols, int NullDim)
    ode_strength_helper(double Sx, int Sp, int Sj, int nrows, double x, double y, 
        double b, int BDBCols, int NullDim)
    """
  return _multigridtools.ode_strength_helper(*args)

def min_blocks(*args):
  """
    min_blocks(int n_blocks, int blocksize, float Sx, float Tx)
    min_blocks(int n_blocks, int blocksize, double Sx, double Tx)
    """
  return _multigridtools.min_blocks(*args)

def incomplete_matmat(*args):
  """
    incomplete_matmat(int Ap, int Aj, float Ax, int Bp, int Bj, float Bx, 
        int Sp, int Sj, float Sx, int dimen)
    incomplete_matmat(int Ap, int Aj, double Ax, int Bp, int Bj, double Bx, 
        int Sp, int Sj, double Sx, int dimen)
    """
  return _multigridtools.incomplete_matmat(*args)

def block_gauss_seidel(*args):
  """
    block_gauss_seidel(int Ap, int Aj, float Ax, float x, float b, int row_start, 
        int row_stop, int row_step, int blocksize)
    block_gauss_seidel(int Ap, int Aj, double Ax, double x, double b, int row_start, 
        int row_stop, int row_step, int blocksize)
    block_gauss_seidel(int Ap, int Aj, npy_cfloat_wrapper Ax, npy_cfloat_wrapper x, 
        npy_cfloat_wrapper b, int row_start, 
        int row_stop, int row_step, int blocksize)
    block_gauss_seidel(int Ap, int Aj, npy_cdouble_wrapper Ax, npy_cdouble_wrapper x, 
        npy_cdouble_wrapper b, int row_start, 
        int row_stop, int row_step, int blocksize)
    """
  return _multigridtools.block_gauss_seidel(*args)

def gauss_seidel(*args):
  """
    gauss_seidel(int Ap, int Aj, float Ax, float x, float b, int row_start, 
        int row_stop, int row_step)
    gauss_seidel(int Ap, int Aj, double Ax, double x, double b, int row_start, 
        int row_stop, int row_step)
    gauss_seidel(int Ap, int Aj, npy_cfloat_wrapper Ax, npy_cfloat_wrapper x, 
        npy_cfloat_wrapper b, int row_start, 
        int row_stop, int row_step)
    gauss_seidel(int Ap, int Aj, npy_cdouble_wrapper Ax, npy_cdouble_wrapper x, 
        npy_cdouble_wrapper b, int row_start, 
        int row_stop, int row_step)
    """
  return _multigridtools.gauss_seidel(*args)

def jacobi(*args):
  """
    jacobi(int Ap, int Aj, float Ax, float x, float b, float temp, 
        int row_start, int row_stop, int row_step, 
        float omega)
    jacobi(int Ap, int Aj, double Ax, double x, double b, double temp, 
        int row_start, int row_stop, int row_step, 
        double omega)
    jacobi(int Ap, int Aj, npy_cfloat_wrapper Ax, npy_cfloat_wrapper x, 
        npy_cfloat_wrapper b, npy_cfloat_wrapper temp, 
        int row_start, int row_stop, int row_step, 
        npy_cfloat_wrapper omega)
    jacobi(int Ap, int Aj, npy_cdouble_wrapper Ax, npy_cdouble_wrapper x, 
        npy_cdouble_wrapper b, npy_cdouble_wrapper temp, 
        int row_start, int row_stop, 
        int row_step, npy_cdouble_wrapper omega)
    """
  return _multigridtools.jacobi(*args)

def gauss_seidel_indexed(*args):
  """
    gauss_seidel_indexed(int Ap, int Aj, float Ax, float x, float b, int Id, 
        int row_start, int row_stop, int row_step)
    gauss_seidel_indexed(int Ap, int Aj, double Ax, double x, double b, int Id, 
        int row_start, int row_stop, int row_step)
    """
  return _multigridtools.gauss_seidel_indexed(*args)

def kaczmarz_jacobi(*args):
  """
    kaczmarz_jacobi(int Ap, int Aj, float Ax, float x, float b, float Tx, 
        float temp, int row_start, int row_stop, int row_step, 
        float omega)
    kaczmarz_jacobi(int Ap, int Aj, double Ax, double x, double b, double Tx, 
        double temp, int row_start, int row_stop, 
        int row_step, double omega)
    kaczmarz_jacobi(int Ap, int Aj, npy_cfloat_wrapper Ax, npy_cfloat_wrapper x, 
        npy_cfloat_wrapper b, npy_cfloat_wrapper Tx, 
        npy_cfloat_wrapper temp, int row_start, 
        int row_stop, int row_step, npy_cfloat_wrapper omega)
    kaczmarz_jacobi(int Ap, int Aj, npy_cdouble_wrapper Ax, npy_cdouble_wrapper x, 
        npy_cdouble_wrapper b, npy_cdouble_wrapper Tx, 
        npy_cdouble_wrapper temp, int row_start, 
        int row_stop, int row_step, npy_cdouble_wrapper omega)
    """
  return _multigridtools.kaczmarz_jacobi(*args)

def kaczmarz_gauss_seidel(*args):
  """
    kaczmarz_gauss_seidel(int Ap, int Aj, float Ax, float x, float b, int row_start, 
        int row_stop, int row_step, float Tx)
    kaczmarz_gauss_seidel(int Ap, int Aj, double Ax, double x, double b, int row_start, 
        int row_stop, int row_step, double Tx)
    kaczmarz_gauss_seidel(int Ap, int Aj, npy_cfloat_wrapper Ax, npy_cfloat_wrapper x, 
        npy_cfloat_wrapper b, int row_start, 
        int row_stop, int row_step, npy_cfloat_wrapper Tx)
    kaczmarz_gauss_seidel(int Ap, int Aj, npy_cdouble_wrapper Ax, npy_cdouble_wrapper x, 
        npy_cdouble_wrapper b, int row_start, 
        int row_stop, int row_step, npy_cdouble_wrapper Tx)
    """
  return _multigridtools.kaczmarz_gauss_seidel(*args)

def apply_householders(*args):
  """
    apply_householders(float z, float B, int n, int start, int stop, int step)
    apply_householders(double z, double B, int n, int start, int stop, int step)
    apply_householders(npy_cfloat_wrapper z, npy_cfloat_wrapper B, int n, 
        int start, int stop, int step)
    apply_householders(npy_cdouble_wrapper z, npy_cdouble_wrapper B, int n, 
        int start, int stop, int step)
    """
  return _multigridtools.apply_householders(*args)

def householder_hornerscheme(*args):
  """
    householder_hornerscheme(float z, float B, float y, int n, int start, int stop, 
        int step)
    householder_hornerscheme(double z, double B, double y, int n, int start, int stop, 
        int step)
    householder_hornerscheme(npy_cfloat_wrapper z, npy_cfloat_wrapper B, npy_cfloat_wrapper y, 
        int n, int start, int stop, int step)
    householder_hornerscheme(npy_cdouble_wrapper z, npy_cdouble_wrapper B, npy_cdouble_wrapper y, 
        int n, int start, int stop, 
        int step)
    """
  return _multigridtools.householder_hornerscheme(*args)

def apply_givens(*args):
  """
    apply_givens(float B, float x, int n, int nrot)
    apply_givens(double B, double x, int n, int nrot)
    apply_givens(npy_cfloat_wrapper B, npy_cfloat_wrapper x, int n, 
        int nrot)
    apply_givens(npy_cdouble_wrapper B, npy_cdouble_wrapper x, int n, 
        int nrot)
    """
  return _multigridtools.apply_givens(*args)

def maximal_independent_set_parallel(*args):
  """
    maximal_independent_set_parallel(int num_rows, int Ap, int Aj, int active, int C, int F, 
        int x, double y, int max_iters=-1) -> int
    maximal_independent_set_parallel(int num_rows, int Ap, int Aj, int active, int C, int F, 
        int x, double y) -> int
    """
  return _multigridtools.maximal_independent_set_parallel(*args)

def maximal_independent_set_k_parallel(*args):
  """
    maximal_independent_set_k_parallel(int num_rows, int Ap, int Aj, int k, int x, double y, 
        int max_iters=-1)
    maximal_independent_set_k_parallel(int num_rows, int Ap, int Aj, int k, int x, double y)
    """
  return _multigridtools.maximal_independent_set_k_parallel(*args)

def bellman_ford(*args):
  """
    bellman_ford(int num_rows, int Ap, int Aj, int Ax, int x, int y)
    bellman_ford(int num_rows, int Ap, int Aj, float Ax, float x, int y)
    bellman_ford(int num_rows, int Ap, int Aj, double Ax, double x, 
        int y)
    """
  return _multigridtools.bellman_ford(*args)

def lloyd_cluster(*args):
  """
    lloyd_cluster(int num_rows, int Ap, int Aj, int Ax, int num_seeds, 
        int x, int y, int z)
    lloyd_cluster(int num_rows, int Ap, int Aj, float Ax, int num_seeds, 
        float x, int y, int z)
    lloyd_cluster(int num_rows, int Ap, int Aj, double Ax, int num_seeds, 
        double x, int y, int z)
    """
  return _multigridtools.lloyd_cluster(*args)

def fit_candidates(*args):
  """
    fit_candidates(int n_row, int n_col, int K1, int K2, int Ap, int Ai, 
        float Ax, float B, float R, float tol)
    fit_candidates(int n_row, int n_col, int K1, int K2, int Ap, int Ai, 
        double Ax, double B, double R, double tol)
    fit_candidates(int n_row, int n_col, int K1, int K2, int Ap, int Ai, 
        npy_cfloat_wrapper Ax, npy_cfloat_wrapper B, 
        npy_cfloat_wrapper R, float tol)
    fit_candidates(int n_row, int n_col, int K1, int K2, int Ap, int Ai, 
        npy_cdouble_wrapper Ax, npy_cdouble_wrapper B, 
        npy_cdouble_wrapper R, double tol)
    """
  return _multigridtools.fit_candidates(*args)

