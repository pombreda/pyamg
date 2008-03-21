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

def sa_get_aggregates(*args):
  """sa_get_aggregates(int n_row, int Ap, int Aj, int x) -> int"""
  return _multigridtools.sa_get_aggregates(*args)

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


def rs_strong_connections(*args):
  """
    rs_strong_connections(int n_row, float theta, int Ap, int Aj, float Ax, int Sp, 
        int Sj, float Sx)
    rs_strong_connections(int n_row, double theta, int Ap, int Aj, double Ax, 
        int Sp, int Sj, double Sx)
    """
  return _multigridtools.rs_strong_connections(*args)

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

def sa_strong_connections(*args):
  """
    sa_strong_connections(int n_row, float epsilon, int Ap, int Aj, float Ax, 
        int Sp, int Sj, float Sx)
    sa_strong_connections(int n_row, double epsilon, int Ap, int Aj, double Ax, 
        int Sp, int Sj, double Sx)
    """
  return _multigridtools.sa_strong_connections(*args)

def block_gauss_seidel(*args):
  """
    block_gauss_seidel(int Ap, int Aj, float Ax, float x, float b, int row_start, 
        int row_stop, int row_step, int blocksize)
    block_gauss_seidel(int Ap, int Aj, double Ax, double x, double b, int row_start, 
        int row_stop, int row_step, int blocksize)
    """
  return _multigridtools.block_gauss_seidel(*args)

def gauss_seidel(*args):
  """
    gauss_seidel(int Ap, int Aj, float Ax, float x, float b, int row_start, 
        int row_stop, int row_step)
    gauss_seidel(int Ap, int Aj, double Ax, double x, double b, int row_start, 
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

def maximal_independent_set_parallel(*args):
  """
    maximal_independent_set_parallel(int num_rows, int Ap, int Aj, int active, int C, int F, 
        int x, double y, int max_iters=-1) -> int
    maximal_independent_set_parallel(int num_rows, int Ap, int Aj, int active, int C, int F, 
        int x, double y) -> int
    """
  return _multigridtools.maximal_independent_set_parallel(*args)

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

