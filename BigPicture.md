(This page is in construction)

# Where PyAMG fits in #

  * PyAMG solves large linear systems with algebraic multigrid methods.
  * PyAMG benefits from a wealth of scientific libraries in Python
    * SciPy provides sparse linear algebra support (`scipy.sparse`) and IO
    * NumPy provides dense linear algebra (`numpy.linalg`) and sophisticated array processing facilities
    * Matplotlib provides extensive plotting functionality
  * Conformance to Python conventions and modularity of design allows libraries, like PyAMG, to leverage these components without creating fragile dependencies


<table align='center' border='0'>
<tr><td>
<img src='http://pyamg.googlecode.com/svn/wiki/Images/SoftwareEcosystem.png' />
</td></tr>
</table>




# Application Anatomy #

<table align='center' border='0'>
<tr><td>
<img src='http://pyamg.googlecode.com/svn/wiki/Images/ApplicationAnatomy.png' />
</td></tr></table>


# Hybrid Language Approach with Python #

  * Python itself is not designed for HPC
  * However, performance-sensitive portions of most applications are a small fraction of the total codebase (Pareto principle)
    * E.g. IO, handling configuration files, passing pointers around
    * Illustrate the point with some function in PyAMG
  * Suggests hybrid strategy
    * [80%] Flexible and friendly language (Python) by default
    * [20%] Natively compiled language (C\C++\Fortran) for performance sensitive parts
  * Python is good "glue"
    * Interfaces with C/C++/Fortran and many others
    * Freely available and widely deployed
  * Python comes with PIL, XML, email, networking protocols and file formats, etc.


# Related Software #

## Core Libraries ##
  * NumPy
  * SciPy

## Finite Elements ##
  * FiPy
  * Dolfin
  * sfepy

## Graphs ##
  * Networks
  * PyMETIS

## Visualization ##
  * Matplotlib

# Further Reading #