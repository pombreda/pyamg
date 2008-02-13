"Matrix Gallery for Multigrid Solvers"

from info import __doc__

from laplacian import *
from elasticity import *
from example import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from scipy.testing.pkgtester import Tester
test = Tester().test
