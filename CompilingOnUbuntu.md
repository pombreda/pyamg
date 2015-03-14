# Prerequisites #

Refer to [Installing](Installing.md) for the list of required packages.  If the required packages have already been installed then skip the corresponding commands below.

# (Step 0) Essential Compilers and Libraries #

The instructions below assume that the following Ubuntu packages have been installed:
```
sudo apt-get install build-essential gfortran libatlas-sse2-dev 
sudo apt-get install python-all-dev ipython
sudo apt-get install subversion
```

When using other Fortran compilers (e.g. f77) or BLAS/LAPACK libraries, ensure that the combination is compatible.


# (Step 1) Download Source Code #

More detailed instructions for NumPy and SciPy are available [here](http://www.scipy.org/Installing_SciPy).

## (Recommended) Obtain official versions of nose, NumPy, SciPy, and PyAMG ##

```
wget http://python-nose.googlecode.com/files/nose-0.10.1.tar.gz
wget http://superb-east.dl.sourceforge.net/sourceforge/numpy/numpy-1.2.1.tar.gz
wget http://voxel.dl.sourceforge.net/sourceforge/scipy/scipy-0.7.0.tar.gz
wget http://pyamg.googlecode.com/files/pyamg-1.0.0.tar.gz
tar xvfz nose-0.10.1.tar.gz
tar xvfz numpy-1.2.1.tar.gz
tar xvfz scipy-0.7.0.tar.gz
tar xvfz pyamg-1.0.0.tar.gz
```

Alternative Links: [Nose 0.10.1](http://code.google.com/p/python-nose/downloads/list)
[Numpy 1.2.1](http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103)
[SciPy 0.7](http://sourceforge.net/project/showfiles.php?group_id=27747&package_id=19531&release_id=660191)
[PyAMG 1.0](http://code.google.com/p/pyamg/downloads/list)

## (Alternative) Obtain development versions of nose, NumPy, SciPy, and PyAMG ##
```
svn checkout http://python-nose.googlecode.com/svn/trunk/ nose
svn checkout http://svn.scipy.org/svn/numpy/trunk/ numpy
svn checkout http://svn.scipy.org/svn/scipy/trunk/ scipy
svn checkout http://pyamg.googlecode.com/svn/trunk/ pyamg
```



# (Step 2) Install Packages #

## Install nose ##
```
cd nose
sudo python setup.py install
cd ..
```

## Install NumPy ##
```
cd numpy
python setup.py build
sudo python setup.py install
cd ..
```

## Install SciPy ##
If you have umfpack installed from apt-get install suitesparse suitesparse-dev, put the following site.cfg file in your scipy directory:
```
[amd]
library_dirs = /usr/lib
include_dirs = /usr/include/suitesparse
amd_libs = amd

[umfpack]
library_dirs = /usr/lib
include_dirs = /usr/include/suitesparse
umfpack_libs = umfpack
```

Then follow with
```
cd scipy
python setup.py build
sudo python setup.py install
cd ..
```

## Install PyAMG ##
```
cd pyamg
python setup.py build
sudo python setup.py install
```



# (Step 3) Test Installation #

Using Python or IPython (recommended) to run the following commands

```
import pyamg
pyamg.test()
```

should have output similar to

```
Running unit tests for pyamg
NumPy version 1.2.1
NumPy is installed in /usr/lib/python2.5/site-packages/numpy
SciPy version 0.7.0
SciPy is installed in /usr/lib/python2.5/site-packages/scipy
Python version 2.5.2 (r252:60911, Jul 31 2008, 17:31:22) [GCC 4.2.3 (Ubuntu 4.2.3-2ubuntu7)]
nose version 0.10.4
PyAMG version 1.0.0
PyAMG is installed in /usr/lib/python2.5/site-packages/pyamg
..........................................................................................................................
----------------------------------------------------------------------
Ran 122 tests in 27.593s

OK
```


Testing NumPy and SciPy is done in similar fashion:
```
import numpy
numpy.test()

import scipy
scipy.test()
```