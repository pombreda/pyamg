# Prerequisites #

Refer to [Installing](Installing.md) for the list of required packages.  If the required packages have already been installed then skip the corresponding commands below.

This installation process outlined here targets Mac OS X 10.5 Leopard (Intel), but will likely work for 10.4 Tiger as well.

# (Step 0) Essential Compilers and Libraries #

FORTRAN and C compilers are required.  Mac Xcode Tools http://developer.apple.com/TOOLS/xcode/ provides a set of Gnu C/C++ compilers, while http://r.research.att.com/tools/ provides binaries for gfortran (g95).

It is recommended that the stock version of Numpy (1.1.1) be removed prior to the following installation.  It is located in /Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Python.framework/Versions/2.5/Extras/lib/python/numpy

# (Step 1) Download Source Code #

More detailed instructions for NumPy and SciPy are available [here](http://www.scipy.org/Installing_SciPy). Installation info pages:
  * [Nose](http://code.google.com/p/python-nose/)
  * [Numpy](http://www.scipy.org/Installing_SciPy/Mac_OS_X)
  * [Scipy](http://www.scipy.org/Installing_SciPy/Mac_OS_X)
  * [Matplotlib](http://matplotlib.sourceforge.net/faq/installing_faq.html#install-svn)

## (Recommended) Obtain official versions of nose, NumPy, SciPy, and PyAMG ##

```
curl -O http://python-nose.googlecode.com/files/nose-0.10.1.tar.gz
curl -O http://superb-east.dl.sourceforge.net/sourceforge/numpy/numpy-1.2.1.tar.gz
curl -O http://voxel.dl.sourceforge.net/sourceforge/scipy/scipy-0.7.0.tar.gz
curl -O http://pyamg.googlecode.com/files/pyamg-1.0.0.tar.gz
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
```
cd scipy
export MACOSX_DEPLOYMENT_TARGET=10.5
python setup.py build_src build_clib --fcompiler=gnu95 build_ext --fcompiler=gnu95 build
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