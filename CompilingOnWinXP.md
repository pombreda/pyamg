# Prerequisites #

Refer to [Installing](Installing.md) for the list of required packages.  If the required packages have already been installed then skip the corresponding commands below.

# (Step 0) Essential Compilers and Libraries #

## Install MinGW ##

Download the [Automated MinGW Installer](http://sourceforge.net/project/showfiles.php?group_id=2435).  Run the installer and ensure that compilers for C, C++, and Fortran are selected.

Now add the directory where the MinGW compilers are stored (usually C:\MinGW\bin)to your PATH.

Permanent method: Control Panel->System->Advanced->Environment Variables->System Variables->Path

Temporary method: From the command line run: `set PATH=c:\mingw\bin\;%PATH%`

## Install Python 2.5 ##
Install Python 2.5 using the official [installer](http://www.python.org/ftp/python/2.5.2/python-2.5.2.msi).


## LAPACK ##
[LAPACK](http://www.netlib.org/lapack/lapack-lite-3.1.1.tgz) (including BLAS)

Extract LAPACK to LIBS\lapack-3.1.1 (or whatever the current version is)

Set environment variables for subsequent steps
```
set BLAS_SRC=..\lapack-3.1.1\BLAS\SRC
set LAPACK_SRC=..\lapack-3.1.1\SRC
```


## Install Subversion ##
We'll need Subversion to download the latest versions of NumPy, SciPy, and PyAMG.  There are [several](http://subversion.tigris.org/getting.html#windows) binary distributions available.  The one provided by [Silk SVN](http://www.sliksvn.com/en/download) is easy to install.


# (Step 1) Download Source Code #

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
python setup.py install
cd ..
```

## Install NumPy ##
```
cd numpy
python setup.py config --compiler=mingw32 build --compiler=mingw32 install
cd ..
```

> Note: Simply rerun the previous command a few times if the following error occurs.
```
g77.exe: C:/MinGW/bin/../lib/gcc/mingw32/3.4.5/specs: Too many open files
error: Command "C:\MinGW\bin\g77.exe <some arguments> failed with exit status 1
Try rerunning setup command until build succeeds.
```


## Install SciPy ##
```
cd scipy
python setup.py config --compiler=mingw32 build --compiler=mingw32 install
cd ..
```

## Install PyAMG ##
```
cd pyamg
python setup.py config --compiler=mingw32 build --compiler=mingw32 install
cd..
```

(Optional) Build binary installer with:
```
python setup.py  build --compiler=mingw32 bdist_wininst
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