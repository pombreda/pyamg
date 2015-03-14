On Linux systems we usually install PyAMG with
```
cd pyamg
python setup.py build
sudo python setup.py install
```
which places PyAMG in the standard site-packages directory (e.g. `/usr/local/lib/python2.5/site-packages/pyamg`).  However, if you lack super-user access on your system (i.e. you cannot `sudo`) then you'll need to install PyAMG to a local directory.  For example
```
mkdir $HOME/.local
cd pyamg
python setup.py build
python setup.py install --home=$HOME/.local
export PYTHONPATH=$HOME/.local/lib/python2.6/site-packages/
```
installs PyAMG to a `$HOME/.local/lib/python2.6/site-packages/` and adds the local site-packages directory to the list of directories to search for python packages.

See http://docs.python.org/install/index.html#alternate-installation-the-home-scheme for more details.