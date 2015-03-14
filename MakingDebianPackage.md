General distribution instruction can be found here: http://docs.python.org/distutils/builtdist.html


# Debian baed binaries (.deb) #

## Preliminaries ##

  * install stdeb:
```
apt-get install python-stdeb
```
  * install distutils and setuptools:
```
apt-get install python-distutils-extra python-setuptools
```
  * install debian helper files and python dev files:
```
apt-get install fakeroot python-all-dev debhelper
```

## Procedure ##
(note that dist/pyamg-1.1.0 may be called dist/pyamg-1.0.0.dev700 or similar)

```
cd pyamg
python setup.py sdist
cd dist
py2dsc pyamg-1.1.0.tar.gz
cd deb_dist/pyamg-1.1.0
dpkg-buildpackage -rfakeroot -uc -us
cd ..
```

## Install ##
Then install with
```
sudo dpkg -i python-pyamg_1.0.0-1_i386.deb
```
(note: the name of the .deb file may be different depending on the architecture)

# Fedora binaries (.rpm) #

## Procedure ##
```
python setup.py bdist_rpm
cd dist
```

## Install ##
(to test on Ubuntu)
```
sudo apt-get install alien
alien pyamg-1.1.0-1.i386.rpm
```