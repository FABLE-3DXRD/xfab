Welcome to xfab 
==================
 ***- a python package for crystallographic computations!***
________________________________________________________

[![macos-latest, setup.py install and pytest](https://github.com/FABLE-3DXRD/xfab/actions/workflows/install-and-test-macos-py39.yml/badge.svg)](https://github.com/FABLE-3DXRD/xfab/actions/workflows/install-and-test-macos-py39.yml)
[![ubuntu-latest, setup.py install and pytest](https://github.com/FABLE-3DXRD/xfab/actions/workflows/install-and-test-ubuntu-py39.yml/badge.svg)](https://github.com/FABLE-3DXRD/xfab/actions/workflows/install-and-test-ubuntu-py39.yml)
[![On pypi release, ubuntu-latest, pip install and pytest](https://github.com/FABLE-3DXRD/xfab/actions/workflows/verify-pypi-release-ubuntu-py39.yml/badge.svg)](https://github.com/FABLE-3DXRD/xfab/actions/workflows/verify-pypi-release-ubuntu-py39.yml)
[![On pypi release, macos-latest, pip install and pytest](https://github.com/FABLE-3DXRD/xfab/actions/workflows/verify-pypi-release-macos-py39.yml/badge.svg)](https://github.com/FABLE-3DXRD/xfab/actions/workflows/verify-pypi-release-macos-py39.yml)
________________________________________________________

xfab was originally a python library for the totalcryst program with crystallographic computations inside. Mostly written by Henning Sorensen and Jette Oddershede at Riso.

Today xfab is maintaned by the [FABLE-3DXRD organization](https://github.com/FABLE-3DXRD) which features many open source libraries for crystallography and diffraction analysis.

Installation
----------------------------
For general use you may want to get the [latest release from pypi](https://pypi.org/project/xfab/) using pip
    
    pip install xfab
 
 If you want the smokingly hot and fresh current code, you could clone [the repo](https://github.com/FABLE-3DXRD/xfab) and build locally from the source
 
     git clone https://github.com/FABLE-3DXRD/xfab.git
     cd xfab
     pip install .

If you want to make a contribution to xfab you will want to make a fork, clone it, and install in editable mode, using

     pip install --editable .
