#!/usr/bin/env python
import numpy as np
from numpy.distutils.core import setup,Extension
from setuptools import find_packages

minos = Extension(name = 'europa_seismo.minos', sources = ['src/minos_bran.f'])
rayleigh_python = Extension(name='europa_seismo.rayleigh_python',sources=['src/rayleigh_python.f', \
                            'src/ltable_python.f'])

setup(name='europa_seismo',
      version='0.1',
      description='Package for seismic modeling of Europa',
      author='Ross Maguire',
      author_email='rmaguire@umd.edu',
      url='www.github.com/romaguir/europa_seismo',
      packages=find_packages(),
      install_requires=['pymc'],
      scripts=['bin/find_accel_max'],
      license='GNU',
      ext_modules = [minos,rayleigh_python])
