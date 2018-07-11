#!/usr/bin/env python
from setuptools import setup
import numpy as np

setup(name='europa_seismo',
      version='0.1',
      description='Package for seismic modeling of Europa',
      author='Ross Maguire',
      author_email='rmaguire@umd.edu',
      url='www.github.com/romaguir/europa_seismo',
      packages=['europa_seismo'],
      install_requires=['pymc'],
      scripts=['bin/find_accel_max'],
      license='GNU')
