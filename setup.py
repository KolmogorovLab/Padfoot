import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from padfoot.__version__ import __version__


setup(name='padfoot',
      version=__version__,
      description='Annotation tool for Severus and Wakhan',
      url='https://github.com/KolmogorovLab/Padfoot',
      author='Ayse Keskus',
      author_email = 'aysegokce.keskus@nih.gov',
      license='BSD-3-Clause',
      packages=['python'],
      package_data={'python': ['beds/*']},
      entry_points={'console_scripts': ['padfoot = padfoot.main:main']},
      )
