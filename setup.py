from os import path
from setuptools import setup, find_packages
from hotspots_framework import __version__


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()


setup(
    name='hotspots_framework',
    version=__version__,
    description='BBGLab tool',
    packages=find_packages(),
    install_requires=required
)