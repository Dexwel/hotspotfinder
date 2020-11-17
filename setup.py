from os import path
from setuptools import setup, find_packages
from hotspots_framework import __version__


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

with open(path.join(directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='hotspots_framework',
    version=__version__,
    description='BBGLab tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=required,
    url="https://bitbucket.org/carnedo/hotspots_framework",
    author="Claudia Arnedo",
    author_email="bbglab@irbbarcelona.org",
)