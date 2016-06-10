from setuptools import setup
from mdgenerate import __version__


setup(
    name='mdgenerate',
    description='Collection of python utilities for generation of gromacs simulations.',
    author_email='niels.mueller@physik.tu-darmstadt.de',
    packages=['mdgenerate'],
    version=__version__,
    requires=['numpy', 'mdevaluate', 'yaml', 'jinja2'],
    package_data={'mdgenerate': ['templates/*']},
)
