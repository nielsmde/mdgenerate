from setuptools import setup


def get_version(module):
    version = ''
    with open(module) as f:
        for line in f:
            if '__version__' in line:
                version = line.split('=')[-1].strip("' \n\t")
                break
    return version


__version__ = get_version('mdgenerate/__init__.py')


setup(
    name='mdgenerate',
    description='Collection of python utilities for generation of gromacs simulations.',
    author_email='niels.mueller@physik.tu-darmstadt.de',
    packages=['mdgenerate'],
    version=__version__,
    requires=['numpy', 'mdevaluate', 'yaml', 'jinja2'],
    package_data={'mdgenerate': ['templates/*']},
    entry_points={
        'console_scripts': [
            'mdprocess = mdgenerate.cli:run',
            'npt2nvt = mdgenerate.cli:generate_nvt',
            'makeshort = mdgenerate.cli:make_short'
        ]
    },
    zip_safe=False,
)
