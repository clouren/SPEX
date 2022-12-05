from setuptools import find_packages, setup

setup(
    name='SPEX',
    packages=find_packages(include=['backslash', 'utilities', '_connect']),
    version='0.1.0',
    description='Python interface for SPEX',
    author='Lorena Mejia Domenzain',
    license='GNU',
)

