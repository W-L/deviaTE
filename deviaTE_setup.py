from setuptools import setup

setup(
    name='deviaTE',
    version='0.1',
    description='polymorphic TEs',
    license='BSD',
    author='Lukas W',
    # author_email='foomail@foo.com',
    # url="http://www.foopackage.com/",
    packages=['bin'],  # same as name
    install_requires=['pysam', 'pandas'],  # subprocess, itertools
    scripts=['deviate_analyse.py']
)
