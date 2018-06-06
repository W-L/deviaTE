from setuptools import setup

setup(
    name='deviaTE',
    version='0.1',
    description='polymorphic TEs',
    license='BSD',
    author='Lukas W',
    url="https://github.com/W-L/deviaTE.git",
    packages=['bin'],
    install_requires=['pysam', 'pandas'],
    scripts=['deviaTE_analyse.py',
    		 'deviaTE_plot.R',
    		 'deviaTE_prep.sh',
    		 'deviaTE.py']
)
