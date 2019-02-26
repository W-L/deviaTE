from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='deviaTE',
    version='0.3.3',
    description='Analysis and visualization of mobile genetic element composition',
    long_description=readme(),
    license='LICENSE',
    author='Lukas W',
    url="https://github.com/W-L/deviaTE.git",
    packages=['deviaTE'],
    include_package_data=True,
    install_requires=['pysam', 'pandas'],
    scripts=['bin/deviaTE_analyse',
             'bin/deviaTE_plot',
             'bin/deviaTE_prep',
             'bin/deviaTE_fuse',
             'bin/deviaTE',
             'bin/deviaTE_trim.pl',
             'bin/deviaTE_trim.pm']
)
