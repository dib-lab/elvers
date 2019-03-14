from setuptools import setup, find_packages
import glob
import os

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from elvers import __version__, _program


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = _program,
    version = __version__,
    test_suite='pytest.collector',
    tests_require=['pytest'],
    description="snakemake automated workflow system",
    url="https://github.com/dib-lab/elvers",
    author="N.T. Pierce and C. Titus Brown",
    author_email="ntpierce@ucdavis.edu",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points="""
    [console_scripts]
    {program} = {program}.__main__:main
      """.format(program = _program),
    install_requires = required,
    include_package_data=True,
    package_data = { "elvers": ["Snakefile", "*.yaml"] },
    zip_safe=False
)
