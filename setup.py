from setuptools import setup, find_packages

from elvers import __version__, _program

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = _program,
    version = __version__,
    description="automated RNAseq analysis",
    url="https://github.com/bluegenes/elvers",
    author="N. Tessa Pierce and C. Titus Brown",
    author_email="ntpierce@gmail.com",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points = {'console_scripts': [
        'elvers  = elvers.__main__:main'
        ]
    },
    include_package_data=True,
    package_data = { "elvers": ["elvers.snakefile", "*.yaml", "*.yml"] },
    setup_requires = [ "setuptools>=38.6.0",
                       'setuptools_scm', 'setuptools_scm_git_archive' ],
    use_scm_version = {"write_to": "elvers/version.py"},
    install_requires = ['snakemake>=5.26', 'click>=7', 'pandas>=1.1.2']
)
