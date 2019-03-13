from setuptools import setup, find_packages

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
    name = 'elvers',
    version = "0.1",
    description="snakemake automated workflow system",
    url="https://github.com/dib-lab/eelpond",
    author="N.T. Pierce and C. Titus Brown",
    author_email="ntpierce@ucdavis.edu",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points = {'console_scripts': [
        'elvers  = elvers.__main__:main'
        ]
    },
    include_package_data=True,
    package_data = { "elvers": ["Snakefile", "*.yaml"] },
    install_requires = [
        'yaml',
        'snakemake',
        'numpy',
        'pandas',
        'pytest',
        'graphviz',
        'networkx',
        'pygraphviz',
        'psutil',]
)
