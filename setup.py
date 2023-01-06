"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject

pip releasing a next version:
1. vim pcDataLoader/VERSION.py  # increase version number in file
2. git add pcDataLoader/VERSION.py
3. git commit -m'@ pcDataLoader : next release.'
4. git tag -a v0.0.0 -m'version 0.0.0'
5. python3 -m build --sdist  # make source distribution
6. python3 -m build --wheel  # make binary distribution python wheel
7. twine upload dist/* --verbose
8. git push origin master
9. git push --tag
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Get the version number from the VERSION.py file
exec(open('./pcDataLoader/VERSION.py').read())

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # the basics
    name="pcDataLoader",  # Required
    version=__version__,  # Required

    # description
    description="pcDataLoader provides a platform independent, python3 based, pip installable interface to load output, generated with the PhysiCell agent based modeling framework, into python3.",  # Optional
    long_description=long_description,  # Optional
    long_description_content_type="text/markdown",  # Optional (see note above)

    # the project's main homepage.
    url="https://github.com/elmbeech/pcDataLoader",  # Optional
    author="Elmar Bucher",  # Optional
    author_email="epbucher@iu.edu",  # Optional
    #author_email="heiland@iu.edu",  # Optional

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish
        "License :: OSI Approved :: BSD License",

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3 :: Only",
    ],

    # What does your project relate to?
    keywords="physicell, python3, data analysis",  # Optional

    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    #package_dir={"": "src"},  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(where=".", exclude=(), include=('*',)),  # Required


    # Specify which Python versions you support.
    python_requires=">=3.6, <4",

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/discussions/install-requires-vs-requirements/
    install_requires=["anndata", "numpy", "pandas", "scipy"],  # Optional

    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    #extras_require={  # Optional
    #    "dev": ["check-manifest"],
    #    "test": ["coverage"],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.
    package_data={  # Optional
        "pcDataLoader": [
            "data_timeseries_2d/*.mat",
            "data_timeseries_2d/*.svg",
            "data_timeseries_2d/*.txt",
            "data_timeseries_2d/*.xml"
            "data_timeseries_3d/*.mat",
            "data_timeseries_3d/*.svg",
            "data_timeseries_3d/*.txt",
            "data_timeseries_3d/*.xml"
        ],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/distutils/setupscript.html#installing-additional-files
    #
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[("my_data", ["data/data_file"])],  # Optional

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    #
    # For example, the following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    #entry_points={  # Optional
    #    "console_scripts": [
    #        "sample=sample:main",
    #    ],
    #},

    # List additional URLs that are relevant to your project as a dict.
    #
    # This field corresponds to the "Project-URL" metadata fields:
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    #
    # Examples listed include a pattern for specifying where the package tracks
    # issues, where the source is hosted, where to say thanks to the package
    # maintainers, and where to support the project financially. The key is
    # what's used to render the link text on PyPI.
    project_urls={  # Optional
        "Bug Reports": "https://github.com/elmbeech/pcDataLoader/issues",
        "Funding": "http://www.mathcancer.org/",
        "Say Thanks!": "http://physicell.org/",
        "Source": "https://github.com/elmbeech/pcDataLoader",
    },
)
