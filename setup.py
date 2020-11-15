"""openmoltools: Tools for Small Molecules, Antechamber, OpenMM, and More.
"""

DOCLINES = __doc__.split("\n")

from setuptools import setup

##########################
VERSION = "0.0.0dev0"
ISRELEASED = False
__version__ = VERSION
##########################


CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

extensions = []

setup(
    name="openmoltools",
    author="Kyle A. Beauchamp",
    author_email="kyleabeauchamp@gmail.com",
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=__version__,
    license="MIT",
    url="http://github.com/choderalab/openmoltools",
    platforms=["Linux", "Mac OS-X", "Unix"],
    classifiers=CLASSIFIERS.splitlines(),
    packages=["openmoltools", "openmoltools.tests"],
    zip_safe=False,
    scripts=["scripts/generate_example_data.py", "scripts/processAmberForceField.py"],
    ext_modules=extensions,
    # Install all data directories of the form testsystems/data/X/
    package_data={"openmoltools": ["chemicals/*/*", "parameters/*"]},
)
