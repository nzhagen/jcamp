# chardet's setup.py
from distutils.core import setup

setup(
    name = "jcamp",
    py_modules = ["jcamp"],
    version = "1.0",
    description = "JCAMP-DX file reader",
    author = "Nathan Hagen",
    author_email = "nhagen@optics.arizona.edu",
    url = "https://github.com/nzhagen/jcamp",
    download_url = "https://github.com/nzhagen/jcamp",
    keywords = ["jcamp", "jcamp-dx", "spectra"],
    classifiers = [
        "Programming Language :: Python :: 2.7",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        ],
    long_description = """\
A reader for JCAMP-DX format spectral data files.
-----------------------------------------------------------

"jcamp" reads in any JCAMP-DX format file and builds a dictionary containing the file
header metadata. The spectrum itself is saved as a pair of Numpy arrays, accessed as
keys "x" and "y" in the returned dictionary.

This version requires Python 2.7 or later.
"""
)
