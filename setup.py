from setuptools import setup
from codecs import open  ## To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
    
setup(
    name="jcamp",
    py_modules=["jcamp"],
    version="1.1",
    description="JCAMP-DX file reader",
    long_description = long_description,         ## from above
    author="Nathan Hagen",
    author_email="nhagen@optics.arizona.edu",
    url="https://github.com/nzhagen/jcamp",
    download_url="https://github.com/nzhagen/jcamp",
    install_requires=['numpy'],
    test_suite='tests',
    tests_require=['tox', 'coverage'],
    keywords=["jcamp", "jcamp-dx", "spectra"],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        ],
    include_package_data=True
)
