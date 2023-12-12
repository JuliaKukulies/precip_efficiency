import os
import setuptools 
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "microphysics",
    version = "0.0.1",
    author = "Julia Kukulies",
    author_email = "kukulies@ucar.edu",
    description = ("A package for the analysis of microphysical processes in model output."),
    license = "BSD-3-Clause License", 
    packages = ["microphysics"],
    long_description=read('README.md'),
    zip_safe = False,
)
