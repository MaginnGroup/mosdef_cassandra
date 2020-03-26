from setuptools import setup, find_packages

#####################################
VERSION = "0.1.1"
ISRELEASED = True
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + ".dev0"
#####################################

requirements = [
    "numpy",
    "parmed",
    "networkx",
    "mbuild >=0.10.5",
]

setup(
    name="mosdef_cassandra",
    version=__version__,
    packages=find_packages(),
    license="MIT",
    author="Ryan S. DeFever",
    author_email="rdefever@nd.edu",
    url="https://github.com/MaginnGroup/mosdef_cassandra",
    install_requires=requirements,
    python_requires=">=3.6, <4",
)
