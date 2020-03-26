from setuptools import setup, find_packages

#####################################
VERSION = "0.1.1"
ISRELEASED = False
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################

requirements = [
    'numpy',
    'parmed',
    'networkx',
    'mbuild >= 0.10.5',
    'cassandra >= 1.2.2',
]

setup(
    name='gmso',
    version=__version__,
    packages=find_packages(),
    zip_safe=True,
    author='Matthew W Thompson, Justin Gilmer',
    author_email='matt.thompson@vanderbilt.edu, justin.b.gilmer@vanderbilt.edu',
    url='https://github.com/mosdef-hub/gmso',
    download_url='https://github.com/mosdef-hub/gmso/tarball/{}'.format(__version__),
    license="MIT",
    keywords='gmso',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
    install_requires=requirements,
    python_requires='>=3.6, <4',
)

setup(
    name="mosdef_cassandra",
    version=__version__,
    packages=find_packages(),
    license="MIT",
    author="Ryan S. DeFever",
    author_email="rdefever@nd.edu",
    url="https://github.com/MaginnGroup/mosdef_cassandra",
    install_requires=requirements,
    python_requires='>=3.6, <4',
)
