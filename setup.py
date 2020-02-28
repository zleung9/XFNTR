# All .ui files and .so files are added through keyword: package_data, because setuptools doesn't include them automatically.
import sys
import os
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = "xfntr",
    version = "0.2.1",
    author = "Zhu Liang",
    author_email = "zliang8@uic.edu",
    description = "A software that analyzes xfntr data",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/zhul9311/XFNTR.git",
    packages = find_packages(),
    package_dir = {'':'.'},
    package_data = {
        '' : ['xr_ref.cpython-37m-darwin.so','GUI/*']
    },
    exclude_package_data = {
        '' : ['.git/','.setup.py.swp']
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires = '>=3.6',
    install_requires = [
        'pyqt5',
        'scipy',
        'matplotlib',
        'lmfit',
        'periodictable'
    ],
    entry_points = { # create scripts and add to sys.PATH
        'console_scripts':[
            'xfntr1 = xfntr.main:main'
        ],
        'gui_scripts': [
            'xfntr = xfntr.main:main'
        ]
    },
)
