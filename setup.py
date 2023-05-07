from setuptools import setup, find_packages

from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.1.0"

ext_modules = [Pybind11Extension('endfpy._records', ["src/endfpy/_records.cpp"])]

setup(
    name="endfpy",
    version=__version__,
    author="Paul Romano",
    author_email="paul.k.romano@gmail.com",
    url="https://github.com/pybind/python_example",
    description="ENDF Parser",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
