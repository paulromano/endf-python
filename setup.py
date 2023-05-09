from setuptools import setup

from pybind11.setup_helpers import build_ext, intree_extensions


ext_modules = intree_extensions(["src/endf/_records.cpp"])

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
