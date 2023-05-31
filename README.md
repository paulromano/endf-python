# ENDF Python Interface

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/pypi/v/endf?label=PyPI)](https://pypi.org/project/endf)

`endf` is a Python package for reading and interpreting
[ENDF-6](https://doi.org/10.2172/1425114) and
[ACE](https://github.com/NuclearData/ACEFormat) format nuclear data files.
Compared to other packages that provide functionality for working with ENDF and
ACE files, this package has numerous advantages:

- Easily installable through `pip`
- Thoughtful API design targeting Python first
- Offers both a low-level interface for working with the raw data in an ENDF
  file as well as a more intuitive high-level interface
- Fast file-loading performance thanks to optimized read routines
- Fully documented, tested, and type-hinted

## Installation

The `endf` Python package can be installed from the [Python Package
Index](https://pypi.org/project/endf/) using `pip`:

```sh
python -m pip install endf
```

## Documentation

For information on how to use the `endf` package and a full API reference,
please see the [online documentation](https://endf-python.readthedocs.io).
