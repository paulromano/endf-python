# ENDF Python Interface

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

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

The `endf` Python package can be install from the [Python Package
Index](https://pypi.org/project/endf/) using `pip`:

```sh
python -m pip install endf
```

## Usage

To read an ENDF file with a single material, simply pass the filename to the
`Material` class:
```python
import endf

mat = endf.Material('n-092_U_235.endf')
```

The `section_text` attribute holds a dictionary mapping (MF, MT) to the section
of text from the ENDF file:
```python
>>> print(mat.section_text[12, 75])
 9.223500+4 2.330248+2          2          2         25          0922812 75
 4.386000+5 0.000000+0          0          0          3          1922812 75
 0.000000+0 1.000000+0 1.000000+0                                 922812 75
```

The `section_data` attribute also holds a dictionary with keys that are (MF, MT)
pairs, but the values are dictionaries that hold the individual pieces of data
from the ENDF file:
```python
>>> mat.section_data[3, 16]
{'ZA': 92235,
 'AWR': 233.0248,
 'QM': -5298000.0,
 'QI': -5298000.0,
 'LR': 0,
 'sigma': <Tabulated1D: 39 points, 1 regions>}
 ```

While this form is more useful, it still may be a little too "raw". The
`Material` class has an `interpret()` method that returns a class based on the
ENDF sublibrary type (for example, incident-neutron data will result in an instance of the `IncidentNeutron` class):
```python
```
