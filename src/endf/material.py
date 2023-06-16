# SPDX-FileCopyrightText: 2023 OpenMC contributors and Paul Romano
# SPDX-License-Identifier: MIT

"""Module for parsing and manipulating data from ENDF files.

All the classes and functions in this module are based on the ENDF-102 report
titled "ENDF-6 Formats Manual: Data Formats and Procedures for the Evaluated
Nuclear Data Files". The version from January 2018 can be found at
https://doi.org/10.2172/1425114.

"""
import io
from typing import List, Tuple, Any, Union, TextIO, Optional
from warnings import warn

import endf
from .fileutils import PathLike
from .mf1 import parse_mf1_mt451, parse_mf1_mt452, parse_mf1_mt455, \
    parse_mf1_mt458, parse_mf1_mt460
from .mf2 import parse_mf2
from .mf3 import parse_mf3
from .mf4 import parse_mf4
from .mf5 import parse_mf5
from .mf6 import parse_mf6
from .mf7 import parse_mf7_mt2, parse_mf7_mt4, parse_mf7_mt451
from .mf8 import parse_mf8, parse_mf8_mt454, parse_mf8_mt457
from .mf9 import parse_mf9_mf10
from .mf12 import parse_mf12
from .mf13 import parse_mf13
from .mf14 import parse_mf14
from .mf15 import parse_mf15
from .mf23 import parse_mf23
from .mf26 import parse_mf26
from .mf27 import parse_mf27
from .mf28 import parse_mf28
from .mf33 import parse_mf33
from .mf34 import parse_mf34
from .mf40 import parse_mf40


_LIBRARY = {
    0: 'ENDF/B',
    1: 'ENDF/A',
    2: 'JEFF',
    3: 'EFF',
    4: 'ENDF/B High Energy',
    5: 'CENDL',
    6: 'JENDL',
    17: 'TENDL',
    18: 'ROSFOND',
    21: 'SG-23',
    31: 'INDL/V',
    32: 'INDL/A',
    33: 'FENDL',
    34: 'IRDF',
    35: 'BROND',
    36: 'INGDB-90',
    37: 'FENDL/A',
    38: 'IAEA/PD',
    41: 'BROND'
}

_SUBLIBRARY = {
    0: 'Photo-nuclear data',
    1: 'Photo-induced fission product yields',
    3: 'Photo-atomic data',
    4: 'Radioactive decay data',
    5: 'Spontaneous fission product yields',
    6: 'Atomic relaxation data',
    10: 'Incident-neutron data',
    11: 'Neutron-induced fission product yields',
    12: 'Thermal neutron scattering data',
    19: 'Neutron standards',
    113: 'Electro-atomic data',
    10010: 'Incident-proton data',
    10011: 'Proton-induced fission product yields',
    10020: 'Incident-deuteron data',
    10030: 'Incident-triton data',
    20030: 'Incident-helion (3He) data',
    20040: 'Incident-alpha data'
}



class Material:
    """ENDF material with multiple files/sections

    Parameters
    ----------
    filename_or_obj
        Path to ENDF file to read or an open file positioned at the start of an
        ENDF material
    encoding
        Encoding of the ENDF-6 formatted file

    Attributes
    ----------
    MAT
        ENDF material number
    sections
        List of (MF, MT) sections
    section_text
        Dictionary mapping (MF, MT) to corresponding section of the ENDF file.
    section_data
        Dictionary mapping (MF, MT) to a dictionary representing the
        corresponding section of the ENDF file.

    """
    # TODO: Remove need to list properties here
    MAT: int
    sections: List[Tuple[int, int]]
    section_text: dict
    section_data: dict

    def __init__(self, filename_or_obj: Union[PathLike, TextIO], encoding: Optional[str] = None):
        if isinstance(filename_or_obj, PathLike.__args__):
            fh = open(str(filename_or_obj), 'r', encoding=encoding)
            need_to_close = True
        else:
            fh = filename_or_obj
            need_to_close = False
        self.section_text = {}

        # Skip TPID record. Evaluators sometimes put in TPID records that are
        # ill-formated because they lack MF/MT values or put them in the wrong
        # columns.
        if fh.tell() == 0:
            fh.readline()
        MF = 0

        # Determine MAT number for this material
        while MF == 0:
            position = fh.tell()
            line = fh.readline()
            MF = int(line[70:72])
        self.MAT = int(line[66:70])
        fh.seek(position)

        while True:
            # Find next section
            while True:
                position = fh.tell()
                line = fh.readline()
                MAT = int(line[66:70])
                MF = int(line[70:72])
                MT = int(line[72:75])
                if MT > 0 or MAT == 0:
                    fh.seek(position)
                    break

            # If end of material reached, exit loop
            if MAT == 0:
                fh.readline()
                break

            section_text = ''
            while True:
                line = fh.readline()
                if line[72:75] == '  0':
                    break
                else:
                    section_text += line
            self.section_text[MF, MT] = section_text

        if need_to_close:
            fh.close()

        self.section_data = {}
        for (MF, MT), text in self.section_text.items():
            file_obj = io.StringIO(text)
            if MF == 1 and MT == 451:
                self.section_data[MF, MT] = parse_mf1_mt451(file_obj)
            elif MF == 1 and MT in (452, 456):
                self.section_data[MF, MT] = parse_mf1_mt452(file_obj)
            elif MF == 1 and MT == 455:
                self.section_data[MF, MT] = parse_mf1_mt455(file_obj)
            elif MF == 1 and MT == 458:
                self.section_data[MF, MT] = parse_mf1_mt458(file_obj)
            elif MF == 1 and MT == 460:
                self.section_data[MF, MT] = parse_mf1_mt460(file_obj)
            elif MF == 2 and MT == 151:
                self.section_data[MF, MT] = parse_mf2(file_obj)
            elif MF == 3:
                self.section_data[MF, MT] = parse_mf3(file_obj)
            elif MF == 4:
                self.section_data[MF, MT] = parse_mf4(file_obj)
            elif MF == 5:
                self.section_data[MF, MT] = parse_mf5(file_obj)
            elif MF == 6:
                self.section_data[MF, MT] = parse_mf6(file_obj)
            elif MF == 7 and MT == 2:
                self.section_data[MF, MT] = parse_mf7_mt2(file_obj)
            elif MF == 7 and MT == 4:
                self.section_data[MF, MT] = parse_mf7_mt4(file_obj)
            elif MF == 7 and MT == 451:
                self.section_data[MF, MT] = parse_mf7_mt451(file_obj)
            elif MF == 8 and MT in (454, 459):
                self.section_data[MF, MT] = parse_mf8_mt454(file_obj)
            elif MF == 8 and MT == 457:
                self.section_data[MF, MT] = parse_mf8_mt457(file_obj)
            elif MF == 8:
                self.section_data[MF, MT] = parse_mf8(file_obj)
            elif MF in (9, 10):
                self.section_data[MF, MT] = parse_mf9_mf10(file_obj, MF)
            elif MF == 12:
                self.section_data[MF, MT] = parse_mf12(file_obj)
            elif MF == 13:
                self.section_data[MF, MT] = parse_mf13(file_obj)
            elif MF == 14:
                self.section_data[MF, MT] = parse_mf14(file_obj)
            elif MF == 15:
                self.section_data[MF, MT] = parse_mf15(file_obj)
            elif MF == 23:
                self.section_data[MF, MT] = parse_mf23(file_obj)
            elif MF == 26:
                self.section_data[MF, MT] = parse_mf26(file_obj)
            elif MF == 27:
                self.section_data[MF, MT] = parse_mf27(file_obj)
            elif MF == 28:
                self.section_data[MF, MT] = parse_mf28(file_obj)
            elif MF == 33:
                self.section_data[MF, MT] = parse_mf33(file_obj)
            elif MF == 34:
                self.section_data[MF, MT] = parse_mf34(file_obj, MT)
            elif MF == 40:
                self.section_data[MF, MT] = parse_mf40(file_obj)
            else:
                warn(f"{MF=}, {MT=} ignored")

    def __contains__(self, mf_mt: Tuple[int, int]) -> bool:
        return mf_mt in self.section_data

    def __getitem__(self, mf_mt: Tuple[int, int]) -> dict:
        return self.section_data[mf_mt]

    def __setitem__(self, key: Tuple[int, int], value):
        self.section_data[key] = value

    def __repr__(self) -> str:
        metadata = self.section_data[1, 451]
        name = metadata['ZSYMAM'].replace(' ', '')
        return '<{} for {} {}>'.format(_SUBLIBRARY[metadata['NSUB']], name,
                                       _LIBRARY[metadata['NLIB']])

    @property
    def sections(self) -> List[Tuple[int, int]]:
        return list(self.section_text.keys())

    def interpret(self) -> Any:
        """Get high-level interface class for the ENDF material

        Returns
        -------
        Instance of a high-level interface class, e.g.,
        :class:`endf.IncidentNeutron`.

        """
        NSUB = self.section_data[1, 451]['NSUB']
        if NSUB == 10:
            return endf.IncidentNeutron.from_endf(self)
        else:
            raise NotImplementedError(f"No class implemented for {NSUB=}")


def get_materials(filename: PathLike, encoding: Optional[str] = None) -> List[Material]:
    """Return a list of all materials within an ENDF file.

    Parameters
    ----------
    filename
        Path to ENDF-6 formatted file
    encoding
        Encoding of the ENDF-6 formatted file

    Returns
    -------
    A list of ENDF materials

    """
    materials = []
    with open(str(filename), 'r', encoding=encoding) as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if line[66:70] == '  -1':
                break
            fh.seek(pos)
            materials.append(Material(fh))
    return materials
