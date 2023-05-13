# SPDX-FileCopyrightText: 2023 International Atomic Energy Agency
# SPDX-License-Identifier: MIT

"""Module for parsing and manipulating data from ENDF files.

All the classes and functions in this module are based on the ENDF-102 report
titled "ENDF-6 Formats Manual: Data Formats and Procedures for the Evaluated
Nuclear Data Files". The version from January 2018 can be found at
https://doi.org/10.2172/1425114.

"""
import io
from typing import List, Tuple, Any, Union, TextIO

import endf
from .data import gnds_name
from .fileutils import PathLike
from .mf1 import parse_mf1_mt451, parse_mf1_mt452, parse_mf1_mt455, \
    parse_mf1_mt458, parse_mf1_mt460
from .mf2 import parse_mf2
from .mf3 import parse_mf3
from .mf4 import parse_mf4
from .mf5 import parse_mf5
from .mf6 import parse_mf6
from .mf7 import parse_mf7_mt2, parse_mf7_mt4
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

SUM_RULES = {
    1: [2, 3],
    3: [4, 5, 11, 16, 17, 22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 35,
        36, 37, 41, 42, 44, 45, 152, 153, 154, 156, 157, 158, 159, 160,
        161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
        173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185,
        186, 187, 188, 189, 190, 194, 195, 196, 198, 199, 200],
    4: list(range(50, 92)),
    16: list(range(875, 892)),
    18: [19, 20, 21, 38],
    27: [18, 101],
    101: [102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114,
        115, 116, 117, 155, 182, 191, 192, 193, 197],
    103: list(range(600, 650)),
    104: list(range(650, 700)),
    105: list(range(700, 750)),
    106: list(range(750, 800)),
    107: list(range(800, 850))
}


class Material:
    """ENDF material with multiple files/sections

    Parameters
    ----------
    filename
        Path to ENDF file to read or an open file positioned at the start of an
        ENDF material

    Attributes
    ----------
    MAT : int
        ENDF material number
    sections : list of tuple
        List of (MF, MT) sections
    section_text : dict
        Dictionary mapping (MF, MT) to corresponding section of the ENDF file.
    section_data : dict
        Dictionary mapping (MF, MT) to a dictionary representing the
        corresponding section of the ENDF file.

    """
    def __init__(self, filename_or_obj: Union[PathLike, TextIO]):
        if isinstance(filename_or_obj, PathLike.__args__):
            fh = open(str(filename_or_obj), 'r')
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
            else:
                pass

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
    def gnds_name(self) -> str:
        metadata = self[1, 451]
        Z, A = divmod(metadata['ZA'], 1000)
        m = metadata['LISO']
        return gnds_name(Z, A, m)

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
            return endf.IncidentNeutron(self)
        else:
            raise NotImplementedError(f"No class implemented for {NSUB=}")


def get_materials(filename: PathLike) -> List[Material]:
    """Return a list of all materials within an ENDF file.

    Parameters
    ----------
    filename
        Path to ENDF-6 formatted file

    Returns
    -------
    A list of :class:`Material` instances.

    """
    materials = []
    with open(str(filename), 'r') as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if line[66:70] == '  -1':
                break
            fh.seek(pos)
            materials.append(Material(fh))
    return materials
