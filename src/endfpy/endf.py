"""Module for parsing and manipulating data from ENDF evaluations.

All the classes and functions in this module are based on document
ENDF-102 titled "Data Formats and Procedures for the Evaluated Nuclear
Data File ENDF-6". The latest version from June 2009 can be found at
http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

"""
import io
from pathlib import PurePath
from typing import TextIO
from warnings import warn

import numpy as np

from .data import gnds_name
from .energy_distribution import ArbitraryTabulated, GeneralEvaporation, \
    MaxwellEnergy, Evaporation,  WattEnergy, MadlandNix
from .records import get_head_record, get_text_record, get_cont_record, \
    get_tab1_record, get_list_record, get_tab2_record


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

SUM_RULES = {1: [2, 3],
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
             107: list(range(800, 850))}


def get_evaluations(filename):
    """Return a list of all evaluations within an ENDF file.

    Parameters
    ----------
    filename : str
        Path to ENDF-6 formatted file

    Returns
    -------
    list
        A list of :class:`openmc.data.endf.Evaluation` instances.

    """
    evaluations = []
    with open(str(filename), 'r') as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if line[66:70] == '  -1':
                break
            fh.seek(pos)
            evaluations.append(Evaluation(fh))
    return evaluations


def parse_mf3(file_obj: TextIO) -> dict:
    # Generate cross section
    ZA, AWR, *_ = get_head_record(file_obj)
    params, xs = get_tab1_record(file_obj)
    return {
        'ZA': ZA,
        'AWR': AWR,
        'QM': params[0],
        'QI': params[1],
        'LR': params[3],
        'sigma': xs
    }


def parse_mf4(file_obj: TextIO) -> dict:
    # Read first two records
    ZA, AWR, LVT, LTT, _, _ = get_head_record(file_obj)
    _, _, LI, LCT, NK, NM = get_cont_record(file_obj)

    # initialize dictionary for angular distribution
    data = {'ZA': ZA, 'AWR': AWR, 'LTT': LTT, 'LI': LI, 'LCT': LCT}

    # Check for obsolete energy transformation matrix. If present, just skip
    # it and keep reading
    if LVT > 0:
        warn('Obsolete energy transformation matrix in MF=4 angular distribution.')
        for _ in range((NK + 5)//6):
            file_obj.readline()

    def legendre_data(file_obj):
        data = {}
        params, data['E_int'] = get_tab2_record(file_obj)
        n_energy = params[5]

        energy = np.zeros(n_energy)
        a_l = []
        for i in range(n_energy):
            items, al = get_list_record(file_obj)
            data['T'] = items[0]
            energy[i] = items[1]
            data['LT'] = items[2]
            coefficients = np.array(al)
            a_l.append(coefficients)
        data['a_l'] = a_l
        data['E'] = energy
        return data

    def tabulated_data(file_obj):
        data = {}
        params, data['E_int'] = get_tab2_record(file_obj)
        n_energy = params[5]

        energy = np.zeros(n_energy)
        mu = []
        for i in range(n_energy):
            params, f = get_tab1_record(file_obj)
            data['T'] = params[0]
            energy[i] = params[1]
            data['LT'] = params[2]
            mu.append(f)
        data['E'] = energy
        data['mu'] = mu
        return data

    if LTT == 0 and LI == 1:
        # Purely isotropic
        pass

    elif LTT == 1 and LI == 0:
        # Legendre polynomial coefficients
        data['legendre'] = legendre_data(file_obj)

    elif LTT == 2 and LI == 0:
        # Tabulated probability distribution
        data['tabulated'] = tabulated_data(file_obj)

    elif LTT == 3 and LI == 0:
        # Legendre for low energies / tabulated for high energies
        data['legendre'] = legendre_data(file_obj)
        data['tabulated'] = tabulated_data(file_obj)

    return data


def parse_mf5(file_obj: TextIO) -> dict:
    ZA, AWR, _, _, NK, _ = get_head_record(file_obj)

    data = {'ZA': ZA, 'AWR': AWR, 'NK': NK}
    data['subsections'] = []
    for _ in range(NK):
        subsection = {}
        params, applicability = get_tab1_record(file_obj)
        subsection['LF'] = LF = params[3]
        subsection['p'] = applicability
        if LF == 1:
            dist = ArbitraryTabulated.dict_from_endf(file_obj, params)
        elif LF == 5:
            return GeneralEvaporation.dict_from_endf(file_obj, params)
        elif LF == 7:
            return MaxwellEnergy.dict_from_endf(file_obj, params)
        elif LF == 9:
            return Evaporation.dict_from_endf(file_obj, params)
        elif LF == 11:
            return WattEnergy.dict_from_endf(file_obj, params)
        elif LF == 12:
            return MadlandNix.dict_from_endf(file_obj, params)

        subsection['distribution'] = dist
        data['subsections'].append(subsection)

    return data


class Evaluation:
    """ENDF material evaluation with multiple files/sections

    Parameters
    ----------
    filename_or_obj : str or file-like
        Path to ENDF file to read or an open file positioned at the start of an
        ENDF material

    Attributes
    ----------
    info : dict
        Miscellaneous information about the evaluation.
    target : dict
        Information about the target material, such as its mass, isomeric state,
        whether it's stable, and whether it's fissionable.
    projectile : dict
        Information about the projectile such as its mass.
    section : dict
        Dictionary mapping (MF, MT) to corresponding section of the ENDF file.

    """
    def __init__(self, filename_or_obj):
        if isinstance(filename_or_obj, (str, PurePath)):
            fh = open(str(filename_or_obj), 'r')
            need_to_close = True
        else:
            fh = filename_or_obj
            need_to_close = False
        self.section = {}

        # Skip TPID record. Evaluators sometimes put in TPID records that are
        # ill-formated because they lack MF/MT values or put them in the wrong
        # columns.
        if fh.tell() == 0:
            fh.readline()
        MF = 0

        # Determine MAT number for this evaluation
        while MF == 0:
            position = fh.tell()
            line = fh.readline()
            MF = int(line[70:72])
        self.material = int(line[66:70])
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
            self.section[MF, MT] = section_text

        if need_to_close:
            fh.close()

        self.section_data = {}

        for (mf, mt), text in self.section.items():
            file_obj = io.StringIO(text)
            if mf == 1 and mt == 451:
                self.section_data[mf, mt] = self._read_mf1_mt451(file_obj)
            elif mf == 3:
                self.section_data[mf, mt] = parse_mf3(file_obj)
            elif mf == 4:
                self.section_data[mf, mt] = parse_mf4(file_obj)
            elif mf == 5:
                self.section_data[mf, mt] = parse_mf5(file_obj)

    def __repr__(self):
        metadata = self.section_data[1, 451]
        name = metadata['ZSYMAM'].replace(' ', '')
        return '<{} for {} {}>'.format(_SUBLIBRARY[metadata['NSUB']], name,
                                       _LIBRARY[metadata['NLIB']])

    def _read_mf1_mt451(self, file_obj: TextIO) -> dict:
        # Information about target/projectile
        ZA, AWR, LRP, LFI, NLIB, NMOD = get_head_record(file_obj)
        data = {
            'ZA': ZA, 'AWR': AWR, 'LRP': LRP,
            'LFI': LFI, 'NLIB': NLIB, 'NMOD': NMOD
        }

        # Control record 1
        ELIS, STA, LIS, LISO, _, NFOR = get_cont_record(file_obj)
        data['ELIS'] = ELIS
        data['STA'] = STA
        data['LIS'] = LIS
        data['LISO'] = LISO
        data['NFOR'] = NFOR

        # Control record 2
        AWI, EMAX, LREL, _, NSUB, NVER = get_cont_record(file_obj)
        data['AWI'] = AWI
        data['EMAX'] = EMAX
        data['LREL'] = LREL
        data['NSUB'] = NSUB
        data['NVER'] = NVER

        # Control record 3
        TEMP, _, LDRV, _, NWD, NXC = get_cont_record(file_obj)
        data['TEMP'] = TEMP
        data['LDRV'] = LDRV
        data['NWD'] = NWD
        data['NXC'] = NXC

        # Text records
        text = [get_text_record(file_obj) for i in range(NWD)]
        if len(text) >= 5:
            data['ZSYMAM'] = text[0][0:11]
            data['ALAB'] = text[0][11:22]
            data['EDATE'] = text[0][22:32]
            data['AUTH'] = text[0][32:66]
            data['REF'] = text[1][1:22]
            data['DDATE'] = text[1][22:32]
            data['RDATE'] = text[1][33:43]
            data['ENDATE'] = text[1][55:63]
            data['HSUB'] = text[2:5]
            data['description'] = text[5:]
        else:
            data['ZSYMAM'] = None

        # File numbers, reaction designations, and number of records
        data['section_list'] = []
        for _ in range(NXC):
            _, _, mf, mt, nc, mod = get_cont_record(file_obj, skip_c=True)
            data['section_list'].append((mf, mt, nc, mod))

        return data

    @property
    def gnds_name(self):
        return gnds_name(self.target['atomic_number'],
                         self.target['mass_number'],
                         self.target['isomeric_state'])
