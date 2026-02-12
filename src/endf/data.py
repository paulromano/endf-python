# SPDX-FileCopyrightText: 2023-2025 Paul Romano
# SPDX-License-Identifier: MIT

import re
from typing import Tuple

# Dictionary to give element symbols from IUPAC names
# (and some common mispellings)
ELEMENT_SYMBOL = {
    'neutron': 'n', 'hydrogen': 'H', 'helium': 'He',
    'lithium': 'Li', 'beryllium': 'Be', 'boron': 'B',
    'carbon': 'C', 'nitrogen': 'N', 'oxygen': 'O', 'fluorine': 'F',
    'neon': 'Ne', 'sodium': 'Na', 'magnesium': 'Mg',
    'aluminium': 'Al', 'aluminum': 'Al', 'silicon': 'Si',
    'phosphorus': 'P', 'sulfur': 'S', 'sulphur': 'S',
    'chlorine': 'Cl', 'argon': 'Ar', 'potassium': 'K',
    'calcium': 'Ca', 'scandium': 'Sc', 'titanium': 'Ti',
    'vanadium': 'V', 'chromium': 'Cr', 'manganese': 'Mn',
    'iron': 'Fe', 'cobalt': 'Co', 'nickel': 'Ni', 'copper': 'Cu',
    'zinc': 'Zn', 'gallium': 'Ga', 'germanium': 'Ge',
    'arsenic': 'As', 'selenium': 'Se', 'bromine': 'Br',
    'krypton': 'Kr', 'rubidium': 'Rb', 'strontium': 'Sr',
    'yttrium': 'Y', 'zirconium': 'Zr', 'niobium': 'Nb',
    'molybdenum': 'Mo', 'technetium': 'Tc', 'ruthenium': 'Ru',
    'rhodium': 'Rh', 'palladium': 'Pd', 'silver': 'Ag',
    'cadmium': 'Cd', 'indium': 'In', 'tin': 'Sn', 'antimony': 'Sb',
    'tellurium': 'Te', 'iodine': 'I', 'xenon': 'Xe',
    'caesium': 'Cs', 'cesium': 'Cs', 'barium': 'Ba',
    'lanthanum': 'La', 'cerium': 'Ce', 'praseodymium': 'Pr',
    'neodymium': 'Nd', 'promethium': 'Pm', 'samarium': 'Sm',
    'europium': 'Eu', 'gadolinium': 'Gd', 'terbium': 'Tb',
    'dysprosium': 'Dy', 'holmium': 'Ho', 'erbium': 'Er',
    'thulium': 'Tm', 'ytterbium': 'Yb', 'lutetium': 'Lu',
    'hafnium': 'Hf', 'tantalum': 'Ta', 'tungsten': 'W',
    'wolfram': 'W', 'rhenium': 'Re', 'osmium': 'Os',
    'iridium': 'Ir', 'platinum': 'Pt', 'gold': 'Au',
    'mercury': 'Hg', 'thallium': 'Tl', 'lead': 'Pb',
    'bismuth': 'Bi', 'polonium': 'Po', 'astatine': 'At',
    'radon': 'Rn', 'francium': 'Fr', 'radium': 'Ra',
    'actinium': 'Ac', 'thorium': 'Th', 'protactinium': 'Pa',
    'uranium': 'U', 'neptunium': 'Np', 'plutonium': 'Pu',
    'americium': 'Am', 'curium': 'Cm', 'berkelium': 'Bk',
    'californium': 'Cf', 'einsteinium': 'Es', 'fermium': 'Fm',
    'mendelevium': 'Md', 'nobelium': 'No', 'lawrencium': 'Lr',
    'rutherfordium': 'Rf', 'dubnium': 'Db', 'seaborgium': 'Sg',
    'bohrium': 'Bh', 'hassium': 'Hs', 'meitnerium': 'Mt',
    'darmstadtium': 'Ds', 'roentgenium': 'Rg', 'copernicium': 'Cn',
    'nihonium': 'Nh', 'flerovium': 'Fl', 'moscovium': 'Mc',
    'livermorium': 'Lv', 'tennessine': 'Ts', 'oganesson': 'Og'
}

ATOMIC_SYMBOL = {
    0: 'n', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',
    7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al',
    14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K',
    20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn',
    26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga',
    32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb',
    38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
    44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
    50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs',
    56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm',
    62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho',
    68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta',
    74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
    80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
    86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa',
    92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk',
    98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No',
    103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh',
    108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Cn',
    113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts',
    118: 'Og'
}
ATOMIC_NUMBER = {value: key for key, value in ATOMIC_SYMBOL.items()}

# Boltzmann constant in [eV/K]
K_BOLTZMANN = 8.617333262e-5

EV_PER_MEV = 1.0e6

# Regex for GNDS nuclide names (used in zam function)
_GNDS_NAME_RE = re.compile(r'([A-Zn][a-z]*)(\d+)((?:_[em]\d+)?)')

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
    107: list(range(800, 850)),
    501: [502, 504, 516, 522],
    516: [515, 517],
    522: list(range(534, 573)),    
}

def gnds_name(Z: int, A: int, m: int = 0) -> str:
    """Return nuclide name using GNDS convention

    Parameters
    ----------
    Z
        Atomic number
    A
        Mass number
    m
        Metastable state

    Returns
    -------
    Nuclide name in GNDS convention, e.g., 'Am242_m1'

    """
    if m > 0:
        return f'{ATOMIC_SYMBOL[Z]}{A}_m{m}'
    return f'{ATOMIC_SYMBOL[Z]}{A}'


def zam(name: str) -> Tuple[int, int, int]:
    """Return tuple of (atomic number, mass number, metastable state)

    Parameters
    ----------
    name
        Name of nuclide using GNDS convention, e.g., 'Am242_m1'

    Returns
    -------
    Atomic number, mass number, and metastable state

    """
    try:
        symbol, A, state = _GNDS_NAME_RE.fullmatch(name).groups()
    except AttributeError:
        raise ValueError(f"'{name}' does not appear to be a nuclide name in "
                         "GNDS format")

    if symbol not in ATOMIC_NUMBER:
        raise ValueError(f"'{symbol}' is not a recognized element symbol")

    metastable = int(state[2:]) if state else 0
    return (ATOMIC_NUMBER[symbol], int(A), metastable)


def temperature_str(T: float) -> str:
    """Return temperature as a string

    Parameters
    ----------
    T
        Temperature in [K]

    Returns
    -------
    String representation of temperature, e.g., '294K'

    """
    return "{}K".format(int(round(T)))
