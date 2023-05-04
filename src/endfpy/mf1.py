from typing import TextIO

from .records import get_head_record, get_cont_record, get_text_record, \
    get_list_record, get_tab1_record, get_tab2_record


def parse_mf1_mt451(file_obj: TextIO) -> dict:
    """Parse descriptive data and directory from MF=1, MT=451

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Descriptive data

    """
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


def parse_mf1_mt452(file_obj: TextIO) -> dict:
    """Parse number of neutrons per fission from MF=1, MT=452/456

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Data on number of neutrons per fission

    """
    # Determine representation from HEAD record
    ZA, AWR, _, LNU, _, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LNU': LNU}

    if LNU == 1:
        # Polynomial representation
        _, data['C'] = get_list_record(file_obj)
    elif LNU == 2:
        # Tabulated representation
        _, data['nu'] = get_tab1_record(file_obj)

    return data


def parse_mf1_mt455(file_obj: TextIO) -> dict:
    """Parse delayed neutron data from MF=1, MT=455

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Delayed neutron data

    """
    ZA, AWR, LDG, LNU, _, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LDG': LDG, 'LNU': LNU}

    if LDG == 0:
        # Delayed-group constants energy independent
        _, data['lambda'] = get_list_record(file_obj)
    elif LDG == 1:
        # Delayed-group constants energy dependent
        params, data['E_int'] = get_tab2_record(file_obj)
        NE = params[5]

        data['constants'] = []
        for _ in range(NE):
            (_, E, *_), values = get_list_record(file_obj)
            data['constants'].append({
                'E': E, 'lambda': values[::2], 'alpha': values[1::2]
            })

    # In MF=1, MT=455, the delayed-group abundances are actually not
    # specified if the group constants are energy-independent. In this case,
    # the abundances must be inferred from MF=5, MT=455 where multiple
    # energy distributions are given.
    if LNU == 1:
        # Nu represented as polynomial
        _, data['C'] = get_list_record(file_obj)
    elif LNU == 2:
        # Nu represented by tabulation
        _, data['nu'] = get_tab1_record(file_obj)

    return data
