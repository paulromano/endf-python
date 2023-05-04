from typing import TextIO

from .records import get_head_record, get_cont_record, get_text_record


def parse_mf1_mt451(file_obj: TextIO) -> dict:
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
