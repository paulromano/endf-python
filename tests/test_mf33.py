# SPDX-FileCopyrightText: 2026 OpenMC contributors
# SPDX-License-Identifier: MIT

from io import StringIO

from pytest import approx
from endf.mf33 import parse_mf33


def _field(x):
    return f'{x:11d}' if isinstance(x, int) else f'{x:11.4e}'


def _cont(*items):
    return ''.join(_field(x) for x in items) + '\n'


def _list(c1, c2, l1, l2, n2, values):
    lines = _cont(c1, c2, l1, l2, len(values), n2)
    for i in range(0, len(values), 6):
        lines += ''.join(_field(x) for x in values[i:i + 6]) + '\n'
    return lines


def _mf33_section():
    text = _cont(2.6056e4, 55.3667, 0, 0, 0, 1)  # HEAD, NL=1
    text += _cont(0.0, 0.0, 0, 102, 2, 3)  # subsection, NC=2, NI=3

    # NC-type sub-subsection, LTY=0: sigma = 1*sigma(MT=4) + 1*sigma(MT=51)
    text += _cont(0.0, 0.0, 0, 0, 0, 0)
    text += _list(1.0e6, 2.0e7, 0, 0, 2, [1.0, 4.0, 1.0, 51.0])

    # NC-type sub-subsection, LTY=1
    text += _cont(0.0, 0.0, 0, 1, 0, 0)
    text += _list(1.0e3, 1.0e6, 9237, 18, 2, [33.0, 0.0, 1.0e3, 0.1, 1.0e6, 0.2])

    # NI-type sub-subsection, LB=1, LT=0: a single Ek-table with NP=3 pairs
    text += _list(0.0, 0.0, 0, 1, 3, [1.0e-5, 0.01, 1.0e6, 0.02, 2.0e7, 0.0])

    # NI-type sub-subsection, LB=3, LT=2: NP - LT = 3 pairs in the Ek-table,
    # LT = 2 pairs in the El-table
    text += _list(0.0, 0.0, 2, 3, 5, [1.0e-5, 0.1, 1.0e6, 0.2, 2.0e7, 0.0,
                                      1.0e-5, 0.3, 2.0e7, 0.0])

    # NI-type sub-subsection, LB=4, LT=1: NP - LT = 2 pairs in the Ek-table,
    # LT = 1 pair in the El-table
    text += _list(0.0, 0.0, 1, 4, 3, [1.0e-5, 0.4, 2.0e7, 0.0, 1.0e-5, 0.5])
    return StringIO(text)


def test_mf33_nc_subsections():
    data = parse_mf33(_mf33_section())
    subsection, = data['subsections']
    nc = subsection['nc_subsections']
    assert len(nc) == subsection['NC'] == 2

    assert nc[0]['LTY'] == 0
    assert nc[0]['E1'] == approx(1.0e6)
    assert nc[0]['E2'] == approx(2.0e7)
    assert nc[0]['NCI'] == 2
    assert nc[0]['CI'] == approx([1.0, 1.0])
    assert nc[0]['XMTI'] == approx([4.0, 51.0])

    assert nc[1]['LTY'] == 1
    assert nc[1]['MATS'] == 9237
    assert nc[1]['MTS'] == 18
    assert nc[1]['XMFS'] == approx(33.0)
    assert nc[1]['EI'] == approx([1.0e3, 1.0e6])
    assert nc[1]['WEI'] == approx([0.1, 0.2])


def test_mf33_ni_lb_0_4_tables():
    data = parse_mf33(_mf33_section())
    subsection, = data['subsections']
    ni = subsection['ni_subsections']
    assert len(ni) == subsection['NI'] == 3

    # LT=0: all NP pairs belong to the Ek-table
    lb1 = ni[0]
    assert (lb1['LB'], lb1['LT'], lb1['NP'], lb1['NT']) == (1, 0, 3, 6)
    assert lb1['Ek'] == approx([1.0e-5, 1.0e6, 2.0e7])
    assert lb1['Fk'] == approx([0.01, 0.02, 0.0])
    assert len(lb1['El']) == len(lb1['Fl']) == 0

    # LT>0: first NP - LT pairs in the Ek-table, last LT pairs in the El-table
    lb3 = ni[1]
    assert (lb3['LB'], lb3['LT'], lb3['NP'], lb3['NT']) == (3, 2, 5, 10)
    assert lb3['Ek'] == approx([1.0e-5, 1.0e6, 2.0e7])
    assert lb3['Fk'] == approx([0.1, 0.2, 0.0])
    assert lb3['El'] == approx([1.0e-5, 2.0e7])
    assert lb3['Fl'] == approx([0.3, 0.0])

    lb4 = ni[2]
    assert (lb4['LB'], lb4['LT'], lb4['NP'], lb4['NT']) == (4, 1, 3, 6)
    assert lb4['Ek'] == approx([1.0e-5, 2.0e7])
    assert lb4['Fk'] == approx([0.4, 0.0])
    assert lb4['El'] == approx([1.0e-5])
    assert lb4['Fl'] == approx([0.5])
