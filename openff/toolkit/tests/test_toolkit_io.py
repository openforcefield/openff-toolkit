#!/usr/bin/env python
# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Tests for I/O functionality of the toolkit wrappers

"""

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import os
import sys
from io import BytesIO, StringIO
import pytest
import pathlib
import tempfile
import numpy as np
from numpy.testing import assert_allclose
from simtk import unit

from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper
from openff.toolkit.utils import exceptions
from openff.toolkit.topology.molecule import Molecule

from openff.toolkit.tests import create_molecules


# ================================================================
# Data records used for testing.
# ================================================================

ETHANOL = create_molecules.create_ethanol()
ETHANOL.name = "ethanol"

# ========================================================
# Various records for caffeine
# ========================================================

# From https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL113/
CAFFEINE_2D_SDF = """\
caffeine
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
   -1.1875   -9.6542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1875   -8.9625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8125  -10.0292    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4167   -8.9625    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4167   -9.6542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8125   -8.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -9.8917    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -8.7625    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1125   -9.3042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0250  -10.0375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8125   -7.8917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8125  -10.7417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0250   -8.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2917   -8.0750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  2  0
  3  1  1  0
  4  5  1  0
  5  3  1  0
  6  2  1  0
  7  1  1  0
  8  2  1  0
  9  7  2  0
 10  5  2  0
 11  6  2  0
 12  3  1  0
 13  4  1  0
 14  8  1  0
  9  8  1  0
  4  6  1  0
M  END

> <chembl_id>
CHEMBL113

> <chembl_pref_name>
CAFFEINE

$$$$
"""

CAFFEINE_2D_COORDS = np.array([
    
    (-1.1875, -9.6542, 0.0000),
    (-1.1875, -8.9625, 0.0000),
    (-1.8125, -10.0292, 0.0000),
    (-2.4167, -8.9625, 0.0000),
    (-2.4167, -9.6542, 0.0000),
    (-1.8125, -8.6000, 0.0000),
    (-0.5000, -9.8917, 0.0000),
    (-0.5000, -8.7625, 0.0000),
    (-0.1125, -9.3042, 0.0000),
    (-3.0250, -10.0375, 0.0000),
    (-1.8125, -7.8917, 0.0000),
    (-1.8125, -10.7417, 0.0000),
    (-3.0250, -8.6000, 0.0000),
    (-0.2917, -8.0750, 0.0000),
    ], np.double) * unit.angstrom


# From https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL113/
CAFFEINE_SMI = "Cn1c(=O)c2c(ncn2C)n(C)c1=O CHEMBL113\n"

# From https://pubchem.ncbi.nlm.nih.gov/compound/2519#section=2D-Structure
CAFFEINE_3D_SDF = """\
2519
  -OEChem-07012107543D

 24 25  0     0  0  0  0  0  0999 V2000
    0.4700    2.5688    0.0006 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1271   -0.4436   -0.0003 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9686   -1.3125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2182    0.1412   -0.0003 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3477    1.0797   -0.0001 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4119   -1.9372    0.0002 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8579    0.2592   -0.0008 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3897   -1.0264   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0307    1.4220   -0.0006 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9061   -0.2495   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5032   -1.1998    0.0003 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4276   -2.6960    0.0008 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1926    1.2061    0.0003 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2969    2.1881    0.0007 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5163   -1.5787    0.0008 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0451   -3.1973   -0.8937 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5186   -2.7596    0.0011 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0447   -3.1963    0.8957 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.1992    0.7801    0.0002 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0468    1.8092   -0.8992 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0466    1.8083    0.9004 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8087    3.1651   -0.0003 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9322    2.1027    0.8881 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9346    2.1021   -0.8849 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  9  2  0  0  0  0
  2 10  2  0  0  0  0
  3  8  1  0  0  0  0
  3 10  1  0  0  0  0
  3 12  1  0  0  0  0
  4  7  1  0  0  0  0
  4 11  1  0  0  0  0
  4 13  1  0  0  0  0
  5  9  1  0  0  0  0
  5 10  1  0  0  0  0
  5 14  1  0  0  0  0
  6  8  1  0  0  0  0
  6 11  2  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
 11 15  1  0  0  0  0
 12 16  1  0  0  0  0
 12 17  1  0  0  0  0
 12 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
 13 21  1  0  0  0  0
 14 22  1  0  0  0  0
 14 23  1  0  0  0  0
 14 24  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
2519

> <PUBCHEM_CONFORMER_RMSD>
0.4

> <PUBCHEM_CONFORMER_DIVERSEORDER>
1

> <PUBCHEM_MMFF94_PARTIAL_CHARGES>
15
1 -0.57
10 0.69
11 0.04
12 0.3
13 0.26
14 0.3
15 0.15
2 -0.57
3 -0.42
4 0.05
5 -0.42
6 -0.57
7 -0.24
8 0.29
9 0.71

> <PUBCHEM_EFFECTIVE_ROTOR_COUNT>
0

> <PUBCHEM_PHARMACOPHORE_FEATURES>
5
1 1 acceptor
1 2 acceptor
3 4 6 11 cation
5 4 6 7 8 11 rings
6 3 5 7 8 9 10 rings

> <PUBCHEM_HEAVY_ATOM_COUNT>
14

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
0

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
1

> <PUBCHEM_CACTVS_TAUTO_COUNT>
1

> <PUBCHEM_CONFORMER_ID>
000009D700000001

> <PUBCHEM_MMFF94_ENERGY>
22.901

> <PUBCHEM_FEATURE_SELFOVERLAP>
25.487

> <PUBCHEM_SHAPE_FINGERPRINT>
10967382 1 18338799025773621285
11132069 177 18339075025094499008
12524768 44 18342463625094026902
13140716 1 17978511158789908153
16945 1 18338517550775811621
193761 8 15816500986559935910
20588541 1 18339082691204868851
21501502 16 18338796715286957384
22802520 49 18128840606503503494
2334 1 18338516344016692929
23402539 116 18270382932679789735
23552423 10 18262240993325675966
23559900 14 18199193898169584358
241688 4 18266458702623303353
2748010 2 18266180539182415717
5084963 1 17698433339235542986
528886 8 18267580380709240570
53812653 166 18198902694142226312
66348 1 18339079396917369615

> <PUBCHEM_SHAPE_MULTIPOLES>
256.45
4.01
2.83
0.58
0.71
0.08
0
-0.48
0
-0.81
0
0.01
0
0

> <PUBCHEM_SHAPE_SELFOVERLAP>
550.88

> <PUBCHEM_SHAPE_VOLUME>
143.9

> <PUBCHEM_COORDINATE_TYPE>
2
5
10

$$$$
"""

CAFFEINE_3D_COORDS = np.array([
    (0.4700,    2.5688,    0.0006),
    (-3.1271,   -0.4436,   -0.0003),
    (-0.9686,   -1.3125,    0.0000),
    (2.2182,    0.1412,   -0.0003),
    (-1.3477,    1.0797,   -0.0001),
    (1.4119,   -1.9372,    0.0002),
    (0.8579,    0.2592,   -0.0008),
    (0.3897,   -1.0264,   -0.0004),
    (0.0307,    1.4220,   -0.0006),
    (-1.9061,   -0.2495,   -0.0004),
    (2.5032,   -1.1998,    0.0003),
    (-1.4276,   -2.6960,    0.0008),
    (3.1926,    1.2061,    0.0003),
    (-2.2969,    2.1881,    0.0007),
    (3.5163,   -1.5787,    0.0008),
    (-1.0451,   -3.1973,   -0.8937),
    (-2.5186,   -2.7596,    0.0011),
    (-1.0447,   -3.1963,    0.8957),
    (4.1992,    0.7801,    0.0002),
    (3.0468,    1.8092,   -0.8992),
    (3.0466,    1.8083,    0.9004),
    (-1.8087,    3.1651,   -0.0003),
    (-2.9322,    2.1027,    0.8881),
    (-2.9346,    2.1021,   -0.8849),
    ], np.double) * unit.angstrom

# ========================================================
# Various records for aspirin
# ========================================================


# From https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25/
ASPIRIN_2D_SDF = """\
aspirin

     RDKit          2D

 13 13  0  0  0  0  0  0  0  0999 V2000
    8.8810   -2.1206    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8798   -2.9479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5946   -3.3607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3110   -2.9474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3081   -2.1170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5928   -1.7078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0210   -1.7018    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.7369   -2.1116    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.0260   -3.3588    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.0273   -4.1837    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.7423   -4.5949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3136   -4.5972    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.0178   -0.8769    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  5  7  1  0
  3  4  2  0
  7  8  2  0
  4  9  1  0
  4  5  1  0
  9 10  1  0
  2  3  1  0
 10 11  1  0
  5  6  2  0
 10 12  2  0
  6  1  1  0
  7 13  1  0
M  END

> <chembl_id>
CHEMBL25

> <chembl_pref_name>
ASPIRIN

$$$$
"""

# From https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25/
ASPIRIN_SMI = "CC(=O)Oc1ccccc1C(=O)O ASPIRIN\n"

# From https://pubchem.ncbi.nlm.nih.gov/compound/2244
ASPIRIN_3D_SDF = """\
2244
  -OEChem-07012107583D

 21 21  0     0  0  0  0  0  0999 V2000
    1.2333    0.5540    0.7792 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6952   -2.7148   -0.7502 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7958   -2.1843    0.8685 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7813    0.8105   -1.4821 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0857    0.6088    0.4403 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7927   -0.5515    0.1244 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7288    1.8464    0.4133 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1426   -0.4741   -0.2184 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0787    1.9238    0.0706 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7855    0.7636   -0.2453 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1409   -1.8536    0.1477 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1094    0.6715   -0.3113 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5305    0.5996    0.1635 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1851    2.7545    0.6593 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7247   -1.3605   -0.4564 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5797    2.8872    0.0506 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8374    0.8238   -0.5090 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7290    1.4184    0.8593 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2045    0.6969   -0.6924 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7105   -0.3659    0.6426 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2555   -3.5916   -0.7337 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1 12  1  0  0  0  0
  2 11  1  0  0  0  0
  2 21  1  0  0  0  0
  3 11  2  0  0  0  0
  4 12  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  2  0  0  0  0
  6  8  2  0  0  0  0
  6 11  1  0  0  0  0
  7  9  1  0  0  0  0
  7 14  1  0  0  0  0
  8 10  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  2  0  0  0  0
  9 16  1  0  0  0  0
 10 17  1  0  0  0  0
 12 13  1  0  0  0  0
 13 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
2244

> <PUBCHEM_CONFORMER_RMSD>
0.6

> <PUBCHEM_CONFORMER_DIVERSEORDER>
1
11
10
3
15
17
13
5
16
7
14
9
8
4
18
6
12
2

> <PUBCHEM_MMFF94_PARTIAL_CHARGES>
18
1 -0.23
10 -0.15
11 0.63
12 0.66
13 0.06
14 0.15
15 0.15
16 0.15
17 0.15
2 -0.65
21 0.5
3 -0.57
4 -0.57
5 0.08
6 0.09
7 -0.15
8 -0.15
9 -0.15

> <PUBCHEM_EFFECTIVE_ROTOR_COUNT>
3

> <PUBCHEM_PHARMACOPHORE_FEATURES>
5
1 2 acceptor
1 3 acceptor
1 4 acceptor
3 2 3 11 anion
6 5 6 7 8 9 10 rings

> <PUBCHEM_HEAVY_ATOM_COUNT>
13

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
0

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
1

> <PUBCHEM_CACTVS_TAUTO_COUNT>
1

> <PUBCHEM_CONFORMER_ID>
000008C400000001

> <PUBCHEM_MMFF94_ENERGY>
39.5952

> <PUBCHEM_FEATURE_SELFOVERLAP>
25.432

> <PUBCHEM_SHAPE_FINGERPRINT>
1 1 18265615372930943622
100427 49 16967750034970055351
12138202 97 18271247217817981012
12423570 1 16692715976000295083
12524768 44 16753525617747228747
12716758 59 18341332292274886536
13024252 1 17968377969333732145
14181834 199 17830728755827362645
14614273 12 18262232214645093005
15207287 21 17703787037639964108
15775835 57 18340488876329928641
16945 1 18271533103414939405
193761 8 17907860604865584321
20645476 183 17677348215414174190
20871998 184 18198632231250704846
21040471 1 18411412921197846465
21501502 16 18123463883164380929
23402539 116 18271795865171824860
23419403 2 13539898140662769886
23552423 10 18048876295495619569
23559900 14 18272369794190581304
241688 4 16179044415907240795
257057 1 17478316999871287486
2748010 2 18339085878070479087
305870 269 18263645056784260212
528862 383 18117272558388284091
53812653 8 18410289211719108569
7364860 26 17910392788380644719
81228 2 18050568744116491203

> <PUBCHEM_SHAPE_MULTIPOLES>
244.06
3.86
2.45
0.89
1.95
1.58
0.15
-1.85
0.38
-0.61
-0.02
0.29
0.01
-0.33

> <PUBCHEM_SHAPE_SELFOVERLAP>
513.037

> <PUBCHEM_SHAPE_VOLUME>
136

> <PUBCHEM_COORDINATE_TYPE>
2
5
10

$$$$
"""


# ========================================================
# CHEBI:1148 triggers 'allow_undefined_stereo' exceptions
# ========================================================

CHEBI_1148_SDF = """\
CHEBI:1148
  Marvin  09120817212D          

  7  6  0  0  0  0            999 V2000
    1.4289   -0.1650    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.2475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.1650    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -0.7145    0.2475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289   -0.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.0725    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.9900    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  6  2  2  0  0  0  0
  3  2  1  0  0  0  0
  7  3  1  0  0  0  0
  4  3  1  0  0  0  0
  5  4  1  0  0  0  0
M  END
> <ChEBI ID>
CHEBI:1148

> <ChEBI Name>
2-hydroxybutyric acid

> <Star>
3

$$$$
"""

# ========================================================
# Used to test that _cls is passed correctly
# ========================================================

class SingingMolecule(Molecule):
    def sing(self):
        return "The hills are alive with sound of music!"

# ========================================================
# Manage the input files
# ========================================================

class FilenameDescriptor:
    def __init__(self, content):
        self.content = content

    def __set_name__(self, owner, name):
        self.private_name = f"_{name}_filename"
        basename, _, ext = name.rpartition("_")
        self.filename = f"{basename}.{ext}"

    def __get__(self, obj, objtype = None):
        # Have we been called before?
        filename = getattr(obj, self.private_name, None)
        if filename is None:
            # No.
            # Write the content to the named file in the manager's directory
            file_path = obj.dir_path / self.filename
            file_path.write_text(self.content)
            # Save the filename to obj's private_name.
            filename = str(file_path)
            setattr(obj, self.private_name, filename)
            
        # Return the path as a string
        return filename


class FileManager:
    def __init__(self, temporary_directory):
        # Use Python's object lifetime to manage directory cleanup
        self.temporary_directory = temporary_directory
        
        self.dir_path = pathlib.Path(temporary_directory.name)

    caffeine_2d_sdf = FilenameDescriptor(CAFFEINE_2D_SDF)
    caffeine_3d_sdf = FilenameDescriptor(CAFFEINE_3D_SDF)
    caffeine_smi = FilenameDescriptor(CAFFEINE_SMI)
    
    aspirin_2d_sdf = FilenameDescriptor(ASPIRIN_2D_SDF)
    aspirin_3d_sdf = FilenameDescriptor(ASPIRIN_3D_SDF)
    aspirin_smi = FilenameDescriptor(ASPIRIN_SMI)

    chebi_1148_sdf = FilenameDescriptor(CHEBI_1148_SDF)

    
file_manager = FileManager(tempfile.TemporaryDirectory())

class FileobjDescriptor:
    def __init__(self, content):
        self.content = content.encode("utf8")

    def __get__(self, obj, objtype = None):
        # Have we been called before?
        return BytesIO(self.content)

class FileObjManager:
    caffeine_2d_sdf = FileobjDescriptor(CAFFEINE_2D_SDF)
    caffeine_3d_sdf = FileobjDescriptor(CAFFEINE_3D_SDF)
    caffeine_smi = FileobjDescriptor(CAFFEINE_SMI)
    
    aspirin_2d_sdf = FileobjDescriptor(ASPIRIN_2D_SDF)
    aspirin_3d_sdf = FileobjDescriptor(ASPIRIN_3D_SDF)
    aspirin_smi = FileobjDescriptor(ASPIRIN_SMI)

    chebi_1148_sdf = FileobjDescriptor(CHEBI_1148_SDF)
    
fileobj_manager = FileObjManager()

# ========================================================
# Base class to test from_file() and from_file_obj()
# ========================================================


class BaseFromFileIO:
    # == Test variations of "sdf" and "mol"

    @pytest.mark.parametrize("file_format", ("SDF", "sdf", "sdF", "MOL", "mol"))
    def test_from_file_sdf_ignores_file_format_case(self, file_format):
        mols = self.toolkit_wrapper.from_file(file_manager.caffeine_2d_sdf, file_format)
        self._test_from_sdf_ignores_file_format_case(mols)

    @pytest.mark.parametrize("file_format", ("SDF", "sdf", "sdF", "MOL", "mol"))
    def test_from_file_obj_sdf_ignores_file_format_case(self, file_format):
        mols = self.toolkit_wrapper.from_file_obj(fileobj_manager.caffeine_2d_sdf, file_format)
        self._test_from_sdf_ignores_file_format_case(mols)

    def _test_from_sdf_ignores_file_format_case(self, mols):
        assert len(mols) == 1
        mol = mols[0]
        assert mol.name == "caffeine"

    # == Test variations of "smi" format
        
    @pytest.mark.parametrize("file_format", ("SMI", "smi", "sMi"))
    def test_from_file_smi_ignores_file_format_case(self, file_format):
        mols = self.toolkit_wrapper.from_file(file_manager.caffeine_smi, file_format)
        self._test_from_smi_ignores_file_format_case(mols)

    @pytest.mark.parametrize("file_format", ("SMI", "smi", "sMi"))
    def test_from_fileobj_smi_ignores_file_format_case(self, file_format):
        mols = self.toolkit_wrapper.from_file_obj(fileobj_manager.caffeine_smi, file_format)
        self._test_from_smi_ignores_file_format_case(mols)

    def _test_from_smi_ignores_file_format_case(self, mols):
        assert len(mols) == 1
        mol = mols[0]
        assert mol.name == "CHEMBL113"
        assert mol.n_atoms == 24
        assert mol.n_bonds == 25
        assert mol.n_conformers == 0

    # == Test format "qwe" raises an exception

    def test_from_file_qwe_format_raises_exception(self):
        with pytest.raises(ValueError, match = "Unsupported file format: QWE"):
            mols = self.toolkit_wrapper.from_file(file_manager.caffeine_2d_sdf, file_format = "qwe")
            
    def test_from_file_obj_qwe_format_raises_exception(self):
        with pytest.raises(ValueError, match = "Unsupported file format: QWE"):
            mols = self.toolkit_wrapper.from_file_obj(fileobj_manager.caffeine_2d_sdf, file_format = "qwe")
            
        
    # == Test reading SDF 2D coordinates
        
    def test_from_file_2D_sdf_coords(self):
        mol = self.toolkit_wrapper.from_file(file_manager.caffeine_2d_sdf, "sdf")[0]
        self._test_2D_sdf_coords(mol)

    def test_from_file_obj_2D_sdf_coords(self):
        mol = self.toolkit_wrapper.from_file_obj(fileobj_manager.caffeine_2d_sdf, "sdf")[0]
        self._test_2D_sdf_coords(mol)
        
    def _test_2D_sdf_coords(self, mol):
        assert mol.n_atoms == 24
        assert mol.n_bonds == 25
        assert mol.n_conformers == 1
        conformer = mol.conformers[0]
        assert conformer.shape == (24, 3)
        assert conformer.unit == unit.angstrom

        # Beyond this are the hydrogens, which are added by algorithm.
        assert_allclose(conformer[:CAFFEINE_2D_COORDS.shape[0]], CAFFEINE_2D_COORDS)
        
    # == Test reading SDF 3D coordinates
    
    def test_from_file_3D_sdf_keeps_hydrogens(self):
        mol = self.toolkit_wrapper.from_file(file_manager.caffeine_3d_sdf, "sdf")[0]
        self._test_3D_sdf_keeps_hydrogens(mol)
        
    def test_from_file_obj_3D_sdf_keeps_hydrogens(self):
        mol = self.toolkit_wrapper.from_file_obj(fileobj_manager.caffeine_3d_sdf, "sdf")[0]
        self._test_3D_sdf_keeps_hydrogens(mol)

    def _test_3D_sdf_keeps_hydrogens(self, mol):
        assert mol.n_atoms == 24
        assert mol.n_bonds == 25
        assert mol.n_conformers == 1
        conformer = mol.conformers[0]
        assert conformer.shape == (24, 3)
        assert conformer.unit == unit.angstrom

        # All hydrogens are explicit in the file
        assert_allclose(conformer, CAFFEINE_3D_COORDS)

    # == Test allow_undefined_stereo

    def test_from_file_with_undefined_stereo(self):
        with pytest.raises(
                exceptions.UndefinedStereochemistryError,
                match = f"Unable to make OFFMol from {self.tk_mol_name}: {self.tk_mol_name} has unspecified stereochemistry",
                ):
            self.toolkit_wrapper.from_file(file_manager.chebi_1148_sdf, "sdf")

    def test_from_file_obj_with_undefined_stereo(self):
        with pytest.raises(
                exceptions.UndefinedStereochemistryError,
                match = f"Unable to make OFFMol from {self.tk_mol_name}: {self.tk_mol_name} has unspecified stereochemistry",
                ):
            self.toolkit_wrapper.from_file_obj(fileobj_manager.chebi_1148_sdf, "sdf")

    def test_from_file_with_undefined_stereo_allowed(self):
        self.toolkit_wrapper.from_file(file_manager.chebi_1148_sdf, "sdf", allow_undefined_stereo=True)[0]

    def test_from_file_obj_with_undefined_stereo_allowed(self):
        self.toolkit_wrapper.from_file_obj(fileobj_manager.chebi_1148_sdf, "sdf", allow_undefined_stereo=True)[0]


    # == Test passing in a user-defined _cls

    @pytest.mark.parametrize("name,file_format", [("caffeine_2d_sdf", "SDF"), ("caffeine_smi", "SMI")])
    def test_from_file_handles_cls(self, name, file_format):
        filename = getattr(file_manager, name)
        mol = self.toolkit_wrapper.from_file(filename, file_format, _cls = SingingMolecule)[0]
        mol.sing()

    @pytest.mark.parametrize("name,file_format", [("caffeine_2d_sdf", "SDF"), ("caffeine_smi", "SMI")])
    def test_from_file_handles_cls(self, name, file_format):
        file_obj = getattr(fileobj_manager, name)
        mol = self.toolkit_wrapper.from_file_obj(file_obj, file_format, _cls = SingingMolecule)[0]
        mol.sing()
        
    
    ## def to_file_obj(self):
    ##     molecule, file_obj, file_format

    ## def to_file(self):
    ##     pass

@pytest.fixture(scope="class")
def init_toolkit(request):
    request.cls.toolkit_wrapper = request.cls.toolkit_wrapper_class()

@pytest.mark.usefixtures("init_toolkit")
class TestOpenEyeToolkitFromFileIO(BaseFromFileIO):
    toolkit_wrapper_class = OpenEyeToolkitWrapper
    tk_mol_name = "OEMol"

@pytest.mark.usefixtures("init_toolkit")
class TestRDKitToolkitFromFileIO(BaseFromFileIO):
    toolkit_wrapper_class = RDKitToolkitWrapper
    tk_mol_name = "RDMol"

    def test_from_file_obj_smi_supports_stringio(self):
        # Backwards compability. Older versions of OpenFF supported
        # string-based file objects, but not file-based ones.
        with open(file_manager.caffeine_smi) as file_obj:
            mol = self.toolkit_wrapper.from_file_obj(file_obj, "SMI")[0]
        assert mol.name == "CHEMBL113"

# ========================================================
# Base class to test to_file() and to_file_obj()
# ========================================================

@pytest.fixture(scope="class")
def tmpdir(request):
    request.cls.tmpdir = tmpdir = tempfile.TemporaryDirectory()
    with tmpdir:
        yield
    request.cls.tmpdir = None


def assert_is_ethanol_sdf(f):
    assert f.readline() == "ethanol\n" # title line
    f.readline() # ignore next two lines
    f.readline()
    assert f.readline()[:6] == "  9  8" # check 9 atoms, 8 bonds

def assert_is_ethanol_smi(f):
    line = f.readline()
    assert "ethanol" in line
    assert_is_ethanol_smiles(line.split()[0])
    assert not f.readline(), "should only have one line"

def assert_is_ethanol_smiles(smiles):
    assert smiles.count("[H]") == 6
    assert smiles.count("O") == 1
    assert smiles.count("C") == 2
    
        
class BaseToFileIO:
    def get_tmpfile(self, name):
        return os.path.join(self.tmpdir.name, name)
    
    # == Test variations of "sdf" and "mol"

    @pytest.mark.parametrize("format_name", ["SDF", "sdf", "sDf", "mol", "MOL"])
    def test_to_file_sdf(self, format_name):
        filename = self.get_tmpfile("abc.xyz")
        self.toolkit_wrapper.to_file(ETHANOL, filename, format_name)
        with open(filename) as f:
            assert_is_ethanol_sdf(f)

    @pytest.mark.parametrize("format_name", ["SDF", "sdf", "sDf", "mol", "MOL"])
    def test_to_file_obj_sdf(self, format_name):
        f = StringIO()
        self.toolkit_wrapper.to_file_obj(ETHANOL, f, format_name)
        f.seek(0)
        assert_is_ethanol_sdf(f)

    def test_to_file_obj_sdf_with_bytesio(self):
        f = BytesIO()
        with pytest.raises(ValueError, match = "Need a text mode file object like StringIO or a file opened with mode 't'"):
            self.toolkit_wrapper.to_file_obj(ETHANOL, f, "sdf")

    # === Test variations of "smi"

    @pytest.mark.parametrize("format_name", ["SMI", "smi", "sMi"])
    def test_to_file_smi(self, format_name):
        filename = self.get_tmpfile("abc.xyz")
        self.toolkit_wrapper.to_file(ETHANOL, filename, format_name)
        with open(filename) as f:
            assert_is_ethanol_smi(f)

    @pytest.mark.parametrize("format_name", ["SMI", "smi", "sMi"])
    def test_to_file_obj_smi(self, format_name):
        f = StringIO()
        self.toolkit_wrapper.to_file_obj(ETHANOL, f, format_name)
        f.seek(0)
        assert_is_ethanol_smi(f)

    def test_to_file_obj_smi_with_bytesio(self):
        f = BytesIO()
        with pytest.raises(ValueError, match = "Need a text mode file object like StringIO or a file opened with mode 't'"):
            self.toolkit_wrapper.to_file_obj(ETHANOL, f, "smi")

    # === Test unsupported format

    def test_to_file_qwe_format_raises_exception(self):
        with tempfile.NamedTemporaryFile(suffix=".smi") as fileobj:
            with pytest.raises(ValueError, match=f"Unsupported file format: QWE"):
                self.toolkit_wrapper.to_file(ETHANOL, fileobj.name, "QWE")

@pytest.mark.usefixtures("init_toolkit", "tmpdir")
class TestOpenEyeToolkitToFileIO(BaseToFileIO):
    toolkit_wrapper_class = OpenEyeToolkitWrapper

@pytest.mark.usefixtures("init_toolkit", "tmpdir")
class TestRDKitToolkitToFileIO(BaseToFileIO):
    toolkit_wrapper_class = RDKitToolkitWrapper


# ========================================================
# Base class to test SMILES parsing
# ========================================================

class BaseSmiles:
    def test_parse_methane_with_implicit_Hs(self):
        mol = self.toolkit_wrapper.from_smiles("C")
        # adds hydrogens
        assert mol.n_atoms == 5
        assert mol.n_bonds == 4
        assert mol.partial_charges is None
        
    def test_parse_methane_with_implicit_Hs_but_said_they_are_explicit(self):
        with pytest.raises(
                ValueError,
                match = (
                    "'hydrogens_are_explicit' was specified as True, but (OpenEye|RDKit) [Tt]oolkit interpreted SMILES "
                    "[^ ]+ as having implicit hydrogen. If this SMILES is intended to express all explicit hydrogens "
                    f"in the molecule, then you should construct the desired molecule as an {self.tk_mol_name}"
                    )):
            mol = self.toolkit_wrapper.from_smiles("C", hydrogens_are_explicit=True)
    
    def test_parse_methane_with_explicit_Hs(self):
        mol = self.toolkit_wrapper.from_smiles("[C]([H])([H])([H])([H])")
        # add hydrogens
        assert mol.n_atoms == 5
        assert mol.n_bonds == 4
        assert molecule.partial_charges is None
    
    def test_parse_methane_with_explicit_Hs(self):
        mol = self.toolkit_wrapper.from_smiles("[C]([H])([H])([H])([H])", hydrogens_are_explicit=True)
        # add hydrogens
        assert mol.n_atoms == 5
        assert mol.n_bonds == 4

    def test_parse_bad_smiles(self):
        with pytest.raises(ValueError, match="Unable to parse the SMILES string"):
            mol = self.toolkit_wrapper.from_smiles("QWERT")

    ### Copied from test_toolkits.py

    @pytest.mark.parametrize(
        "title, smiles",
        [("unspec_chiral_smiles", r"C\C(F)=C(/F)CC(C)(Cl)Br"),
         ("spec_chiral_smiles", r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"),
         ("unspec_db_smiles", r"CC(F)=C(F)C[C@@](C)(Cl)Br"),
         ("spec_db_smiles", r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"),])
    def test_smiles_missing_stereochemistry(self, title, smiles):
        if "unspec" in title:
            with pytest.raises(exceptions.UndefinedStereochemistryError):
                self.toolkit_wrapper.from_smiles(smiles)
            mol = self.toolkit_wrapper.from_smiles(smiles, allow_undefined_stereo=True)
        else:
            mol = self.toolkit_wrapper.from_smiles(smiles)
        assert mol.n_atoms == 18

    @pytest.mark.parametrize(
        "smiles, expected_map", [("[Cl:1][H]", {0: 1}), ("[Cl:1][H:2]", {0: 1, 1: 2})]
    )
    def test_from_smiles_with_atom_map(self, smiles, expected_map):
        mol = self.toolkit_wrapper.from_smiles(smiles)
        assert mol.properties["atom_map"] == expected_map
        
    def test_ethanol_to_smiles(self):
        smiles = self.toolkit_wrapper.to_smiles(ETHANOL)
        assert_is_ethanol_smiles(smiles)

    def test_round_trip_ethanol_to_smiles_adds_hydrogens(self):
        mol = self.toolkit_wrapper.from_smiles("CCO")
        smiles = self.toolkit_wrapper.to_smiles(mol)
        assert_is_ethanol_smiles(smiles)


            
@pytest.mark.usefixtures("init_toolkit")
class TestOpenEyeToolkitSmiles(BaseSmiles):
    toolkit_wrapper_class = OpenEyeToolkitWrapper
    tk_mol_name = "OEMol"

@pytest.mark.usefixtures("init_toolkit")
class TestRDKitToolkitSmiles(BaseSmiles):
    toolkit_wrapper_class = RDKitToolkitWrapper
    tk_mol_name = "RDMol"

    
    
if __name__ == "__main__":
    sys.exit(pytest.main([__file__] + sys.argv[1:]))
