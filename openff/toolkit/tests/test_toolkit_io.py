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
import pathlib
import sys
import tempfile
from io import BytesIO, StringIO

import numpy as np
import pytest
from numpy.testing import assert_allclose
from simtk import unit

from openff.toolkit.tests import create_molecules
from openff.toolkit.tests.utils import requires_openeye, requires_rdkit
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper, exceptions

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

CAFFEINE_2D_COORDS = (
    np.array(
        [
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
        ],
        np.double,
    )
    * unit.angstrom
)


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

CAFFEINE_3D_COORDS = (
    np.array(
        [
            (0.4700, 2.5688, 0.0006),
            (-3.1271, -0.4436, -0.0003),
            (-0.9686, -1.3125, 0.0000),
            (2.2182, 0.1412, -0.0003),
            (-1.3477, 1.0797, -0.0001),
            (1.4119, -1.9372, 0.0002),
            (0.8579, 0.2592, -0.0008),
            (0.3897, -1.0264, -0.0004),
            (0.0307, 1.4220, -0.0006),
            (-1.9061, -0.2495, -0.0004),
            (2.5032, -1.1998, 0.0003),
            (-1.4276, -2.6960, 0.0008),
            (3.1926, 1.2061, 0.0003),
            (-2.2969, 2.1881, 0.0007),
            (3.5163, -1.5787, 0.0008),
            (-1.0451, -3.1973, -0.8937),
            (-2.5186, -2.7596, 0.0011),
            (-1.0447, -3.1963, 0.8957),
            (4.1992, 0.7801, 0.0002),
            (3.0468, 1.8092, -0.8992),
            (3.0466, 1.8083, 0.9004),
            (-1.8087, 3.1651, -0.0003),
            (-2.9322, 2.1027, 0.8881),
            (-2.9346, 2.1021, -0.8849),
        ],
        np.double,
    )
    * unit.angstrom
)

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

TWO_MOLS_SDF = CAFFEINE_2D_SDF + ASPIRIN_2D_SDF
TWO_MOLS_SMI = CAFFEINE_SMI + ASPIRIN_SMI

# Force invalid records to ensure the invalid one is skipped
THREE_MOLS_SDF = (
    CAFFEINE_2D_SDF + CAFFEINE_3D_SDF.replace("24 25", "20 29") + ASPIRIN_2D_SDF
)
THREE_MOLS_SMI = CAFFEINE_SMI + "Q" + CAFFEINE_SMI + ASPIRIN_SMI

# ========================================================
# Used to test for records with disconnected molecules
# ========================================================

TABLE_SALT_SDF = """\
table salt
     RDKit

  2  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Na  0  0  0  0  0 15  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0 15  0  0  0  0  0  0
M  CHG  2   1   1   2  -1
M  END
$$$$
"""

# ========================================================
# Used to test for records with an R-group
# ========================================================

METHYL_GROUP_SDF = """\
methyl group
     RDKit

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 R   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""

# ========================================================
# Used to test for records with an unexpected bond type
# ========================================================

CHEBI_52729_SDF = """\

Marvin  07090914422D          

 78 91  0  0  0  0            999 V2000
    4.4196   -7.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1341   -6.8652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8486   -7.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5631   -6.8652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2775   -7.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5631   -6.0402    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1341   -6.0402    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7522   -6.7927    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    3.0848   -7.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3397   -8.0623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1647   -8.0623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.9920   -6.8652    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    7.2775   -8.1027    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    7.9920   -7.6902    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.6693   -5.4838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0818   -4.7694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6693   -4.0549    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0818   -3.3404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6693   -2.6259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9068   -3.3404    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9068   -4.7694    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1543   -6.1513    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.6693   -6.8187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8848   -6.5638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8848   -5.7388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0818   -1.9115    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.8443   -2.6259    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.2568   -1.9115    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    8.9428   -5.4589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5303   -4.7444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9428   -4.0299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5303   -3.3155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9428   -2.6010    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7053   -3.3155    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7053   -4.7444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4579   -6.1263    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    8.9428   -6.7937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.7274   -6.5388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.7274   -5.7138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5303   -1.8865    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    9.7678   -2.6010    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    9.3553   -1.8865    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    5.7161    1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0016    1.2375    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.0016    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7161    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.4305    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4305    1.2375    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.7161    2.4750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7161    4.1251    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0016    3.7125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0016    2.8875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4305    2.8875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4305    3.7125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7161    4.9501    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.4305    5.3626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4305    6.1876    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0016    5.3626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2871    4.9501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4972    0.2380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0268   -0.4103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5260   -1.0966    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3049   -0.8248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.2871    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3407    1.0772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8561    0.4095    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1922   -0.3440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1612    0.9915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.9349    0.2380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4053   -0.4103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.9062   -1.0966    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1272   -0.8248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.1450    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.0915    1.0772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5760    0.4095    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2400   -0.3440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2710    0.9915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8929   -4.3018    0.0000 Eu  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  4  6  1  0  0  0  0
  2  7  2  0  0  0  0
 11  1  2  0  0  0  0
  8  1  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
  5 12  1  0  0  0  0
  5 13  1  0  0  0  0
  5 14  1  0  0  0  0
 15 16  1  0  0  0  0
 25 15  2  0  0  0  0
 22 15  1  0  0  0  0
 16 17  1  0  0  0  0
 16 21  2  0  0  0  0
 17 18  2  0  0  0  0
 18 19  1  0  0  0  0
 18 20  1  0  0  0  0
 19 26  1  0  0  0  0
 19 27  1  0  0  0  0
 19 28  1  0  0  0  0
 22 23  1  0  0  0  0
 23 24  2  0  0  0  0
 24 25  1  0  0  0  0
 29 30  1  0  0  0  0
 39 29  2  0  0  0  0
 36 29  1  0  0  0  0
 30 31  1  0  0  0  0
 30 35  2  0  0  0  0
 31 32  2  0  0  0  0
 32 33  1  0  0  0  0
 32 34  1  0  0  0  0
 33 40  1  0  0  0  0
 33 41  1  0  0  0  0
 33 42  1  0  0  0  0
 36 37  1  0  0  0  0
 37 38  2  0  0  0  0
 38 39  1  0  0  0  0
 43 44  2  0  0  0  0
 43 48  1  0  0  0  0
 44 45  1  0  0  0  0
 45 46  2  0  0  0  0
 46 47  1  0  0  0  0
 47 48  2  0  0  0  0
 43 49  1  0  0  0  0
 52 49  2  0  0  0  0
 49 53  1  0  0  0  0
 50 51  2  0  0  0  0
 50 54  1  0  0  0  0
 51 52  1  0  0  0  0
 53 54  2  0  0  0  0
 50 55  1  0  0  0  0
 55 56  1  0  0  0  0
 56 57  1  0  0  0  0
 55 58  1  0  0  0  0
 58 59  1  0  0  0  0
 47 73  1  0  0  0  0
 45 64  1  0  0  0  0
 60 64  1  0  0  0  0
 61 62  1  0  0  0  0
 62 63  2  0  0  0  0
 63 64  1  0  0  0  0
 60 68  2  0  0  0  0
 61 60  1  0  0  0  0
 67 61  2  0  0  0  0
 65 66  2  0  0  0  0
 65 68  1  0  0  0  0
 66 67  1  0  0  0  0
 69 73  1  0  0  0  0
 69 77  2  0  0  0  0
 70 69  1  0  0  0  0
 70 71  1  0  0  0  0
 76 70  2  0  0  0  0
 71 72  2  0  0  0  0
 72 73  1  0  0  0  0
 74 75  2  0  0  0  0
 74 77  1  0  0  0  0
 75 76  1  0  0  0  0
 20 78  1  0  0  0  0
 34 78  1  0  0  0  0
  6 78  1  0  0  0  0
 21 78  8  0  0  0  0
  7 78  8  0  0  0  0
 35 78  8  0  0  0  0
 46 78  8  0  0  0  0
 72 78  8  0  0  0  0
 63 78  8  0  0  0  0
M  STY  6   1 DAT   2 DAT   3 DAT   4 DAT   5 DAT   6 DAT
M  SAL   1  2  21  78
M  SDT   1 MRV_COORDINATE_BOND_TYPE                              
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 86
M  SAL   2  2   7  78
M  SDT   2 MRV_COORDINATE_BOND_TYPE                              
M  SDD   2     0.0000    0.0000    DR    ALL  0       0  
M  SED   2 87
M  SAL   3  2  35  78
M  SDT   3 MRV_COORDINATE_BOND_TYPE                              
M  SDD   3     0.0000    0.0000    DR    ALL  0       0  
M  SED   3 88
M  SAL   4  2  46  78
M  SDT   4 MRV_COORDINATE_BOND_TYPE                              
M  SDD   4     0.0000    0.0000    DR    ALL  0       0  
M  SED   4 89
M  SAL   5  2  72  78
M  SDT   5 MRV_COORDINATE_BOND_TYPE                              
M  SDD   5     0.0000    0.0000    DR    ALL  0       0  
M  SED   5 90
M  SAL   6  2  63  78
M  SDT   6 MRV_COORDINATE_BOND_TYPE                              
M  SDD   6     0.0000    0.0000    DR    ALL  0       0  
M  SED   6 91
M  END
> <ID>
CHEBI:52729

> <NAME>
Eu(tta)3DEADIT

> <DEFINITION>
A europium coordination entity composed of europium(III) coordinated to 4-[4,6-di(1H-indazol-1-yl)-1,3,5-triazin-2-yl]-N,N-diethylaniline and three 4,4,4-trifluoro-1-(thiophen-2-yl)butane-1,3-dione units.

> <STAR>
3

> <IUPAC_NAME>
{4-[4,6-di(1H-indazol-1-yl-kappaN(2))-1,3,5-triazin-2-yl-kappaN(5)]-N,N-diethylaniline}{tris[4,4,4-trifluoro-3-(hydroxy-kappaO)-1-(thiophen-2-yl)but-2-en-1-onato-kappaO]}europium

> <FORMULA>
C51H36EuF9N8O6S3

> <MASS>
1276.02600

> <MONOISOTOPIC_MASS>
1276.09885

> <CHARGE>
0

> <SMILES>
CCN(CC)c1ccc(cc1)-c1nc([nH]c(n1)-n1[nH]cc2ccccc12)-n1[nH]cc2ccccc12.FC(F)(F)C(\\O[Eu](O\\C(=C/C(=O)c1cccs1)C(F)(F)F)O\\C(=C/C(=O)c1cccs1)C(F)(F)F)=C\\C(=O)c1cccs1

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


# This code manages a temporary directory.
# Data files can be found inside that directory.
# They are created when needed, and stored for later use.

# These descriptors make it possible to ask for
# `file_manager.name_sdf' and get the path to the named file inside
# the directory.


class FilenameDescriptor:
    def __init__(self, content):
        self.content = content

    def __set_name__(self, owner, name):
        # Figure out the descriptor name, and
        # the variable name used to store the filename.
        self.private_name = f"_{name}_filename"
        basename, _, ext = name.rpartition("_")
        self.filename = f"{basename}.{ext}"

    def __get__(self, obj, objtype=None):
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

    # Used to test that file_format overrides any automatic file detection.
    caffeine_not_smi = FilenameDescriptor(CAFFEINE_2D_SDF)
    caffeine_not_sdf = FilenameDescriptor(CAFFEINE_SMI)

    ## aspirin_2d_sdf = FilenameDescriptor(ASPIRIN_2D_SDF)
    ## aspirin_3d_sdf = FilenameDescriptor(ASPIRIN_3D_SDF)
    ## aspirin_smi = FilenameDescriptor(ASPIRIN_SMI)

    two_mols_sdf = FilenameDescriptor(TWO_MOLS_SDF)
    two_mols_smi = FilenameDescriptor(TWO_MOLS_SMI)

    three_mols_sdf = FilenameDescriptor(THREE_MOLS_SDF)
    three_mols_smi = FilenameDescriptor(THREE_MOLS_SMI)

    chebi_1148_sdf = FilenameDescriptor(CHEBI_1148_SDF)


file_manager = FileManager(tempfile.TemporaryDirectory())


# This is similar to the FilenameDescriptor except that it returns a
# BytesIO.


class FileobjDescriptor:
    def __init__(self, content):
        self.content = content.encode("utf8")

    def __get__(self, obj, objtype=None):
        # Have we been called before?
        return BytesIO(self.content)


class FileObjManager:
    caffeine_2d_sdf = FileobjDescriptor(CAFFEINE_2D_SDF)
    caffeine_3d_sdf = FileobjDescriptor(CAFFEINE_3D_SDF)
    caffeine_smi = FileobjDescriptor(CAFFEINE_SMI)

    ## aspirin_2d_sdf = FileobjDescriptor(ASPIRIN_2D_SDF)
    ## aspirin_3d_sdf = FileobjDescriptor(ASPIRIN_3D_SDF)
    ## aspirin_smi = FileobjDescriptor(ASPIRIN_SMI)

    two_mols_sdf = FileobjDescriptor(TWO_MOLS_SDF)
    two_mols_smi = FileobjDescriptor(TWO_MOLS_SMI)

    three_mols_sdf = FileobjDescriptor(THREE_MOLS_SDF)
    three_mols_smi = FileobjDescriptor(THREE_MOLS_SMI)

    chebi_1148_sdf = FileobjDescriptor(CHEBI_1148_SDF)

    table_salt_sdf = FileobjDescriptor(TABLE_SALT_SDF)
    methyl_group_sdf = FileobjDescriptor(METHYL_GROUP_SDF)
    chebi_52729_sdf = FileobjDescriptor(CHEBI_52729_SDF)


file_obj_manager = FileObjManager()

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
        mols = self.toolkit_wrapper.from_file_obj(
            file_obj_manager.caffeine_2d_sdf, file_format
        )
        self._test_from_sdf_ignores_file_format_case(mols)

    def _test_from_sdf_ignores_file_format_case(self, mols):
        assert len(mols) == 1
        mol = mols[0]
        assert mol.name == "caffeine"

    # == Test reading two molecules from an SDF

    def test_from_file_sdf_two_molecules(self):
        mols = self.toolkit_wrapper.from_file(file_manager.two_mols_sdf, "SDF")
        self._test_from_sdf_two_molecules(mols)

    def test_from_file_obj_sdf_two_molecules(self):
        mols = self.toolkit_wrapper.from_file_obj(file_obj_manager.two_mols_sdf, "SDF")
        self._test_from_sdf_two_molecules(mols)

    def _test_from_sdf_two_molecules(self, mols):
        assert len(mols) == 2
        assert mols[0].name == "caffeine"
        assert mols[1].name == "aspirin"

    # == Test skipping an error molecule from an SDF

    def test_from_file_sdf_three_molecules(self):
        mols = self.toolkit_wrapper.from_file(file_manager.three_mols_sdf, "SDF")
        self._test_from_sdf_two_molecules(mols)

    def test_from_file_obj_sdf_three_molecules(self):
        mols = self.toolkit_wrapper.from_file_obj(
            file_obj_manager.three_mols_sdf, "SDF"
        )
        self._test_from_sdf_two_molecules(mols)

    # == Test variations of "smi" format

    @pytest.mark.parametrize("file_format", ("SMI", "smi", "sMi"))
    def test_from_file_smi_ignores_file_format_case(self, file_format):
        mols = self.toolkit_wrapper.from_file(file_manager.caffeine_smi, file_format)
        self._test_from_smi_ignores_file_format_case(mols)

    @pytest.mark.parametrize("file_format", ("SMI", "smi", "sMi"))
    def test_from_fileobj_smi_ignores_file_format_case(self, file_format):
        mols = self.toolkit_wrapper.from_file_obj(
            file_obj_manager.caffeine_smi, file_format
        )
        self._test_from_smi_ignores_file_format_case(mols)

    def _test_from_smi_ignores_file_format_case(self, mols):
        assert len(mols) == 1
        mol = mols[0]
        assert mol.name == "CHEMBL113"
        assert mol.n_atoms == 24
        assert mol.n_bonds == 25
        assert mol.n_conformers == 0

    # == Test reading two molecules from a SMILES file

    def test_from_file_smi_two_molecules(self):
        mols = self.toolkit_wrapper.from_file(file_manager.two_mols_smi, "SMI")
        self._test_from_smi_two_molecules(mols)

    def test_from_file_obj_smi_two_molecules(self):
        mols = self.toolkit_wrapper.from_file_obj(file_obj_manager.two_mols_smi, "SMI")
        self._test_from_smi_two_molecules(mols)

    def _test_from_smi_two_molecules(self, mols):
        assert len(mols) == 2
        assert mols[0].name == "CHEMBL113"
        assert mols[1].name == "ASPIRIN"

    # == Test skipping an error molecule from a SMI

    def test_from_file_smi_three_molecules(self):
        mols = self.toolkit_wrapper.from_file(file_manager.three_mols_smi, "SMI")
        self._test_from_smi_two_molecules(mols)

    def test_from_file_obj_smi_three_molecules(self):
        mols = self.toolkit_wrapper.from_file_obj(
            file_obj_manager.three_mols_smi, "SMI"
        )
        self._test_from_smi_two_molecules(mols)

    # == Test that file_format overrides the extension

    def test_from_file_sdf_ignores_filename_extension(self):
        mols = self.toolkit_wrapper.from_file(file_manager.caffeine_not_smi, "sdf")
        self._test_from_sdf_ignores_file_format_case(mols)

    def test_from_file_smi_ignores_filename_extension(self):
        mols = self.toolkit_wrapper.from_file(file_manager.caffeine_not_sdf, "smi")
        self._test_from_smi_ignores_file_format_case(mols)

    # == Test format "qwe" raises an exception

    def test_from_file_qwe_format_raises_exception(self):
        with pytest.raises(ValueError, match="Unsupported file format: QWE"):
            mols = self.toolkit_wrapper.from_file(
                file_manager.caffeine_2d_sdf, file_format="qwe"
            )

    def test_from_file_obj_qwe_format_raises_exception(self):
        with pytest.raises(ValueError, match="Unsupported file format: QWE"):
            mols = self.toolkit_wrapper.from_file_obj(
                file_obj_manager.caffeine_2d_sdf, file_format="qwe"
            )

    # == Test when a file doesn't exist

    @pytest.mark.parametrize("file_format", ["smi", "sdf", "mol"])
    def test_from_file_when_the_filename_does_not_exist(self, file_format):
        filename = f"/qwelkjpath/to/file/that/does/not/exist/asdflkjasdf.{file_format}"
        with pytest.raises(OSError):
            self.toolkit_wrapper.from_file(filename, file_format=file_format)

    # == Test reading SDF 2D coordinates

    def test_from_file_2D_sdf_coords(self):
        mol = self.toolkit_wrapper.from_file(file_manager.caffeine_2d_sdf, "sdf")[0]
        self._test_2D_sdf_coords(mol)

    def test_from_file_obj_2D_sdf_coords(self):
        mol = self.toolkit_wrapper.from_file_obj(
            file_obj_manager.caffeine_2d_sdf, "sdf"
        )[0]
        self._test_2D_sdf_coords(mol)

    def _test_2D_sdf_coords(self, mol):
        assert mol.n_atoms == 24
        assert mol.n_bonds == 25
        assert mol.n_conformers == 1
        conformer = mol.conformers[0]
        assert conformer.shape == (24, 3)
        assert conformer.unit == unit.angstrom

        # Beyond this are the hydrogens, which are added by algorithm.
        assert_allclose(conformer[: CAFFEINE_2D_COORDS.shape[0]], CAFFEINE_2D_COORDS)

    # == Test reading SDF 3D coordinates

    def test_from_file_3D_sdf_keeps_hydrogens(self):
        mol = self.toolkit_wrapper.from_file(file_manager.caffeine_3d_sdf, "sdf")[0]
        self._test_3D_sdf_keeps_hydrogens(mol)

    def test_from_file_obj_3D_sdf_keeps_hydrogens(self):
        mol = self.toolkit_wrapper.from_file_obj(
            file_obj_manager.caffeine_3d_sdf, "sdf"
        )[0]
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
            match=f"Unable to make OFFMol from {self.tk_mol_name}: {self.tk_mol_name} has unspecified stereochemistry",
        ):
            self.toolkit_wrapper.from_file(file_manager.chebi_1148_sdf, "sdf")

    def test_from_file_obj_with_undefined_stereo(self):
        with pytest.raises(
            exceptions.UndefinedStereochemistryError,
            match=f"Unable to make OFFMol from {self.tk_mol_name}: {self.tk_mol_name} has unspecified stereochemistry",
        ):
            self.toolkit_wrapper.from_file_obj(file_obj_manager.chebi_1148_sdf, "sdf")

    def test_from_file_with_undefined_stereo_allowed(self):
        self.toolkit_wrapper.from_file(
            file_manager.chebi_1148_sdf, "sdf", allow_undefined_stereo=True
        )[0]

    def test_from_file_obj_with_undefined_stereo_allowed(self):
        self.toolkit_wrapper.from_file_obj(
            file_obj_manager.chebi_1148_sdf, "sdf", allow_undefined_stereo=True
        )[0]

    # == Test passing in a user-defined _cls

    @pytest.mark.parametrize(
        "name,file_format", [("caffeine_2d_sdf", "SDF"), ("caffeine_smi", "SMI")]
    )
    def test_from_file_handles_cls(self, name, file_format):
        filename = getattr(file_manager, name)
        mol = self.toolkit_wrapper.from_file(
            filename, file_format, _cls=SingingMolecule
        )[0]
        mol.sing()

    @pytest.mark.parametrize(
        "name,file_format", [("caffeine_2d_sdf", "SDF"), ("caffeine_smi", "SMI")]
    )
    def test_from_file_obj_handles_cls(self, name, file_format):
        file_obj = getattr(file_obj_manager, name)
        mol = self.toolkit_wrapper.from_file_obj(
            file_obj, file_format, _cls=SingingMolecule
        )[0]
        mol.sing()


# Create the appropriate toolkit wrapper for each test class instance.
#
# We could create them at the module level, but this would require
# re-doing the license checks that @requires_openeye and
# @requires_rdkit already do. Instead, let those decorators handle the
# license check, and use pytest's fixture system to create the
# instance given the appropriate wrapper class defined in the
# subclass.
#
# The wrappers are stateless, so specify it as class fixture, instead
# of creating a new instance for each test


@pytest.fixture(scope="class")
def init_toolkit(request):
    request.cls.toolkit_wrapper = request.cls.toolkit_wrapper_class()


@pytest.mark.usefixtures("init_toolkit")
@requires_openeye
class TestOpenEyeToolkitFromFileIO(BaseFromFileIO):
    toolkit_wrapper_class = OpenEyeToolkitWrapper
    tk_mol_name = "OEMol"


@pytest.mark.usefixtures("init_toolkit")
@requires_rdkit
class TestRDKitToolkitFromFileIO(BaseFromFileIO):
    toolkit_wrapper_class = RDKitToolkitWrapper
    tk_mol_name = "RDMol"

    def test_from_file_obj_smi_supports_stringio(self):
        # Test the backwards compatibility support for
        # passing in file objects open in "t"ext mode.
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
    assert f.readline() == "ethanol\n"  # title line
    f.readline()  # ignore next two lines
    f.readline()
    assert f.readline()[:6] == "  9  8"  # check 9 atoms, 8 bonds


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
        with pytest.raises(
            ValueError,
            match="Need a text mode file object like StringIO or a file opened with mode 't'",
        ):
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
        with pytest.raises(
            ValueError,
            match="Need a text mode file object like StringIO or a file opened with mode 't'",
        ):
            self.toolkit_wrapper.to_file_obj(ETHANOL, f, "smi")

    # === Test format "qwe" raises an exception

    def test_to_file_qwe_format_raises_exception(self):
        with tempfile.NamedTemporaryFile(suffix=".smi") as fileobj:
            with pytest.raises(ValueError, match=f"Unsupported file format: QWE"):
                self.toolkit_wrapper.to_file(ETHANOL, fileobj.name, "QWE")

    # === Test writing to a file that does not exist

    @pytest.mark.parametrize("format_name", ["smi", "sdf", "mol"])
    def test_to_file_when_the_file_does_not_exist(self, format_name):
        filename = self.get_tmpfile("does/not/exist.smi")
        with pytest.raises(OSError):
            self.toolkit_wrapper.to_file(ETHANOL, filename, format_name)


@pytest.mark.usefixtures("init_toolkit", "tmpdir")
@requires_openeye
class TestOpenEyeToolkitToFileIO(BaseToFileIO):
    toolkit_wrapper_class = OpenEyeToolkitWrapper


@pytest.mark.usefixtures("init_toolkit", "tmpdir")
@requires_rdkit
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
            match=(
                "'hydrogens_are_explicit' was specified as True, but (OpenEye|RDKit) [Tt]oolkit interpreted SMILES "
                "[^ ]+ as having implicit hydrogen. If this SMILES is intended to express all explicit hydrogens "
                f"in the molecule, then you should construct the desired molecule as an {self.tk_mol_name}"
            ),
        ):
            mol = self.toolkit_wrapper.from_smiles("C", hydrogens_are_explicit=True)

    def test_parse_methane_with_explicit_Hs(self):
        mol = self.toolkit_wrapper.from_smiles("[C]([H])([H])([H])([H])")
        # add hydrogens
        assert mol.n_atoms == 5
        assert mol.n_bonds == 4
        assert mol.partial_charges is None

    def test_parse_methane_with_explicit_Hs_and_say_they_are_explicit(self):
        mol = self.toolkit_wrapper.from_smiles(
            "[C]([H])([H])([H])([H])", hydrogens_are_explicit=True
        )
        # add hydrogens
        assert mol.n_atoms == 5
        assert mol.n_bonds == 4

    def test_parse_bad_smiles(self):
        with pytest.raises(ValueError, match="Unable to parse the SMILES string"):
            mol = self.toolkit_wrapper.from_smiles("QWERT")

    ### Copied from test_toolkits.py

    @pytest.mark.parametrize(
        "title, smiles",
        [
            ("unspec_chiral_smiles", r"C\C(F)=C(/F)CC(C)(Cl)Br"),
            ("spec_chiral_smiles", r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"),
            ("unspec_db_smiles", r"CC(F)=C(F)C[C@@](C)(Cl)Br"),
            ("spec_db_smiles", r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"),
        ],
    )
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
@requires_openeye
class TestOpenEyeToolkitSmiles(BaseSmiles):
    toolkit_wrapper_class = OpenEyeToolkitWrapper
    tk_mol_name = "OEMol"


@pytest.mark.usefixtures("init_toolkit")
@requires_rdkit
class TestRDKitToolkitSmiles(BaseSmiles):
    toolkit_wrapper_class = RDKitToolkitWrapper
    tk_mol_name = "RDMol"


# ========================================================
# Base class to test checks for unsupported chemsitry
# ========================================================


def test_unsupported_chemistry_exception_hierarchy():
    assert issubclass(
        exceptions.DisconnectedMoleculesError, exceptions.UnsupportedChemistryError
    )
    assert issubclass(
        exceptions.UnsupportedAtomTypeError, exceptions.UnsupportedChemistryError
    )
    assert issubclass(
        exceptions.UnsupportedBondTypeError, exceptions.UnsupportedChemistryError
    )


class BaseUnsupportedChemistry:
    def test_smiles_disconnected(self):
        with pytest.raises(
            exceptions.DisconnectedMoleculesError,
            match="OpenFF does not currently support input structures with more than one disconnected molecule",
        ):
            self.toolkit_wrapper.from_smiles("[Na+].[Cl-]")

    def test_smiles_connected_but_with_dot_disconnect(self):
        mol = self.toolkit_wrapper.from_smiles("C1.C1")
        assert mol.n_atoms == 8

    def test_smiles_with_wildcard(self):
        with pytest.raises(
            exceptions.UnsupportedAtomTypeError,
            match="OpenFF does not support atoms with an atomic number of 0",
        ):
            mol = self.toolkit_wrapper.from_smiles("*C")

    def test_inchi_disconnected(self):
        with pytest.raises(
            exceptions.DisconnectedMoleculesError,
            match="OpenFF does not currently support input structures with more than one disconnected molecule",
        ):
            # [Na+].[Cl-]
            self.toolkit_wrapper.from_inchi("InChI=1S/ClH.Na/h1H;/q;+1/p-1")
            # InChI doesn't accept atoms with element number 0 so can't test for that case.

    def test_sdf_disconnected(self):
        with pytest.raises(
            exceptions.DisconnectedMoleculesError,
            match="OpenFF does not currently support input structures with more than one disconnected molecule",
        ):
            self.toolkit_wrapper.from_file_obj(file_obj_manager.table_salt_sdf, "sdf")

    def test_sdf_with_R_group(self):
        with pytest.raises(
            exceptions.UnsupportedAtomTypeError,
            match="OpenFF does not support atoms with an atomic number of 0",
        ):
            self.toolkit_wrapper.from_file_obj(file_obj_manager.methyl_group_sdf, "sdf")


@pytest.mark.usefixtures("init_toolkit")
@requires_openeye
class TestOpenEyeToolkitUnsupportedChemistry(BaseUnsupportedChemistry):
    toolkit_wrapper_class = OpenEyeToolkitWrapper
    tk_mol_name = "OEMol"

    def test_empty_rdmol_is_accepted(self):
        # This verifies the for the number of components is ">1" not "!=1".
        from openeye import oechem

        rdmol = oechem.OEGraphMol()
        mol = self.toolkit_wrapper.from_object(rdmol)
        assert mol.n_atoms == 0

    ## def test_sdf_with_any_bond(self):
    ##     raise AssertionError("Cannot test - OEChem handles it without a problem.)


@pytest.mark.usefixtures("init_toolkit")
@requires_rdkit
class TestRDKitToolkitUnsupportedChemistry(BaseUnsupportedChemistry):
    toolkit_wrapper_class = RDKitToolkitWrapper
    tk_mol_name = "RDMol"

    def test_empty_rdmol_is_accepted(self):
        # This verifies the for the number of components is ">1" not "!=1".
        from rdkit import Chem

        rdmol = Chem.Mol()
        mol = self.toolkit_wrapper.from_object(rdmol)
        assert mol.n_atoms == 0

    def test_sdf_with_any_bond(self):
        with pytest.raises(
            exceptions.UnsupportedBondTypeError,
            match="OpenFF does not currently support RDKit molecules with a bond of type UNSPECIFIED",
        ):
            self.toolkit_wrapper.from_file_obj(file_obj_manager.chebi_52729_sdf, "sdf")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__] + sys.argv[1:]))
