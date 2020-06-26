#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for forcefield class

"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy
import os

from simtk import openmm, unit
import numpy as np
from numpy.testing import assert_almost_equal

import pytest
from tempfile import NamedTemporaryFile

from openforcefield.utils.toolkits import (OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper,
    ToolkitRegistry, ChargeMethodUnavailableError)
from openforcefield.utils import get_data_file_path
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import (ForceField, IncompatibleParameterError, SMIRNOFFSpecError,
    XMLParameterIOHandler, ParameterHandler)


#======================================================================
# GLOBAL CONSTANTS
#======================================================================

# File paths.
TIP3P_SDF_FILE_PATH = get_data_file_path(os.path.join('systems', 'monomers', 'water.sdf'))

XML_FF_GENERICS = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Bonds version="0.3">
    <Bond smirks="[*:1]~[*:2]" id="b1" k="680.0 * kilocalories_per_mole/angstrom**2" length="1.09 * angstrom"/>
  </Bonds>
  <Angles version="0.3">
    <Angle smirks="[*:1]~[*:2]~[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]~[*:2]~[*:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[*:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1"  switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[*:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
"""

simple_xml_ff = str.encode('''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Bonds version="0.3">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0 * kilocalories_per_mole/angstrom**2" length="1.526 * angstrom"/>
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0 * kilocalories_per_mole/angstrom**2" length="1.51 * angstrom"/>
  </Bonds>
  <Angles version="0.3">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5 * degree" id="a2" k="70.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree" periodicity2="2" phase2="180.0 * degree" idivf2="1" k2="0.250 * kilocalories_per_mole" periodicity3="1" phase3="180.0 * degree" idivf3="1" k3="0.200 * kilocalories_per_mole"/>
    <Proper smirks="[*:1]:[#6X4:2]~[#6X4:3]:[*:4]" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157 * kilocalories_per_mole" id="n2" rmin_half="1.4870 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
''')

xml_ff_w_comments = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <!-- SMIRNOFF (SMIRKS Native Open Force Field) template file -->
  <Date>2018-07-14</Date>
  <Author>C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, UC Irvine; D. L. Mobley, UC Irvine</Author>
  <!-- This file is meant for processing via openforcefield.typing.engines.smirnoff -->
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the bond energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Bonds version="0.3">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0 * kilocalories_per_mole/angstrom**2" length="1.526 * angstrom" />
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0 * kilocalories_per_mole/angstrom**2" length="1.51 * angstrom"/>
  </Bonds>
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the angle energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Angles version="0.3">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5 * degree" id="a2" k="70.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree" periodicity2="2" phase2="180.0 * degree" idivf2="1" k2="0.250 * kilocalories_per_mole" periodicity3="1" phase3="180.0 * degree" idivf3="1" k3="0.200 * kilocalories_per_mole"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157 * kilocalories_per_mole" id="n2" rmin_half="1.4870 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" scale15="1" cutoff="9.0 * angstrom" pme_tolerance="0.00001"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
'''


xml_ff_w_cosmetic_elements = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <!-- SMIRNOFF (SMIRKS Native Open Force Field) template file -->
  <Date>MMXVIII-VII-XIV</Date>
  <Author>Alice and Bob</Author>
  <!-- This file is meant for processing via openforcefield.typing.engines.smirnoff -->
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the bond energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Bonds version="0.3">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0 * kilocalories_per_mole/angstrom**2" length="1.526 * angstrom" parameters="k, length" parameterize_eval="blah=blah2"/>
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0 * kilocalories_per_mole/angstrom**2" length="1.51 * angstrom"/>
  </Bonds>
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the angle energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Angles version="0.3" cosmetic_element="why not?">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5 * degree" id="a2" k="70.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree" periodicity2="2" phase2="180.0 * degree" idivf2="1" k2="0.250 * kilocalories_per_mole" periodicity3="1" phase3="180.0 * degree" idivf3="1" k3="0.200 * kilocalories_per_mole"/>
    <Proper smirks="[*:1]:[#6X4:2]~[#6X4:3]:[*:4]" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157 * kilocalories_per_mole" id="n2" rmin_half="1.4870 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" scale15="1" cutoff="9.0 * angstrom" pme_tolerance="0.00001"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
'''

xml_toolkitam1bcc_ff = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
'''

xml_ethanol_library_charges_ff = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]-[#6:2](-[#1:3])(-[#1:4])-[#6:5](-[#1:6])(-[#1:7])-[#8:8]-[#1:9]" charge1="-0.02*elementary_charge" charge2="-0.2*elementary_charge" charge3="-0.02*elementary_charge" charge4="-0.02*elementary_charge" charge5="-0.1*elementary_charge" charge6="-0.01*elementary_charge" charge7="-0.01*elementary_charge" charge8="0.3*elementary_charge" charge9="0.08*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
'''

xml_ethanol_library_charges_in_parts_ff = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <!-- Note that the oxygen is covered twice here. The correct behavior should be to take the charge from the SECOND LibraryCharge, as it should overwrite the first -->
       <LibraryCharge smirks="[#1:1]-[#6:2](-[#1:3])(-[#1:4])-[#6:5](-[#1:6])(-[#1:7])-[#8:8]" charge1="-0.02*elementary_charge" charge2="-0.2*elementary_charge" charge3="-0.02*elementary_charge" charge4="-0.02*elementary_charge" charge5="-0.1*elementary_charge" charge6="-0.01*elementary_charge" charge7="-0.01*elementary_charge" charge8="-999*elementary_charge" />
       <LibraryCharge smirks="[#8:1]-[#1:2]" charge1="0.3*elementary_charge" charge2="0.08*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
'''

xml_ethanol_library_charges_by_atom_ff = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]-[#6]" charge1="-0.02*elementary_charge" />
       <LibraryCharge smirks="[#6X4:1]" charge1="-0.2*elementary_charge" />
       <LibraryCharge smirks="[#1:1]-[#6]-[#8]" charge1="-0.01*elementary_charge" />
       <LibraryCharge smirks="[#6X4:1]-[#8]" charge1="-0.1*elementary_charge" />
       <LibraryCharge smirks="[#8X2:1]" charge1="0.3*elementary_charge" />
       <LibraryCharge smirks="[#1:1]-[#8]" charge1="0.08*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
'''

xml_OH_library_charges_xml = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]" charge1="1.*elementary_charge" />
       <LibraryCharge smirks="[#8:1]" charge1="-2.*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
'''

xml_CH_zeroes_library_charges_xml = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]" charge1="0.*elementary_charge" />
       <LibraryCharge smirks="[#6:1]" charge1="0.*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
'''

xml_spec_docs_ala_library_charges_xml = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge name="ALA" smirks="[NX3:1]([#1:2])([#6])[#6H1:3]([#1:4])([#6:5]([#1:6])([#1:7])[#1:8])[#6:9](=[#8:10])[#7]" charge1="-0.4157*elementary_charge" charge2="0.2719*elementary_charge" charge3="0.0337*elementary_charge" charge4="0.0823*elementary_charge" charge5="-0.1825*elementary_charge" charge6="0.0603*elementary_charge" charge7="0.0603*elementary_charge" charge8="0.0603*elementary_charge" charge9="0.5973*elementary_charge" charge10="-0.5679*elementary_charge"/>
    </LibraryCharges>
</SMIRNOFF>
'''

xml_spec_docs_tip3p_library_charges_xml = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge name="TIP3P" smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]" charge1="0.417*elementary_charge" charge2="-0.834*elementary_charge" charge3="0.417*elementary_charge"/>
    </LibraryCharges>
</SMIRNOFF>
'''

xml_spec_docs_charge_increment_model_xml = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="AM1-Mulliken">
    <!-- A fractional charge can be moved along a single bond -->
    <ChargeIncrement smirks="[#6X4:1]-[#6X3a:2]" charge_increment1="-0.0073*elementary_charge" charge_increment2="0.0073*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1]-[#6X3a:2]-[#7]" charge_increment1="0.0943*elementary_charge" charge_increment2="-0.0943*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge_increment1="-0.0718*elementary_charge" charge_increment2="0.0718*elementary_charge"/>
    <!--- Alternatively, fractional charges can be redistributed among any number of bonded atoms -->
    <ChargeIncrement smirks="[N:1]([H:2])([H:3])" charge_increment1="0.02*elementary_charge" charge_increment2="-0.01*elementary_charge" charge_increment3="-0.01*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>
'''

xml_charge_increment_model_formal_charges = '''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ChargeIncrementModel version="0.3" number_of_conformers="0" partial_charge_method="formal_charge"/>
</SMIRNOFF>
'''

xml_ff_torsion_bo = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]~[#6X3:2]~[#6X3:3]~[*:4]" id="tbo1" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
    <Proper smirks="[*:1]~[#6X4:2]~[#8X2:3]~[*:4]" id="tbo2" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
</SMIRNOFF>
'''

xml_ff_torsion_bo_standard_supersede = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]~[#6X3:2]~[#6X3:3]~[*:4]" id="tbo1" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
    <Proper smirks="[*:1]~[#6X4:2]~[#8X2:3]~[*:4]" id="tbo2" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
    <Proper smirks="[*:1]~[#6X4:2]-[#8X2:3]~[#1:4]" id="t1" periodicity1="2" phase1="0.0 * degree" k1="1.20*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
</SMIRNOFF>
'''
#======================================================================
# TEST UTILITY FUNCTIONS
#======================================================================

def round_charge(xml):
    """Round charge fields in a serialized OpenMM system to 2 decimal places"""
    # Example Particle line: 				<Particle eps=".4577296" q="-.09709000587463379" sig=".1908"/>
    xmlsp = xml.split(' q="')
    for index, chunk in enumerate(xmlsp):
        # Skip file before first q=
        if index == 0:
            continue
        chunksp = chunk.split('" sig')
        chunksp[0] = f"{float(chunksp[0]):.2f}"
        chunk = '" sig'.join(chunksp)
        xmlsp[index] = chunk
    return ' q="'.join(xmlsp)


def create_cis_1_2_dichloroethene():
    """
    Creates an openforcefield.topology.Molecule representation of cis-1,2-dichloroethene
    without the use of a cheminformatics toolkit.
    """

    cis_dichloroethene = Molecule()
    cis_dichloroethene.add_atom(17, 0, False)
    cis_dichloroethene.add_atom(6, 0, False)
    cis_dichloroethene.add_atom(6, 0, False)
    cis_dichloroethene.add_atom(17, 0, False)
    cis_dichloroethene.add_atom(1, 0, False)
    cis_dichloroethene.add_atom(1, 0, False)
    cis_dichloroethene.add_bond(0, 1, 1, False)
    cis_dichloroethene.add_bond(1, 2, 2, False, 'Z')
    cis_dichloroethene.add_bond(2, 3, 1, False)
    cis_dichloroethene.add_bond(1, 4, 1, False)
    cis_dichloroethene.add_bond(2, 5, 1, False)
    return cis_dichloroethene


def create_ethanol():
    """
    Creates an openforcefield.topology.Molecule representation of
    ethanol without the use of a cheminformatics toolkit
    """
    # Create an ethanol molecule without using a toolkit
    ethanol = Molecule()
    ethanol.add_atom(6, 0, False)  # C0
    ethanol.add_atom(6, 0, False)  # C1
    ethanol.add_atom(8, 0, False)  # O2
    ethanol.add_atom(1, 0, False)  # H3
    ethanol.add_atom(1, 0, False)  # H4
    ethanol.add_atom(1, 0, False)  # H5
    ethanol.add_atom(1, 0, False)  # H6
    ethanol.add_atom(1, 0, False)  # H7
    ethanol.add_atom(1, 0, False)  # H8
    ethanol.add_bond(0, 1, 1, False, fractional_bond_order=1.33)  # C0 - C1
    ethanol.add_bond(1, 2, 1, False, fractional_bond_order=1.23)  # C1 - O2
    ethanol.add_bond(0, 3, 1, False, fractional_bond_order=1)  # C0 - H3
    ethanol.add_bond(0, 4, 1, False, fractional_bond_order=1)  # C0 - H4
    ethanol.add_bond(0, 5, 1, False, fractional_bond_order=1)  # C0 - H5
    ethanol.add_bond(1, 6, 1, False, fractional_bond_order=1)  # C1 - H6
    ethanol.add_bond(1, 7, 1, False, fractional_bond_order=1)  # C1 - H7
    ethanol.add_bond(2, 8, 1, False, fractional_bond_order=1)  # O2 - H8
    charges = unit.Quantity(np.array([-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4]), unit.elementary_charge)
    ethanol.partial_charges = charges

    return ethanol


def create_reversed_ethanol():
    """
    Creates an openforcefield.topology.Molecule representation of
    ethanol without the use of a cheminformatics toolkit. This function
    reverses the atom indexing of create_ethanol
    """
    # Create an ethanol molecule without using a toolkit
    ethanol = Molecule()
    ethanol.add_atom(1, 0, False)  # H0
    ethanol.add_atom(1, 0, False)  # H1
    ethanol.add_atom(1, 0, False)  # H2
    ethanol.add_atom(1, 0, False)  # H3
    ethanol.add_atom(1, 0, False)  # H4
    ethanol.add_atom(1, 0, False)  # H5
    ethanol.add_atom(8, 0, False)  # O6
    ethanol.add_atom(6, 0, False)  # C7
    ethanol.add_atom(6, 0, False)  # C8
    ethanol.add_bond(8, 7, 1, False, fractional_bond_order=1.33)  # C8 - C7
    ethanol.add_bond(7, 6, 1, False, fractional_bond_order=1.23)  # C7 - O6
    ethanol.add_bond(8, 5, 1, False, fractional_bond_order=1)  # C8 - H5
    ethanol.add_bond(8, 4, 1, False, fractional_bond_order=1)  # C8 - H4
    ethanol.add_bond(8, 3, 1, False, fractional_bond_order=1)  # C8 - H3
    ethanol.add_bond(7, 2, 1, False, fractional_bond_order=1)  # C7 - H2
    ethanol.add_bond(7, 1, 1, False, fractional_bond_order=1)  # C7 - H1
    ethanol.add_bond(6, 0, 1, False, fractional_bond_order=1)  # O6 - H0
    charges = unit.Quantity(np.array([0.4, 0.3, 0.2, 0.1, 0.00001, -0.1, -0.2, -0.3, -0.4]), unit.elementary_charge)
    ethanol.partial_charges = charges
    return ethanol

def create_benzene_no_aromatic():
    """
    Creates an openforcefield.topology.Molecule representation of benzene through the API with aromatic bonds
    not defied, used to test the levels of isomorphic matching.
    """
    benzene = Molecule()
    benzene.add_atom(6, 0, False)  # C0
    benzene.add_atom(6, 0, False)  # C1
    benzene.add_atom(6, 0, False)  # C2
    benzene.add_atom(6, 0, False)  # C3
    benzene.add_atom(6, 0, False)  # C4
    benzene.add_atom(6, 0, False)  # C5
    benzene.add_atom(1, 0, False)  # H6
    benzene.add_atom(1, 0, False)  # H7
    benzene.add_atom(1, 0, False)  # H8
    benzene.add_atom(1, 0, False)  # H9
    benzene.add_atom(1, 0, False)  # H10
    benzene.add_atom(1, 0, False)  # H11
    benzene.add_bond(0, 5, 1, False)  # C0 - C5
    benzene.add_bond(0, 1, 1, False)  # C0 - C1
    benzene.add_bond(1, 2, 1, False)  # C1 - C2
    benzene.add_bond(2, 3, 1, False)  # C2 - C3
    benzene.add_bond(3, 4, 1, False)  # C3 - C4
    benzene.add_bond(4, 5, 1, False)  # C4 - C5
    benzene.add_bond(0, 6, 1, False)  # C0 - H6
    benzene.add_bond(1, 7, 1, False)  # C1 - C7
    benzene.add_bond(2, 8, 1, False)  # C2 - C8
    benzene.add_bond(3, 9, 1, False)  # C3 - C9
    benzene.add_bond(4, 10, 1, False)  # C4 - C10
    benzene.add_bond(5, 11, 1, False)  # C5 - C11
    return benzene

def create_acetaldehyde():
    """
    Creates an openforcefield.topology.Molecule representation of acetaldehyde through the API
    """
    acetaldehyde = Molecule()
    acetaldehyde.add_atom(6, 0, False)  # C0
    acetaldehyde.add_atom(6, 0, False)  # C1
    acetaldehyde.add_atom(8, 0, False)  # O2
    acetaldehyde.add_atom(1, 0, False)  # H3
    acetaldehyde.add_atom(1, 0, False)  # H4
    acetaldehyde.add_atom(1, 0, False)  # H5
    acetaldehyde.add_atom(1, 0, False)  # H6
    acetaldehyde.add_bond(0, 1, 1, False)  # C0 - C1
    acetaldehyde.add_bond(1, 2, 2, False)  # C1 = O2
    acetaldehyde.add_bond(0, 3, 1, False)  # C0 - H3
    acetaldehyde.add_bond(0, 4, 1, False)  # C0 - H4
    acetaldehyde.add_bond(0, 5, 1, False)  # C0 - H5
    acetaldehyde.add_bond(1, 6, 1, False)  # C1 - H6
    charges = unit.Quantity(np.array([0, 0, 0, 0, 0, 0, 0]), unit.elementary_charge)
    acetaldehyde.partial_charges = charges
    return acetaldehyde

def create_acetate():
    """
    Creates an openforcefield.topology.Molecule representation of
    acetate without the use of a cheminformatics toolkit
    """
    # Create an acetate molecule without using a toolkit
    acetate = Molecule()
    acetate.add_atom(6, 0, False)  # C0
    acetate.add_atom(6, 0, False)  # C1
    acetate.add_atom(8, 0, False)  # O2
    acetate.add_atom(8, -1, False) # O3
    acetate.add_atom(1, 0, False)  # H4
    acetate.add_atom(1, 0, False)  # H5
    acetate.add_atom(1, 0, False)  # H6
    acetate.add_bond(0, 1, 1, False)  # C0 - C1
    acetate.add_bond(1, 2, 2, False)  # C1 = O2
    acetate.add_bond(1, 3, 1, False)  # C1 - O3[-1]
    acetate.add_bond(0, 4, 1, False)  # C0 - H4
    acetate.add_bond(0, 5, 1, False)  # C0 - H5
    acetate.add_bond(0, 6, 1, False)  # C0 - H6
    return acetate

def create_cyclohexane():
    """
    Creates an openforcefield.topology.Molecule representation of
    cyclohexane without the use of a cheminformatics toolkit
    """
    cyclohexane = Molecule()
    cyclohexane.add_atom(6, 0, False)  # C0
    cyclohexane.add_atom(6, 0, False)  # C1
    cyclohexane.add_atom(6, 0, False)  # C2
    cyclohexane.add_atom(6, 0, False)  # C3
    cyclohexane.add_atom(6, 0, False)  # C4
    cyclohexane.add_atom(6, 0, False)  # C5
    cyclohexane.add_atom(1, 0, False)  # H6
    cyclohexane.add_atom(1, 0, False)  # H7
    cyclohexane.add_atom(1, 0, False)  # H8
    cyclohexane.add_atom(1, 0, False)  # H9
    cyclohexane.add_atom(1, 0, False)  # H10
    cyclohexane.add_atom(1, 0, False)  # H11
    cyclohexane.add_atom(1, 0, False)  # H12
    cyclohexane.add_atom(1, 0, False)  # H13
    cyclohexane.add_atom(1, 0, False)  # H14
    cyclohexane.add_atom(1, 0, False)  # H15
    cyclohexane.add_atom(1, 0, False)  # H16
    cyclohexane.add_atom(1, 0, False)  # H17
    cyclohexane.add_bond(0, 1, 1, False)  # C0 - C1
    cyclohexane.add_bond(1, 2, 1, False)  # C1 - C2
    cyclohexane.add_bond(2, 3, 1, False)  # C2 - C3
    cyclohexane.add_bond(3, 4, 1, False)  # C3 - C4
    cyclohexane.add_bond(4, 5, 1, False)  # C4 - C5
    cyclohexane.add_bond(5, 0, 1, False)  # C5 - C0
    cyclohexane.add_bond(0, 6, 1, False)  # C0 - H6
    cyclohexane.add_bond(0, 7, 1, False)  # C0 - H7
    cyclohexane.add_bond(1, 8, 1, False)  # C1 - H8
    cyclohexane.add_bond(1, 9, 1, False)  # C1 - H9
    cyclohexane.add_bond(2, 10, 1, False)  # C2 - H10
    cyclohexane.add_bond(2, 11, 1, False)  # C2 - H11
    cyclohexane.add_bond(3, 12, 1, False)  # C3 - H12
    cyclohexane.add_bond(3, 13, 1, False)  # C3 - H13
    cyclohexane.add_bond(4, 14, 1, False)  # C4 - H14
    cyclohexane.add_bond(4, 15, 1, False)  # C4 - H15
    cyclohexane.add_bond(5, 16, 1, False)  # C5 - H16
    cyclohexane.add_bond(5, 17, 1, False)  # C5 - H17
    return cyclohexane



nonbonded_resolution_matrix = [
    {'vdw_method': 'cutoff', 'electrostatics_method': 'Coulomb', 'has_periodic_box': True,
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'Coulomb', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'reaction-field', 'has_periodic_box': True,
     'omm_force': None, 'exception': SMIRNOFFSpecError, 'exception_match': 'reaction-field'},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'reaction-field', 'has_periodic_box': False,
     'omm_force': None, 'exception': SMIRNOFFSpecError, 'exception_match': 'reaction-field'},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'PME', 'has_periodic_box': True,
     'omm_force': openmm.NonbondedForce.PME, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'PME', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},

    {'vdw_method': 'PME', 'electrostatics_method': 'Coulomb', 'has_periodic_box': True,
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'Coulomb', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'reaction-field', 'has_periodic_box': True,
     'omm_force': None, 'exception': SMIRNOFFSpecError, 'exception_match': 'reaction-field'},
    {'vdw_method': 'PME', 'electrostatics_method': 'reaction-field', 'has_periodic_box': False,
     'omm_force': None, 'exception': SMIRNOFFSpecError, 'exception_match': 'reaction-field'},
    {'vdw_method': 'PME', 'electrostatics_method': 'PME', 'has_periodic_box': True,
     'omm_force': openmm.NonbondedForce.LJPME, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'PME', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},
     ]

partial_charge_method_resolution_matrix = [
    {'toolkit': AmberToolsToolkitWrapper,
     'partial_charge_method': 'AM1-Mulliken',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': AmberToolsToolkitWrapper,
     'partial_charge_method': 'Gasteiger',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': AmberToolsToolkitWrapper,
     'partial_charge_method': 'Madeup-ChargeMethod',
     'exception': ChargeMethodUnavailableError,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'AM1-Mulliken',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'Gasteiger',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'MMFF94',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'am1bcc',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'am1bccnosymspt',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'am1bccelf10',
     'exception': None,
     'exception_match': ''
     },
    {'toolkit': OpenEyeToolkitWrapper,
     'partial_charge_method': 'Madeup-ChargeMethod',
     'exception': ChargeMethodUnavailableError,
     'exception_match': ''
     }
     ]


#=============================================================================================
# TESTS
#=============================================================================================


toolkit_registries = []
if OpenEyeToolkitWrapper.is_available():
    toolkit_registries.append((ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper]), "OE"))
if RDKitToolkitWrapper.is_available() and AmberToolsToolkitWrapper.is_available():
    toolkit_registries.append((ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper]),
                               'RDKit+AmberTools'))


class TestForceField():
    """Test the ForceField class"""

    def test_create_forcefield_no_args(self):
        """Test empty constructor"""
        forcefield = ForceField()

        # Should find BondHandler and AngleHandler, since they're default classes
        forcefield.get_parameter_handler('Bonds')
        forcefield.get_parameter_handler('Angles')

        # Shouldn't find InvalidKey handler, since it doesn't exist
        with pytest.raises(KeyError) as excinfo:
            forcefield.get_parameter_handler('InvalidKey')

    def test_create_forcefield_custom_handler_classes(self):
        """Test constructor given specific classes to register"""
        from openforcefield.typing.engines.smirnoff import BondHandler
        forcefield = ForceField(parameter_handler_classes=[BondHandler])

        # Should find BondHandler, since we registered it
        forcefield.get_parameter_handler('Bonds')

        # Shouldn't find AngleHandler, since we didn't allow that to be registered
        with pytest.raises(KeyError) as excinfo:
            forcefield.get_parameter_handler('Angles')

    def test_create_forcefield_from_file(self):
        """Test basic file loading in constructor"""
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        assert len(forcefield._parameter_handlers['Bonds']._parameters) == 87
        assert len(forcefield._parameter_handlers['Angles']._parameters) == 38
        assert len(forcefield._parameter_handlers['ProperTorsions']._parameters) == 158
        assert len(forcefield._parameter_handlers['ImproperTorsions']._parameters) == 4
        assert len(forcefield._parameter_handlers['vdW']._parameters) == 35

    def test_load_bad_string(self):
        with pytest.raises(IOError) as exception_info:
            ForceField('1234')
        assert 'Source 1234 could not be read.' in str(exception_info.value)
        assert 'syntax error' in str(exception_info.value)

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_file_list(self):
        # These offxml files are located in package data path, which is automatically installed and searched
        file_paths = [smirnoff99Frosst_offxml_file_path, tip3p_offxml_file_path]
        # Create a forcefield from multiple offxml files
        forcefield = ForceField(file_paths)

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_file_path_iterator(self):
        # These offxml files are located in package data path, which is automatically installed and searched
        file_paths = [smirnoff99Frosst_offxml_file_path, tip3p_offxml_file_path]
        # A generator should work as well
        forcefield = ForceField(iter(file_paths))

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_gbsa():
        """Test reading of ffxml files with GBSA support.
        """
        forcefield = ForceField('test_forcefields/Frosst_AlkEthOH_GBSA.offxml')

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_url(self):
        urls = [
            'https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/test_forcefields/smirnoff99Frosst.offxml',
            'https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/test_forcefields/tip3p.offxml'
        ]
        # Test creation with smirnoff99frosst URL
        forcefield = ForceField(urls[0])

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_url_list(self):
        urls = [
            'https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/test_forcefields/smirnoff99Frosst.offxml',
            'https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/test_forcefields/tip3p.offxml'
        ]
        # Test creation with multiple URLs
        forcefield = ForceField(urls)

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_url_iterator(self):
        urls = [
            'https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/test_forcefields/smirnoff99Frosst.offxml',
            'https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/test_forcefields/tip3p.offxml'
        ]
        # A generator should work as well
        forcefield = ForceField(iter(urls))


    def test_create_forcefield_from_xml_string(self):
        forcefield = ForceField(simple_xml_ff)
        assert len(forcefield._parameter_handlers['Bonds']._parameters) == 2
        assert len(forcefield._parameter_handlers['Angles']._parameters) == 2
        assert len(forcefield._parameter_handlers['ProperTorsions']._parameters) == 3
        assert len(forcefield._parameter_handlers['ImproperTorsions']._parameters) == 2
        assert len(forcefield._parameter_handlers['vdW']._parameters) == 2

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_deep_copy(self):
        forcefield = ForceField(smirnoff99Frosst_offxml_file_path)
        # Deep copy
        forcefield2 = copy.deepcopy(cls.forcefield)
        assert_forcefields_equal(cls.forcefield, forcefield2,
                                 "ForceField deep copy does not match original ForceField")


    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    # TODO: This should check the output of forcefield.to_dict
    def test_serialize(self):

        forcefield = ForceField(smirnoff99Frosst_offxml_file_path)
        # Serialize/deserialize
        serialized_forcefield = cls.forcefield.__getstate__()
        forcefield2 = ForceField.__setstate__(serialized_forcefield)
        assert_forcefields_equal(cls.forcefield, forcefield2,
                                 "Deserialized serialized ForceField does not match original ForceField")

    def test_pickle(self):
        """
        Test pickling and unpickling a forcefield
        """
        import pickle
        forcefield_1 = ForceField(simple_xml_ff)
        pickled = pickle.dumps(forcefield_1)
        forcefield_2 = pickle.loads(pickled)
        assert forcefield_1.to_string() == forcefield_2.to_string()

    def test_pickle_with_cosmetic_attributes(self):
        """
        Test pickling and unpickling a forcefield with cosmetic attributes
        """
        import pickle
        forcefield_1 = ForceField(xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)
        pickled = pickle.dumps(forcefield_1)
        forcefield_2 = pickle.loads(pickled)
        assert forcefield_1.to_string() == forcefield_2.to_string()
        # Ensure that the cosmetic attributes stuck around
        assert 'blah=blah2' in forcefield_2.to_string()

    def test_xml_string_roundtrip(self):
        """
        Test writing a ForceField to an XML string
        """
        forcefield_1 = ForceField(simple_xml_ff)
        string_1 = forcefield_1.to_string('XML')
        forcefield_2 = ForceField(string_1)
        string_2 = forcefield_2.to_string('XML')
        assert string_1 == string_2


    def test_xml_string_roundtrip_keep_cosmetic(self):
        """
        Test roundtripping a forcefield to an XML string with and without retaining cosmetic elements
        """
        # Ensure an exception is raised if we try to read the XML string with cosmetic attributes
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg [(]parameters: k, length[)]  passed") as excinfo:
            forcefield = ForceField(xml_ff_w_cosmetic_elements)

        # Create a forcefield from XML successfully, by explicitly permitting cosmetic attributes
        forcefield_1 = ForceField(xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)

        # Convert the forcefield back to XML
        string_1 = forcefield_1.to_string('XML', discard_cosmetic_attributes=False)

        # Ensure that the new XML string has cosmetic attributes in it
        assert 'cosmetic_element="why not?"' in string_1
        assert 'parameterize_eval="blah=blah2"' in string_1
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg [(]parameters: k, length[)]  passed") as excinfo:
            forcefield = ForceField(string_1, allow_cosmetic_attributes=False)

        # Complete the forcefield_1 --> string --> forcefield_2 roundtrip
        forcefield_2 = ForceField(string_1, allow_cosmetic_attributes=True)

        # Ensure that the forcefield remains the same after the roundtrip
        string_2 = forcefield_2.to_string('XML', discard_cosmetic_attributes=False)
        assert string_1 == string_2

        # Discard the cosmetic attributes and ensure that the string is different
        string_3 = forcefield_2.to_string('XML', discard_cosmetic_attributes=True)
        assert string_1 != string_3
        # Ensure that the new XML string does NOT have cosmetic attributes in it
        assert 'cosmetic_element="why not?"' not in string_3
        assert 'parameterize_eval="blah=blah2"' not in string_3


    def test_read_0_1_smirnoff(self):
        """Test reading an 0.1 spec OFFXML file"""
        ff = ForceField('test_forcefields/smirnoff99Frosst_reference_0_1_spec.offxml')


    def test_read_0_1_smirff(self):
        """Test reading an 0.1 spec OFFXML file, enclosed by the legacy "SMIRFF" tag"""
        ff = ForceField('test_forcefields/smirff99Frosst_reference_0_1_spec.offxml')


    def test_read_0_2_smirnoff(self):
        """Test reading an 0.2 spec OFFXML file"""
        ff = ForceField('test_forcefields/smirnoff99Frosst_reference_0_2_spec.offxml')


    @pytest.mark.parametrize('file_path_extension', ['xml', 'XML', 'offxml', 'OFFXML'])
    @pytest.mark.parametrize('specified_format', [None, 'xml', 'XML', '.xml', '.XML',
                                                  'offxml', 'OFFXML', '.offxml', '.OFFXML',
                                                  XMLParameterIOHandler()])
    def test_xml_file_roundtrip(self, file_path_extension, specified_format):
        """
        Test roundtripping a ForceField to and from an XML file
        """
        # These files will be deleted once garbage collection runs (end of this function)
        iofile1 = NamedTemporaryFile(suffix='.' + file_path_extension)
        iofile2 = NamedTemporaryFile(suffix='.' + file_path_extension)
        forcefield_1 = ForceField(simple_xml_ff)
        forcefield_1.to_file(iofile1.name, io_format=specified_format)
        forcefield_2 = ForceField(iofile1.name)
        forcefield_2.to_file(iofile2.name, io_format=specified_format)
        assert open(iofile1.name).read() == open(iofile2.name).read()


    @pytest.mark.parametrize('file_path_extension', ['xml', 'XML', 'offxml', 'OFFXML'])
    @pytest.mark.parametrize('specified_format', [None, 'xml', 'XML', '.xml', '.XML',
                                                  'offxml', 'OFFXML', '.offxml', '.OFFXML',
                                                  XMLParameterIOHandler()])
    def test_xml_file_roundtrip_keep_cosmetic(self, file_path_extension, specified_format):
        """
        Test roundtripping a forcefield to an XML file with and without retaining cosmetic elements
        """
        # These files will be deleted once garbage collection runs (end of this function)
        iofile1 = NamedTemporaryFile(suffix='.' + file_path_extension)
        iofile2 = NamedTemporaryFile(suffix='.' + file_path_extension)
        iofile3 = NamedTemporaryFile(suffix='.' + file_path_extension)

        # Ensure an exception is raised if we try to read the XML string with cosmetic attributes
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg [(]parameters: k, length[)]  passed") as excinfo:
            forcefield = ForceField(xml_ff_w_cosmetic_elements)

        # Create a forcefield from XML successfully
        forcefield_1 = ForceField(xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)

        # Convert the forcefield back to XML, keeping cosmetic attributes
        forcefield_1.to_file(iofile1.name, discard_cosmetic_attributes=False, io_format=specified_format)

        # Ensure that the new XML string has cosmetic attributes in it
        assert 'cosmetic_element="why not?"' in open(iofile1.name).read()
        assert 'parameterize_eval="blah=blah2"' in open(iofile1.name).read()
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg [(]parameters: k, length[)]  passed") as excinfo:
            forcefield = ForceField(iofile1.name, allow_cosmetic_attributes=False)

        # Complete the forcefield_1 --> file --> forcefield_2 roundtrip
        forcefield_2 = ForceField(iofile1.name, allow_cosmetic_attributes=True)

        # Ensure that the forcefield remains the same after the roundtrip
        forcefield_2.to_file(iofile2.name, discard_cosmetic_attributes=False, io_format=specified_format)
        assert open(iofile1.name).read() == open(iofile2.name).read()

        # Discard the cosmetic attributes and ensure that the string is different
        forcefield_2.to_file(iofile3.name, discard_cosmetic_attributes=True, io_format=specified_format)
        assert open(iofile1.name).read() != open(iofile3.name).read()

        # Ensure that the new XML string does NOT have cosmetic attributes in it
        assert 'cosmetic_element="why not?"' not in open(iofile3.name).read()
        assert 'parameterize_eval="blah=blah2"' not in open(iofile3.name).read()

    def test_load_section_without_section_version(self):
        """Ensure that a SMIRNOFFSpecError is raised if we try to load a SMIRNOFF section without a version.
        Section versions are a requirement added in the 0.3 spec."""
        with pytest.raises(SMIRNOFFSpecError, match="Missing version while trying to construct "
                                                    "<class 'openforcefield.typing.engines."
                                                    "smirnoff.parameters.ToolkitAM1BCCHandler'>.") as excinfo:
            ff = ForceField('<?xml version="1.0" encoding="ASCII"?>'
                            '<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">'
                            '  <ToolkitAM1BCC/>'
                            '</SMIRNOFF>')


    def test_load_two_sources(self):
        """Test loading data from two SMIRNOFF data sources"""
        ff = ForceField(simple_xml_ff, xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)
        assert len(ff.get_parameter_handler('Bonds').parameters) == 4


    def test_load_two_sources_authors_dates(self):
        """Test that authors and dates are handled properly"""
        ff = ForceField(xml_ff_w_cosmetic_elements, xml_ff_w_comments, allow_cosmetic_attributes=True)
        xml_str = ff.to_string('XML')
        assert '<Author>Alice and Bob AND C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, ' \
               'UC Irvine; D. L. Mobley, UC Irvine</Author>' in xml_str
        assert '<Date>MMXVIII-VII-XIV AND 2018-07-14</Date>' in xml_str

        # Test property getters
        assert 'Alice and Bob AND C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, ' \
               'UC Irvine; D. L. Mobley, UC Irvine' == ff.author
        assert 'MMXVIII-VII-XIV AND 2018-07-14' == ff.date

        # Test property setters
        ff.author = 'Me'
        ff.date = 'yesteryear'
        xml_str = ff.to_string('XML')
        assert '<Author>Me</Author>' in xml_str
        assert '<Date>yesteryear</Date>' in xml_str

        # Unset both author and date and ensure they don't get written out.
        ff.author = None
        ff.date = None
        xml_str = ff.to_string('XML')
        assert '<Author>' not in xml_str
        assert '<Date>' not in xml_str


    def test_load_two_sources_incompatible_tags(self):
        """Test loading data from two SMIRNOFF data sources which have incompatible physics"""
        # Make an XML forcefield with a modifiedvdW 1-4 scaling factor
        nonstandard_xml_ff = xml_ff_w_comments.replace('scale14="0.5"', 'scale14="1.0"')
        with pytest.raises(IncompatibleParameterError, match="handler value: 0.5, incompatible value: 1.0") as excinfo:
            ff = ForceField(simple_xml_ff, nonstandard_xml_ff)


    def test_gbsahandler_sa_model_none(self):
        """
        Ensure that string values of "None" are correctly interpreted in the GBSAHandler's sa_model field
        """
        gbsa_ff_xml = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <GBSA version="0.3" gb_model="HCT" solvent_dielectric="78.5" solute_dielectric="1" sa_model="None" surface_area_penalty="5.4*calories/mole/angstroms**2" solvent_radius="1.4*angstroms">
          <Atom smirks="[*:1]" radius="0.15*nanometer" scale="0.8"/>
    </GBSA>
</SMIRNOFF>
'''
        from openforcefield.typing.engines.smirnoff import ForceField
        ff = ForceField(gbsa_ff_xml)


    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_ethanol(self, toolkit_registry, registry_description):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        molecules = [create_ethanol()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology, toolkit_registry=toolkit_registry)

    @pytest.fixture()
    def create_circular_handler_dependencies(self):
        from openforcefield.typing.engines.smirnoff.parameters import BondHandler, AngleHandler, ConstraintHandler

        # Modify the BondHandler and AngleHandler classes to depend on the other one running first during
        # system parameterization. Unfortunately, I can't figure out how to do this just to these _instances_
        # of the Handlers, so I modify them at the class level, and then un-modify them at the end of the test.
        orig_bh_depends = copy.deepcopy(BondHandler._DEPENDENCIES)
        orig_ah_depends = copy.deepcopy(AngleHandler._DEPENDENCIES)
        BondHandler._DEPENDENCIES = [ConstraintHandler, AngleHandler]
        AngleHandler._DEPENDENCIES = [ConstraintHandler, BondHandler]

        # The tests run here. Regardless of outcome, the code after `yield` runs after the test completes
        yield

        # Return handler dependencies to their original states
        BondHandler._DEPENDENCIES = orig_bh_depends
        AngleHandler._DEPENDENCIES = orig_ah_depends

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_ethanol_handler_dependency_loop(self, create_circular_handler_dependencies, toolkit_registry, registry_description):
        """Test parameterizing ethanol, but failing because custom handler classes can not resolve
         which order to run in"""
        from simtk.openmm import app
        # from openforcefield.typing.engines.smirnoff.parameters import BondHandler, AngleHandler, ConstraintHandler

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        molecules = [create_ethanol()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        with pytest.raises(RuntimeError, match="Unable to resolve order in which to run ParameterHandlers. "
                                               "Dependencies do not form a directed acyclic graph") as excinfo:
            omm_system = forcefield.create_openmm_system(topology, toolkit_registry=toolkit_registry)


    def test_parameterize_ethanol_missing_torsion(self):
        from simtk.openmm import app
        from openforcefield.typing.engines.smirnoff.parameters import UnassignedProperTorsionParameterException

        forcefield = ForceField('''
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[#99:1]-[#99X4:2]-[#99:3]-[#99:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
  </ProperTorsions>
</SMIRNOFF>
''')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        molecules = [create_ethanol()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        with pytest.raises(UnassignedProperTorsionParameterException,
                           match='- Topology indices [(]5, 0, 1, 6[)]: '
                                 'names and elements [(](H\d+)? H[)], [(](C\d+)? C[)], [(](C\d+)? C[)], [(](H\d+)? H[)],') \
                as excinfo:
            omm_system = forcefield.create_openmm_system(topology)

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_1_cyclohexane_1_ethanol(self, toolkit_registry, registry_description):
        """Test parameterizing a periodic system of two distinct molecules"""
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_cyclohexane_1_ethanol.pdb'))
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [create_ethanol(), create_cyclohexane()]
        # molecules = [Molecule.from_file(get_data_file_path(name)) for name in ('molecules/ethanol.mol2',
        #                                                                      'molecules/cyclohexane.mol2')]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology)

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_1_cyclohexane_1_ethanol_vacuum(self, toolkit_registry, registry_description):
        """Test parametrizing a nonperiodic system of two distinct molecules"""
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_cyclohexane_1_ethanol.pdb'))
        molecules = [create_ethanol(), create_cyclohexane()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        topology.box_vectors = None

        omm_system = forcefield.create_openmm_system(topology)



    @pytest.mark.slow
    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    @pytest.mark.parametrize("box", ['ethanol_water.pdb',
                                     'cyclohexane_water.pdb',
                                     'cyclohexane_ethanol_0.4_0.6.pdb',
                                     'propane_methane_butanol_0.2_0.3_0.5.pdb'])
    def test_parameterize_large_system(self, toolkit_registry, registry_description, box):
        """Test parameterizing a large system of several distinct molecules.
        This test is very slow, so it is only run if the --runslow option is provided to pytest.
        """
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        box_file_path = get_data_file_path(os.path.join('systems', 'packmol_boxes', box))
        pdbfile = app.PDBFile(box_file_path)
        mol_names = ['water', 'cyclohexane', 'ethanol', 'propane', 'methane', 'butanol']
        sdf_files = [get_data_file_path(os.path.join('systems', 'monomers', name+'.sdf')) for name in mol_names]
        molecules = [Molecule.from_file(sdf_file) for sdf_file in sdf_files]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules, )

        omm_system = forcefield.create_openmm_system(topology, toolkit_registry=toolkit_registry)
        # TODO: Add check to ensure system energy is finite


    @pytest.mark.skipif( not(OpenEyeToolkitWrapper.is_available()), reason='Test requires OE toolkit')
    def test_parameterize_ethanol_different_reference_ordering_openeye(self):
        """
        Test parameterizing the same PDB, using reference mol2s that have different atom orderings.
        The results of both should be identical.
        """
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        from simtk.openmm import app
        from simtk.openmm import XmlSerializer

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        # Load the unique molecules with one atom ordering
        molecules1 = [Molecule.from_file(get_data_file_path('molecules/ethanol.sdf'))]
        topology1 = Topology.from_openmm(pdbfile.topology,
                                         unique_molecules=molecules1,
                                         )
        omm_system1 = forcefield.create_openmm_system(topology1,
                                                      toolkit_registry=toolkit_registry)
        # Load the unique molecules with a different atom ordering
        molecules2 = [Molecule.from_file(get_data_file_path('molecules/ethanol_reordered.sdf'))]
        topology2 = Topology.from_openmm(pdbfile.topology,
                                         unique_molecules=molecules2,
                                         )
        omm_system2 = forcefield.create_openmm_system(topology2,
                                                      toolkit_registry=toolkit_registry)

        serialized_1 = XmlSerializer.serialize(omm_system1)
        serialized_2 = XmlSerializer.serialize(omm_system2)

        serialized_1 = round_charge(serialized_1)
        serialized_2 = round_charge(serialized_2)

        assert serialized_1 == serialized_2

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='Test requires RDKit toolkit')
    def test_parameterize_ethanol_different_reference_ordering_rdkit(self):
        """
        Test parameterizing the same PDB, using reference mol2s that have different atom orderings.
        The results of both should be identical.
        """

        from simtk.openmm import app
        from simtk.openmm import XmlSerializer

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper])
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))

        # Load the unique molecules with one atom ordering
        molecules1 = [Molecule.from_file(get_data_file_path('molecules/ethanol.sdf'))]
        topology1 = Topology.from_openmm(pdbfile.topology,
                                         unique_molecules=molecules1,

                                         )
        omm_system1 = forcefield.create_openmm_system(topology1,
                                                      toolkit_registry=toolkit_registry)

        # Load the unique molecules with a different atom ordering
        molecules2 = [Molecule.from_file(get_data_file_path('molecules/ethanol_reordered.sdf'))]
        topology2 = Topology.from_openmm(pdbfile.topology,
                                         unique_molecules=molecules2,
                                         )
        omm_system2 = forcefield.create_openmm_system(topology2,
                                                      toolkit_registry=toolkit_registry)

        serialized_1 = XmlSerializer.serialize(omm_system1)
        serialized_2 = XmlSerializer.serialize(omm_system2)

        serialized_1 = round_charge(serialized_1)
        serialized_2 = round_charge(serialized_2)

        assert serialized_1 == serialized_2


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='Test requires RDKit toolkit')
    def test_parameterize_mol_missing_stereo_rdkit(self):
        """
        Test parameterizing a molecule with undefined stereochemsity using the RDKit/AmberTools backend.
        """

        from openforcefield.topology import Molecule, Topology
        from openforcefield.typing.engines.smirnoff import ForceField
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper])

        molecule = Molecule.from_smiles('CC1CCC(=O)O1', allow_undefined_stereo=True)
        topology = Topology.from_molecules([molecule])

        force_field = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        force_field.create_openmm_system(topology, toolkit_registry=toolkit_registry)


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='Test requires OpenEye toolkit')
    def test_parameterize_mol_missing_stereo_openeye(self):
        """
        Test parameterizing a molecule with undefined stereochemsity using the OpenEye backend.
        """

        from openforcefield.topology import Molecule, Topology
        from openforcefield.typing.engines.smirnoff import ForceField
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])

        molecule = Molecule.from_smiles('CC1CCC(=O)O1', allow_undefined_stereo=True)
        topology = Topology.from_molecules([molecule])

        force_field = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        force_field.create_openmm_system(topology, toolkit_registry=toolkit_registry)

    @pytest.mark.skip(reason="We will not support going directly to ParmEd for now."
                             "We will instead feed OpenMM System objects to ParmEd "
                             "for further processing.")
    def test_parameterize_ethanol_to_parmed(self):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        #toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [ Molecule.from_file(get_data_file_path(name)) for name in ('molecules/ethanol.mol2',) ]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        parmed_system = forcefield.create_parmed_structure(topology, positions=pdbfile.getPositions())




    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_pass_invalid_kwarg_to_create_openmm_system(self, toolkit_registry, registry_description):
        """Test to ensure an exception is raised when an unrecognized kwarg is passed """
        from simtk.openmm import app

        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(file_path)
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        molecules = []
        molecules.append(Molecule.from_smiles('CCO'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        with pytest.raises(ValueError, match=".* not used by any registered force Handler: {'invalid_kwarg'}.*") as e:
            omm_system = forcefield.create_openmm_system(topology, invalid_kwarg='aaa', toolkit_registry=toolkit_registry)


    @pytest.mark.parametrize("inputs", nonbonded_resolution_matrix)
    def test_nonbonded_method_resolution(self,
                                         inputs
                                         ):
        """Test predefined permutations of input options to ensure nonbonded handling is correctly resolved"""
        from simtk.openmm import app

        vdw_method = inputs['vdw_method']
        electrostatics_method = inputs['electrostatics_method']
        has_periodic_box = inputs['has_periodic_box']
        omm_force = inputs['omm_force']
        exception = inputs['exception']
        exception_match= inputs['exception_match']

        molecules = [create_ethanol()]
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        if not(has_periodic_box):
            topology.box_vectors = None

        if exception is None:
            # The method is validated and may raise an exception if it's not supported.
            forcefield.get_parameter_handler('vdW', {}).method = vdw_method
            forcefield.get_parameter_handler('Electrostatics', {}).method = electrostatics_method
            omm_system = forcefield.create_openmm_system(topology)
            nonbond_method_matched = False
            for f_idx in range(omm_system.getNumForces()):
                force = omm_system.getForce(f_idx)
                if isinstance(force, openmm.NonbondedForce):
                    if force.getNonbondedMethod() == omm_force:
                        nonbond_method_matched = True
            assert nonbond_method_matched
        else:
            with pytest.raises(exception, match=exception_match) as excinfo:
                # The method is validated and may raise an exception if it's not supported.
                forcefield.get_parameter_handler('vdW', {}).method = vdw_method
                forcefield.get_parameter_handler('Electrostatics', {}).method = electrostatics_method
                omm_system = forcefield.create_openmm_system(topology)

    def test_parameter_handler_lookup(self):
        """Ensure __getitem__ lookups work"""
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

        handlers_before = sorted(forcefield._parameter_handlers)

        for val in handlers_before:
            looked_up_handler = forcefield[val]
            assert isinstance(looked_up_handler, ParameterHandler)

        handlers_after = sorted(forcefield._parameter_handlers)

        assert handlers_before == handlers_after

    @pytest.mark.parametrize('unregistered_handler', ['LibraryCharges', 'foobar'])
    def test_unregistered_parameter_handler_lookup(self, unregistered_handler):
        """Ensure __getitem__ lookups do not register new handlers"""
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

        assert unregistered_handler not in forcefield._parameter_handlers
        with pytest.raises(KeyError, match=unregistered_handler):
            forcefield[unregistered_handler]
        assert unregistered_handler not in forcefield._parameter_handlers

    def test_lookup_parameter_handler_object(self):
        """Ensure __getitem__ raises NotImplemented when passed a ParameterHandler object"""
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        bonds = forcefield['Bonds']
        with pytest.raises(NotImplementedError):
            forcefield[bonds]
        with pytest.raises(NotImplementedError):
            forcefield[type(bonds)]
        
class TestForceFieldChargeAssignment:

    def generate_monatomic_ions():
        return (('Li+', +1*unit.elementary_charge),
                ('Na+', +1*unit.elementary_charge),
                ('K+', +1*unit.elementary_charge),
                ('Rb+', +1*unit.elementary_charge),
                ('Cs+', +1*unit.elementary_charge),
                ('F-', -1*unit.elementary_charge),
                ('Cl-', -1*unit.elementary_charge),
                ('Br-', -1*unit.elementary_charge),
                ('I-', -1*unit.elementary_charge))

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_charges_from_molecule(self, toolkit_registry, registry_description):
        """Test skipping charge generation and instead getting charges from the original Molecule"""
        # Create an ethanol molecule without using a toolkit
        molecules = [create_ethanol()]

        from simtk.openmm import app, NonbondedForce

        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(file_path)
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        topology = Topology.from_openmm(pdbfile.topology,
                                        unique_molecules=molecules)
        omm_system = forcefield.create_openmm_system(topology,
                                                     charge_from_molecules=molecules,
                                                     toolkit_registry=toolkit_registry)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = ((0, -0.4 * unit.elementary_charge),
                            (1, -0.3 * unit.elementary_charge),
                            (2, -0.2 * unit.elementary_charge),
                            )
        for particle_index, expected_charge in expected_charges:
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

        # In 1_ethanol_reordered.pdb, the first three atoms go O-C-C instead of C-C-O. This part of the test ensures
        # that the charges are correctly mapped according to this PDB in the resulting system.
        pdbfile2 = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol_reordered.pdb'))
        topology2 = Topology.from_openmm(pdbfile2.topology,
                                         unique_molecules=molecules)
        omm_system2 = forcefield.create_openmm_system(topology2,
                                                      charge_from_molecules=molecules,
                                                      toolkit_registry=toolkit_registry)
        nonbondedForce2 = [f for f in omm_system2.getForces() if type(f) == NonbondedForce][0]
        expected_charges2 = ((0, -0.2*unit.elementary_charge),
                             (1, -0.4*unit.elementary_charge),
                             (2, -0.3*unit.elementary_charge),
                            )
        for particle_index, expected_charge in expected_charges2:
            q, sigma, epsilon = nonbondedForce2.getParticleParameters(particle_index)
            assert q == expected_charge


    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_nonintegral_charge_exception(self, toolkit_registry, registry_description):
        """Test skipping charge generation and instead getting charges from the original Molecule"""
        from simtk.openmm import app
        from openforcefield.typing.engines.smirnoff.parameters import NonintegralMoleculeChargeException
        # Create an ethanol molecule without using a toolkit
        ethanol = create_ethanol()
        ethanol.partial_charges[0] = 1. * unit.elementary_charge


        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(file_path)
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        topology = Topology.from_openmm(pdbfile.topology,
                                        unique_molecules=[ethanol])

        # Fail because nonintegral charges aren't allowed
        with pytest.raises(NonintegralMoleculeChargeException,
                           match="Partial charge sum [(]1.40001 e[)] for molecule"):
            omm_system = forcefield.create_openmm_system(topology,
                                                         charge_from_molecules=[ethanol],
                                                         toolkit_registry=toolkit_registry)
        # Pass when the `allow_nonintegral_charges` keyword is included
        omm_system = forcefield.create_openmm_system(topology,
                                                     charge_from_molecules=[ethanol],
                                                     toolkit_registry=toolkit_registry,
                                                     allow_nonintegral_charges=True)


    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_some_charges_from_molecule(self, toolkit_registry, registry_description):
        """
        Test creating an OpenMM system where some charges come from a Molecule, but others come from toolkit
        calculation
        """
        ethanol = create_ethanol()
        cyclohexane = create_cyclohexane()
        molecules = [ethanol, cyclohexane]

        from simtk.openmm import app, NonbondedForce

        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(file_path)
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_cyclohexane_1_ethanol.pdb'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules, )

        omm_system = forcefield.create_openmm_system(topology,
                                                     charge_from_molecules=[ethanol],
                                                     toolkit_registry=toolkit_registry)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = ((18, -0.4 * unit.elementary_charge),
                            (19, -0.3 * unit.elementary_charge),
                            (20, -0.2 * unit.elementary_charge),
                            )
        for particle_index, expected_charge in expected_charges:
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge
        for particle_index in range(topology.n_topology_particles):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q != (0. * unit.elementary_charge)

    def test_library_charges_to_single_water(self):
        """Test assigning charges to one water molecule using library charges"""
        from simtk.openmm import NonbondedForce

        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', 'test_forcefields/tip3p.offxml')
        mol = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers','water.sdf')))
        omm_system = ff.create_openmm_system(mol.to_topology())
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-0.834, 0.417, 0.417] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_parse_library_charges_from_spec_docs(self):
        """Ensure that the examples for librarycharges in the SMIRNOFF spec page are still valid"""
        # TODO: This test is practically useless while the XML strings are hard-coded at the top of this file.
        #       We should implement something like doctests for the XML snippets on the SMIRNOFF spec page.
        ff = ForceField(xml_spec_docs_ala_library_charges_xml)
        ff = ForceField(xml_spec_docs_tip3p_library_charges_xml)

    def test_parse_charge_increment_model_from_spec_docs(self):
        """Ensure that the examples for librarycharges in the SMIRNOFF spec page are still valid"""
        # TODO: This test is practically useless while the XML strings are hard-coded at the top of this file.
        #       We should implement something like doctests for the XML snippets on the SMIRNOFF spec page.
        ff = ForceField(xml_spec_docs_charge_increment_model_xml)

    def test_charge_increment_model_forward_and_reverse_ethanol(self):
        """Test application of ChargeIncrements to the same molecule with different orderings in the topology"""
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge_increment1="-0.05*elementary_charge" charge_increment2="0.05*elementary_charge"/>
            <ChargeIncrement smirks="[C:1][C:2][O:3]" charge_increment1="0.2*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']
        top = Topology.from_molecules([create_ethanol(), create_reversed_ethanol()])
        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0.2, -0.15, -0.05, 0., 0., 0., 0., 0., 0.,
                            0., 0., 0., 0., 0., 0., -0.05, -0.15, 0.2] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge

    def test_charge_increment_model_initialize_with_no_elements(self):
        """Ensure that we can initialize a ForceField object from an OFFXML with a ChargeIncrementModel header, but no
        ChargeIncrement elements"""
        ff = ForceField(xml_charge_increment_model_formal_charges)

    def test_charge_increment_model_net_charge(self):
        """Test application of charge increments on a molecule with a net charge"""
        from simtk import unit
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#6X3:1]-[#8X1-1:2]" charge_increment1="-0.05*elementary_charge" charge_increment2="0.05*elementary_charge"/>
            <ChargeIncrement smirks="[#6X3:1]=[#8X1:2]" charge_increment1="0.2*elementary_charge" charge_increment2="-0.2*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']

        acetate = create_acetate()
        top = acetate.to_topology()
        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0, 0.15, -0.2 , -0.95, 0, 0, 0] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge


    def test_charge_increment_model_deduplicate_symmetric_matches(self):
        """Test that chargeincrementmodelhandler deduplicates symmetric matches"""
        from simtk import unit

        ethanol = create_ethanol()
        top = ethanol.to_topology()

        # Test a charge increment that matches all C-H bonds at once
        # (this should be applied once: C0-H3-H4-H5)
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])([#1:4])" charge_increment1="0.3*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge" charge_increment4="-0.1*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']


        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0.3, 0, 0, -0.1, -0.1, -0.1, 0., 0., 0.] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge

        # Test a charge increment that matches two C-H bonds at a time
        # (this should be applied 3 times: C0-H3-H4, C0-H3-H5, C0-H4-H5)
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])[#6][#8]" charge_increment1="0.1*elementary_charge" charge_increment2="-0.05*elementary_charge" charge_increment3="-0.05*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']


        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0.3, 0, 0, -0.1, -0.1, -0.1, 0., 0., 0.] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge


        # Test a charge increment that matches ONE C-H bond at a time
        # (this should be applied three times: C0-H3, C0-H4, C0-H5)
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#6X4:1]([#1:2])[#6][#8]" charge_increment1="0.1*elementary_charge" charge_increment2="-0.1*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']

        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0.3, 0, 0, -0.1, -0.1, -0.1, 0., 0., 0.] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge

    def test_charge_increment_model_completely_overlapping_matches_override(self):
        """Ensure that DIFFERENT chargeincrements override one another if they apply to the
        same atoms, regardless of order"""
        from simtk import unit
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#1:1]-[#6:2]([#1:3])([#1:4])" charge_increment1="0.123*elementary_charge" charge_increment2="0.369*elementary_charge" charge_increment3="-0.123*elementary_charge" charge_increment4="0.123*elementary_charge"/>
            <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])([#1:4])" charge_increment1="0.3*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge" charge_increment4="-0.1*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']

        ethanol = create_ethanol()
        top = ethanol.to_topology()
        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0.3, 0, 0, -0.1, -0.1, -0.1, 0., 0., 0.] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge


    def test_charge_increment_model_partially_overlapping_matches_both_apply(self):
        """Ensure that DIFFERENT chargeincrements BOTH get applied if they match
        a partially-overlapping set of atoms"""
        from simtk import unit
        test_charge_increment_model_ff = '''
        <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
          <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
          <ChargeIncrementModel version="0.3" number_of_conformers="0" partial_charge_method="formal_charge">
            <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])([#1:4])" charge_increment1="0.3*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge" charge_increment4="-0.1*elementary_charge"/>
            <ChargeIncrement smirks="[#6X4:1][#6X4:2][#8]" charge_increment1="0.05*elementary_charge" charge_increment2="-0.05*elementary_charge"/>
          </ChargeIncrementModel>
        </SMIRNOFF>'''
        file_path = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        ff = ForceField(file_path, test_charge_increment_model_ff)
        del ff._parameter_handlers['ToolkitAM1BCC']

        ethanol = create_ethanol()
        top = ethanol.to_topology()
        sys = ff.create_openmm_system(top)
        nonbonded_force = [force for force in sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        expected_charges = [0.35, -0.05, 0, -0.1, -0.1, -0.1, 0., 0., 0.] * unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert abs(charge - expected_charge) < 1.e-6 * unit.elementary_charge

    @pytest.mark.parametrize("inputs", partial_charge_method_resolution_matrix)
    def test_partial_charge_resolution(self, inputs):
        """Check that the proper partial charge methods are available, and that unavailable partial charge methods
        raise an exception.
        """
        toolkit_wrapper_class = inputs['toolkit']
        if not(toolkit_wrapper_class.is_available()):
            pytest.skip(f"{toolkit_wrapper_class} is not available.")
        toolkit_wrapper = toolkit_wrapper_class()
        partial_charge_method = inputs['partial_charge_method']
        expected_exception = inputs['exception']
        expected_exception_match = inputs['exception_match']
        ethanol = create_ethanol()
        ethanol.generate_conformers()
        if expected_exception is None:
            ethanol.assign_partial_charges(partial_charge_method=partial_charge_method,
                                           toolkit_registry=toolkit_wrapper)
            abs_charge_sum = 0. * unit.elementary_charge

            # Ensure that nonzero charges were assigned
            for pc in ethanol.partial_charges:
                abs_charge_sum += abs(pc)
            assert abs_charge_sum > 0.5 * unit.elementary_charge

        else:
            with pytest.raises(expected_exception, match=expected_exception_match) as excinfo:
                ethanol.assign_partial_charges(partial_charge_method=partial_charge_method,
                                                toolkit_registry=toolkit_wrapper)

    def test_library_charge_hierarchy(self):
        """Test assigning charges to one water molecule using library charges, where two LCs match and the
        assignment is determined by order they are added to the force field"""
        from simtk.openmm import NonbondedForce

        # Test with xml_OH_library_charges_xml loaded last, which should assign dummy partial charges
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                        'test_forcefields/tip3p.offxml',
                        xml_OH_library_charges_xml)
        mol = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers','water.sdf')))
        omm_system = ff.create_openmm_system(mol.to_topology())
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-2., 1., 1.] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

        # Test again, but with tip3p.offxml loaded last (loading the correct partial charges)
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', xml_OH_library_charges_xml, 'test_forcefields/tip3p.offxml', )
        omm_system = ff.create_openmm_system(mol.to_topology())
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-0.834, 0.417, 0.417] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_library_charges_to_two_waters(self):
        """Test assigning charges to two water molecules using library charges"""
        from simtk.openmm import NonbondedForce

        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', 'test_forcefields/tip3p.offxml')
        mol = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers','water.sdf')))
        top = Topology.from_molecules([mol, mol])
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-0.834, 0.417, 0.417, -0.834, 0.417, 0.417] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_library_charges_to_three_ethanols_different_atom_ordering(self):
        """Test assigning charges to three ethanols with different atom orderings"""
        from simtk.openmm import NonbondedForce

        # Define a library charge parameter for ethanol (C1-C2-O3) where C1 has charge -0.2, and its Hs have -0.02,
        # C2 has charge -0.1 and its Hs have -0.01, and O3 has charge 0.3, and its H has charge 0.08

        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', xml_ethanol_library_charges_ff)

        # ethanol.sdf
        #      H5   H8
        #      |    |
        # H6 - C1 - C2 - O3 - H4
        #      |    |
        #      H7   H9
        #
        # ethanol_reordered.sdf (The middle C and O switch indices)
        #      H5   H8
        #      |    |
        # H6 - C1 - C3 - O2 - H4
        #      |    |
        #      H7   H9
        #
        # create_reversed_ethanol()
        #      H5   H2
        #      |    |
        # H4 - C8 - C7 - O6 - H0
        #      |    |
        #      H3   H1


        molecules = [Molecule.from_file(get_data_file_path('molecules/ethanol.sdf')),
                     Molecule.from_file(get_data_file_path('molecules/ethanol_reordered.sdf')),
                     create_reversed_ethanol()]
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-0.2, -0.1, 0.3, 0.08, -0.02, -0.02, -0.02, -0.01, -0.01,
                            -0.2, 0.3, -0.1, 0.08, -0.02, -0.02, -0.02, -0.01, -0.01,
                            0.08, -0.01, -0.01, -0.02, -0.02, -0.02, 0.3, -0.1, -0.2] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge


    @pytest.mark.parametrize('monatomic_ion,formal_charge', generate_monatomic_ions())
    def test_library_charges_monatomic_ions(self, monatomic_ion, formal_charge):
        """Test assigning library charges to each of the monatomic ions in openff-1.1.0.xml"""
        from simtk.openmm import NonbondedForce

        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                        'test_forcefields/ion_charges.offxml')
        mol = Molecule.from_smiles("[{}]".format(monatomic_ion))
        omm_system = ff.create_openmm_system(mol.to_topology())

        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        q, sigma, epsilon = nonbondedForce.getParticleParameters(0)
        assert q == formal_charge

    def test_charge_method_hierarchy(self):
        """Ensure that molecules are parameterized by charge_from_molecules first, then library charges
        if not applicable, then AM1BCC otherwise"""
        from simtk.openmm import NonbondedForce

        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                        xml_CH_zeroes_library_charges_xml,
                        'test_forcefields/tip3p.offxml',
                        xml_charge_increment_model_formal_charges
                        )
        # Cyclohexane will be assigned nonzero charges based on `charge_from_molecules` kwarg
        cyclohexane = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers','cyclohexane.sdf')))
        # Butanol will be assigned zero-valued charges based on `charge_from_molecules` kwarg
        butanol = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers', 'butanol.sdf')))
        # Propane will be assigned charges from AM1BCC
        propane = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers', 'propane.sdf')))
        # Water will be assigned TIP3P librarycharges
        water = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers', 'water.sdf')))
        # Ethanol should be caught by ToolkitAM1BCC
        ethanol = Molecule.from_file(get_data_file_path(os.path.join('systems', 'monomers', 'ethanol.sdf')))
        # iodide should fail to be assigned charges by AM1-BCC, and instead be caught by ChargeIncrementHandler
        # (and assigned formal charge by dummy charge method)
        iodide = Molecule.from_smiles('[I-1]')

        # Assign dummy partial charges to cyclohexane, which we expect to find in the final system since it
        # is included in the charge_from_molecules kwarg to create_openmm_system
        cyclohexane.partial_charges = np.array([-0.2, -0.2, -0.2, -0.2, -0.2, -0.2,
                                             0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                             0.1, 0.1, 0.1, 0.1, 0.1, 0.1]) * unit.elementary_charge

        # There were previously known issues when parameterizing molecules with all zero charges,
        # so test this explicitly with butanol. Since butanol will be in the charge_from_molecules kwarg,
        # we expect to find these charges in the final system.
        butanol.partial_charges = np.array([0.0] * 15) * unit.elementary_charge

        # Add dummy partial charges to propane, which should be IGNORED since it
        # isn't in the charge_from_molecules kwarg
        propane.partial_charges = np.array([99.] * 11) * unit.elementary_charge

        # Add dummy partial charges to water, which should be IGNORED since it
        # isn't in the charge_from_molecules kwarg
        water.partial_charges = np.array([99.] * 3) * unit.elementary_charge

        #            molecule       correct charge method
        molecules = [cyclohexane, # charge_from_molecules kwarg
                     butanol,     # charge_from_molecules kwarg
                     propane,     # library charges
                     water,       # library charges
                     ethanol,     # AM1-BCC
                     iodide]      # charge increment model (formal charge)
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top, charge_from_molecules=[cyclohexane, butanol])
        existing = [f for f in omm_system.getForces() if type(f) == NonbondedForce]

        # Ensure that the handlers do not make multiple NonbondedForce objects
        assert len(existing) == 1
        nonbondedForce = existing[0]
        expected_charges = [# cyclohexane (18 atoms) should have the following values from charge_from_mols
                            -0.2, -0.2, -0.2, -0.2, -0.2, -0.2,
                             0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                             0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                            # butanol (15 atoms) should have the following values from charge_from_mols
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0,
                            # propane (11 atoms) should have the following values from xml_CH_zeroes_library_charges_xml
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0,
                            # water (3 atoms) should have the following charges from tip3p.offxml
                            -0.834, 0.417, 0.417] * unit.elementary_charge

        # Ensure that the first four molecules have exactly the charges we intended
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

        # Ensure the last molecule (ethanol) had _some_ nonzero charge assigned by an AM1BCC implementation
        for particle_index in range(len(expected_charges), top.n_topology_atoms-1):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q != 0 * unit.elementary_charge

        # Ensure that iodine has a charge of -1, specified by charge increment model charge_method="formal charge"
        q, sigma, epsilon = nonbondedForce.getParticleParameters(top.n_topology_atoms-1)
        assert q == -1. * unit.elementary_charge


    def test_assign_charges_to_molecule_in_parts_using_multiple_library_charges(self):
        """Test assigning charges to parts of a molecule using two library charge lines. Note that these LibraryCharge
        SMIRKS have partial overlap, so this also tests that the hierarchy is correctly obeyed."""
        from simtk.openmm import NonbondedForce
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', xml_ethanol_library_charges_in_parts_ff)


        molecules = [Molecule.from_file(get_data_file_path('molecules/ethanol.sdf')),
                     Molecule.from_file(get_data_file_path('molecules/ethanol_reordered.sdf'))]
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-0.2, -0.1, 0.3, 0.08, -0.02, -0.02, -0.02, -0.01, -0.01, -0.2,
                            0.3, -0.1, 0.08, -0.02, -0.02, -0.02, -0.01, -0.01] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_assign_charges_using_library_charges_by_single_atoms(self):
        """Test assigning charges to parts of a molecule using per-atom library charges. Note that these LibraryCharge
        SMIRKS will match multiple atoms, so this is also a test of correct usage of the parameter hierarchy.."""
        from simtk.openmm import NonbondedForce
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', xml_ethanol_library_charges_by_atom_ff)

        molecules = [Molecule.from_file(get_data_file_path('molecules/ethanol.sdf')),
                     Molecule.from_file(get_data_file_path('molecules/ethanol_reordered.sdf'))]
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        expected_charges = [-0.2, -0.1, 0.3, 0.08, -0.02, -0.02, -0.02, -0.01, -0.01, -0.2,
                            0.3, -0.1, 0.08, -0.02, -0.02, -0.02, -0.01, -0.01] * unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_library_charges_dont_parameterize_molecule_because_of_incomplete_coverage(self):
        """Fail to assign charges to a molecule because not all atoms can be assigned"""
        from simtk.openmm import NonbondedForce
        from openforcefield.typing.engines.smirnoff.parameters import UnassignedMoleculeChargeException

        molecules = [Molecule.from_file(get_data_file_path('molecules/toluene.sdf'))]
        top = Topology.from_molecules(molecules)

        # The library charges in the FF should not be able to fully cover toluene
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', xml_ethanol_library_charges_by_atom_ff)
        # Delete the ToolkitAM1BCCHandler so the molecule won't get charges from anywhere
        del ff._parameter_handlers['ToolkitAM1BCC']
        with pytest.raises(UnassignedMoleculeChargeException,
                           match="did not have charges assigned by any ParameterHandler") as excinfo:
            omm_system = ff.create_openmm_system(top)

        # If we do NOT delete the ToolkiAM1BCCHandler, then toluene should be assigned some nonzero partial charges.
        # The exact value will vary by toolkit, so we don't test that here.
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml', xml_ethanol_library_charges_by_atom_ff)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [f for f in omm_system.getForces() if type(f) == NonbondedForce][0]
        for particle_index in range(top.n_topology_atoms):
            q, sigma, epsilon = nonbondedForce.getParticleParameters(particle_index)
            assert q != 0 * unit.elementary_charge




#======================================================================
# TEST CONSTRAINTS
#======================================================================

class TestForceFieldConstraints:
    """Tests that constraints are correctly applied and behave correctly."""

    @classmethod
    def check_molecule_constraints(cls, molecule, system, bond_elements, bond_length):
        """Check that the bonds in the molecule is correctly constrained."""
        for constraint_idx in range(system.getNumConstraints()):
            atom1_idx, atom2_idx, distance = system.getConstraintParameters(constraint_idx)
            atom_elements = {molecule.atoms[atom1_idx].element.symbol,
                             molecule.atoms[atom2_idx].element.symbol}
            assert atom_elements == bond_elements
            assert np.isclose(distance/unit.angstrom, bond_length/unit.angstrom)

    def test_constraints_hbonds(self):
        """Test that hydrogen bonds constraints are applied correctly to a ethane molecule."""
        # Parametrize an ethane molecule.
        ethane = Molecule.from_smiles('CC')
        topology = Topology.from_molecules([ethane])
        ff = ForceField(XML_FF_GENERICS, 'test_forcefields/old/hbonds.offxml')
        system = ff.create_openmm_system(topology)

        # Check that all C-H bonds have been constrained to the FF bond length.
        self.check_molecule_constraints(ethane, system,
                                        bond_elements={'C', 'H'},
                                        bond_length= 1.09 * unit.angstrom)


#======================================================================
# TEST PARAMETER ASSIGNMENT
#======================================================================

def generate_alkethoh_parameters_assignment_cases():
    """Create dynamically all test cases that should be ran for the AlkEthOH set."""
    # These AlkEthOH molecules are always run by test_alkethoh_parameters_assignment.
    fast_test_cases = [
        'r0',
        'r12',
        'r118',
        'c38',
        'c100',
        'c1161',
        'c1266'
    ]

    def extract_id(file_path):
        """Extract the AlkEthOH molecule ID from the file path."""
        # An example of file path is AlkEthOH_tripos/AlkEthOH_chain_filt1/AlkEthOH_c555.crd
        return os.path.splitext(os.path.basename(file_path))[0][9:]

    # Get all the molecules ids from the tarfiles. The tarball is extracted
    # in conftest.py if slow tests are activated.
    import tarfile
    alkethoh_tar_file_path = get_data_file_path(os.path.join('molecules', 'AlkEthOH_tripos.tar.gz'))
    with tarfile.open(alkethoh_tar_file_path, 'r:gz') as tar:
        # Collect all the files discarding the duplicates in the test_filt1 folder.
        slow_test_cases = {extract_id(m.name) for m in tar.getmembers()
                           if 'crd' in m.name and 'test_filt1' not in m.name}

    # Remove fast test cases from slow ones to avoid duplicate tests.
    # Remove also water (c1302), which was reparameterized in AlkEthOH
    # to be TIP3P (not covered by Frosst_AlkEthOH_parmAtFrosst.
    for fast_test_case in fast_test_cases + ['c1302']:
        slow_test_cases.remove(fast_test_case)

    # Mark all slow cases as slow.
    slow_test_cases = [pytest.param(case, marks=pytest.mark.slow)
                       for case in sorted(slow_test_cases)]

    # Isolate the AlkEthOH ID.
    return fast_test_cases + slow_test_cases


def generate_freesolv_parameters_assignment_cases():
    """Create dynamically all test cases that should be ran for the FreeSolv set."""
    import tarfile

    # For these tests, UndefinedStereochemistryError is ignored.
    # The chirality was manually checked (see issue #175).
    ignore_undefined_stereo = {
        '2501588',
        '3266352',
        '7829570'
    }

    # These molecules are always tested by test_freesolv_parameters_assignment().
    # Each test case is (freesolv_id, force_field_version, allow_undefined_stereo).
    fast_test_cases = [
        ('1019269', '0_0_4_fixed', False),
        ('63712', '0_0_2', False),  # The XML was regenerated after fixing the issue described in #179.
        ('1723043', '0_0_2', False),
        ('2501588', '0_0_2', True),  # Test impropers and undefined stereochemistry.
        ('3323117', '0_0_2', False),  # The XML was regenerated after fixing the issue described in #179.
        ('1107178', '0_0_2', False),  # Molecule with iodine
        ('1036761', '0_0_2', False),  # Molecule with primary amine
    ]

    def extract_id(file_path):
        """Extract the FreeSolv ID and force field version from the file subpath."""
        # An example of file path is FreeSolv/xml_0_0_4_fixed/mobley_7913234_vacuum.xml
        freesolv_id = os.path.basename(file_path).split('_')[1]
        force_field_version = os.path.basename(os.path.dirname(file_path))[4:]
        allow_undefined_stereo = freesolv_id in ignore_undefined_stereo
        return (freesolv_id, force_field_version, allow_undefined_stereo)

    # Get all the tarball XML files available. The tarball is extracted
    # in conftest.py if slow tests are activated.
    freesolv_tar_file_path = get_data_file_path(os.path.join('molecules', 'FreeSolv.tar.gz'))
    with tarfile.open(freesolv_tar_file_path, 'r:gz') as tar:
        slow_test_cases = {extract_id(m.name) for m in tar.getmembers() if '.xml' in m.name}

    # Remove fast test cases from slow ones to avoid duplicate tests.
    for fast_test_case in fast_test_cases:
        slow_test_cases.remove(fast_test_case)

    # Mark all slow cases as slow.
    slow_test_cases = [pytest.param(*case, marks=pytest.mark.slow)
                       for case in sorted(slow_test_cases)]

    return fast_test_cases + slow_test_cases


class TestForceFieldParameterAssignment:
    """Regression tests checking that parameters are assigned correctly."""

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    @pytest.mark.parametrize('alkethoh_id', generate_alkethoh_parameters_assignment_cases())
    def test_alkethoh_parameters_assignment(self, alkethoh_id):
        """Test that ForceField assign parameters correctly in the AlkEthOH set.
        The test compares the System parameters of a AlkEthOH molecule
        parameterized with AMBER and Frosst_AlkEthOH_parmAtFrosst.offxml.

        The AMBER files were prepared following the pipeline described here:
            https://github.com/openforcefield/open-forcefield-data/tree/master/Model-Systems/AlkEthOH_distrib/
        They were built for the SMIRNOFF parametrization to yield exact same
        parameters.

        The AlkEthOH set, however, does not have impropers, which should be
        tested separately. Currently, test_freesolv_parameters_assignment
        does the job.

        """
        from openforcefield.tests.utils import get_alkethoh_file_path, compare_amber_smirnoff

        # Obtain the path to the input files.
        alkethoh_name = 'AlkEthOH_' + alkethoh_id
        mol2_filepath, top_filepath, crd_filepath = get_alkethoh_file_path(alkethoh_name, get_amber=True)

        # Load molecule.
        molecule = Molecule.from_file(mol2_filepath)

        # Load forcefield
        forcefield = ForceField('test_forcefields/Frosst_AlkEthOH_parmAtFrosst.offxml')

        # Compare parameters. Skip the energy checks as the parameter check should be
        # sufficient. We test both energies and parameters in the slow test.
        # We ignore the charges for now as they are not included in the force field.
        # TODO: Reactivate the charge check when we'll be able to load charges from files.
        compare_amber_smirnoff(top_filepath, crd_filepath, forcefield, molecule,
                               check_energies=False, ignore_charges=True)


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    def test_multi_alkethoh_parameters_assignment(self):
        """Test that systems with multiple reference molecules are parametrized correctly.

        The test relies on the fact that we have already verified we can
        parametrize correctly single AlkEthOH molecules in
        test_alkethoh_parameters_assignment(). We use ParmEd to merge
        the AMBER files to be used as reference parameters.

        """
        import parmed
        from openforcefield.tests.utils import (get_alkethoh_file_path,
                                                compare_system_parameters,
                                                compare_system_energies)

        # The AlkEthOH molecule ids to mix in the systems.
        alketoh_ids = ['r0', 'c38', 'c1161']

        # Load molecules and structures.
        molecules = []
        structures = []
        for alkethoh_id in alketoh_ids:
            mol2_filepath, top_filepath, crd_filepath = get_alkethoh_file_path(
                'AlkEthOH_'+alkethoh_id, get_amber=True)
            molecules.append(Molecule.from_file(mol2_filepath))
            amber_parm = parmed.load_file(top_filepath, crd_filepath)
            # Convert this into a real structure as mixing AmberParm objects is bugged (see ParmEd#1045).
            structures.append(amber_parm.copy(parmed.Structure))

        # Merge the structures into a single system with two copies of the last molecule.
        structure_mixture = structures[0] + structures[1] + structures[2] + structures[-1]
        amber_system = structure_mixture.createSystem(nonbondedMethod=openmm.app.NoCutoff)

        # Create the OpenFF System through ForceField.
        topology = Topology.from_openmm(structure_mixture.topology, unique_molecules=molecules)
        topology.box_vectors = None
        ff = ForceField('test_forcefields/Frosst_AlkEthOH_parmAtFrosst.offxml')
        off_system = ff.create_openmm_system(topology)

        # Translate the molecules a little to avoid overlapping atoms.
        positions = copy.deepcopy(structure_mixture.positions)
        translate_vectors = [
            np.array([1.0, 0.0, 0.0])*unit.nanometer,
            np.array([0.0, 1.0, 0.0])*unit.nanometer,
            np.array([0.0, 0.0, 1.0])*unit.nanometer,
            # Leave the fourth molecule where it is.
        ]
        current_atom_idx = 0
        for mol_idx, (translate_vector, mol) in enumerate(zip(translate_vectors, molecules)):
            n_mol_atoms = len(mol.atoms)
            positions[current_atom_idx:current_atom_idx+n_mol_atoms] += translate_vector
            current_atom_idx += n_mol_atoms

        # Compare parameters and systems.
        # TODO: Reactivate charges comparison when we'll be able to read them from the file.
        compare_system_parameters(amber_system, off_system,
                                  systems_labels=('AMBER', 'SMIRNOFF'),
                                  ignore_charges=True)
        compare_system_energies(amber_system, off_system, positions,
                                ignore_charges=True)


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    @pytest.mark.parametrize(('freesolv_id', 'forcefield_version', 'allow_undefined_stereo'),
                             generate_freesolv_parameters_assignment_cases())
    def test_freesolv_parameters_assignment(self, freesolv_id, forcefield_version, allow_undefined_stereo):
        """Regression test on parameters assignment based on the FreeSolv set used in the 0.1 paper.

        This, contrarily to the similar AlkEthOH test, checks also constraints
        and improper torsions.

        """
        from openforcefield.tests.utils import get_freesolv_file_path, compare_system_parameters
        mol2_file_path, xml_file_path = get_freesolv_file_path(freesolv_id, forcefield_version)

        # Load molecules.
        molecule = Molecule.from_file(mol2_file_path, allow_undefined_stereo=allow_undefined_stereo)

        # Create OpenFF System with the current toolkit.
        forcefield_file_path = 'test_forcefields/old/test_ff_' + forcefield_version + '_spec_0_2.offxml'
        ff = ForceField(forcefield_file_path, 'test_forcefields/old/hbonds.offxml')
        ff_system = ff.create_openmm_system(molecule.to_topology())

        # Load OpenMM System created with the 0.1 version of the toolkit.
        from simtk import openmm
        with open(xml_file_path, 'r') as f:
            xml_system = openmm.XmlSerializer.deserialize(f.read())

        # Compare parameters. We ignore the improper folds as in 0.0.3 we
        # used a six-fold implementation while we now use a three-fold way.
        # TODO: Reactivate charge comparison once we'll be able to read them from file.
        compare_system_parameters(ff_system, xml_system,
                                  systems_labels=('current OpenFF', 'SMIRNOFF 0.0.4'),
                                  ignore_charges=True, ignore_improper_folds=True)



    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    @pytest.mark.parametrize(('is_periodic'), (False, True))
    @pytest.mark.parametrize(('gbsa_model'), ['HCT', 'OBC1', 'OBC2'])
    @pytest.mark.parametrize(('freesolv_id', 'forcefield_version', 'allow_undefined_stereo'),
                             generate_freesolv_parameters_assignment_cases())
    def test_freesolv_gbsa_energies(self, gbsa_model, is_periodic, freesolv_id, forcefield_version, allow_undefined_stereo):
        """
        Regression test on HCT, OBC1, and OBC2 GBSA models. This test ensures that the
        SMIRNOFF-based GBSA models match the equivalent OpenMM implementations.
        """

        from openforcefield.tests.utils import (get_freesolv_file_path,
                                                compare_system_energies, create_system_from_amber,
                                                get_context_potential_energy
                                                )
        import parmed as pmd
        from simtk import openmm
        from simtk.openmm import Platform

        mol2_file_path, _ = get_freesolv_file_path(freesolv_id, forcefield_version)

        # Load molecules.
        molecule = Molecule.from_file(mol2_file_path, allow_undefined_stereo=allow_undefined_stereo)

        # Give each atom a unique name, otherwise OpenMM will complain
        for idx, atom in enumerate(molecule.atoms):
            atom.name = f'{atom.element.symbol}{idx}'
        positions = molecule.conformers[0]

        off_gbsas  = {'HCT': 'test_forcefields/GBSA_HCT-1.0.offxml',
                      'OBC1': 'test_forcefields/GBSA_OBC1-1.0.offxml',
                      'OBC2': 'test_forcefields/GBSA_OBC2-1.0.offxml'
                      }
        # Create OpenFF System with the current toolkit.
        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                        off_gbsas[gbsa_model])
        off_top = molecule.to_topology()
        if is_periodic:
            off_top.box_vectors = ((30., 0, 0), (0, 30., 0), (0, 0, 30.)) * unit.angstrom

        else:
            off_top.box_vectors = None

        off_omm_system = ff.create_openmm_system(off_top, charge_from_molecules=[molecule])

        off_nonbonded_force = [force for force in off_omm_system.getForces() if
                               isinstance(force, openmm.NonbondedForce)][0]


        omm_top = off_top.to_openmm()
        pmd_struct = pmd.openmm.load_topology(omm_top, off_omm_system, positions)
        prmtop_file = NamedTemporaryFile(suffix='.prmtop')
        inpcrd_file = NamedTemporaryFile(suffix='.inpcrd')
        pmd_struct.save(prmtop_file.name, overwrite=True)
        pmd_struct.save(inpcrd_file.name, overwrite=True)

        openmm_gbsas = {'HCT': openmm.app.HCT,
                        'OBC1': openmm.app.OBC1,
                        'OBC2': openmm.app.OBC2,
                        }

        # The functional form of the nonbonded force will change depending on whether the cutoff
        # is None during initialization. Therefore, we need to figure that out here.

        # WARNING: The NonbondedMethod enums at openmm.app.forcefield and openmm.CustomGBForce
        # aren't necessarily the same, and could be misinterpreted if the wrong one is used. For
        # create_system_from_amber, we must provide the app.forcefield version.

        if is_periodic:
            amber_nb_method = openmm.app.forcefield.CutoffPeriodic
            amber_cutoff = off_nonbonded_force.getCutoffDistance()
        else:
            amber_nb_method = openmm.app.forcefield.NoCutoff
            amber_cutoff = None

        (amber_omm_system,
         amber_omm_topology,
         amber_positions) = create_system_from_amber(prmtop_file.name,
                                                     inpcrd_file.name,
                                                     implicitSolvent=openmm_gbsas[gbsa_model],
                                                     nonbondedMethod=amber_nb_method,
                                                     nonbondedCutoff=amber_cutoff,
                                                     gbsaModel='ACE',
                                                     implicitSolventKappa=0.,
                                                     )

        # Retrieve the GBSAForce from both the AMBER and OpenForceField systems
        off_gbsa_forces = [force for force in off_omm_system.getForces() if
                              (isinstance(force, openmm.GBSAOBCForce) or
                               isinstance(force, openmm.openmm.CustomGBForce))]
        assert len(off_gbsa_forces) == 1
        off_gbsa_force = off_gbsa_forces[0]
        amber_gbsa_forces = [force for force in amber_omm_system.getForces() if
                                 (isinstance(force, openmm.GBSAOBCForce) or
                                  isinstance(force, openmm.openmm.CustomGBForce))]
        assert len(amber_gbsa_forces) == 1
        amber_gbsa_force = amber_gbsa_forces[0]

        # We get radius and screen values from each model's getStandardParameters method
        if gbsa_model == 'HCT':
            gb_params = openmm.app.internal.customgbforces.GBSAHCTForce.getStandardParameters(omm_top)
        elif gbsa_model == 'OBC1':
            gb_params = openmm.app.internal.customgbforces.GBSAOBC1Force.getStandardParameters(omm_top)
        elif gbsa_model == 'OBC2':
            gb_params = openmm.app.internal.customgbforces.GBSAOBC2Force.getStandardParameters(omm_top)

        # Use GB params from OpenMM GBSA classes to populate parameters
        for idx, (radius, screen) in enumerate(gb_params):
            # Keep the charge, but throw out the old radius and screen values
            q, old_radius, old_screen = amber_gbsa_force.getParticleParameters(idx)
            if isinstance(amber_gbsa_force, openmm.GBSAOBCForce):
                # Note that in GBSAOBCForce, the per-particle parameters are separate
                # arguments, while in CustomGBForce they're a single iterable
                amber_gbsa_force.setParticleParameters(idx, q, radius, screen)

            elif isinstance(amber_gbsa_force, openmm.CustomGBForce):
                # !!! WARNING: CustomAmberGBForceBase expects different per-particle parameters
                # depending on whether you use addParticle or setParticleParameters. In
                # setParticleParameters, we have to apply the offset and scale BEFORE setting
                # parameters, whereas in addParticle, it is applied afterwards, and the particle
                # parameters are not set until an auxillary finalize() method is called. !!!
                amber_gbsa_force.setParticleParameters(idx, (q, radius - 0.009, screen * (radius - 0.009)))


        # Put the GBSA force into a separate group so we can specifically compare GBSA energies
        amber_gbsa_force.setForceGroup(1)
        off_gbsa_force.setForceGroup(1)

        # Some manual overrides to get the OFF system's NonbondedForce matched up with the AMBER system
        if is_periodic:
            off_nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            off_nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

        off_nonbonded_force.setReactionFieldDielectric(1.0)

        # Not sure if zeroing the switching width is essential -- This might only make a difference
        # in the energy if we tested on a molecule larger than the 9A cutoff
        #off_nonbonded_force.setSwitchingDistance(0)

        # Create Contexts
        integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        platform = Platform.getPlatformByName('Reference')
        amber_context = openmm.Context(amber_omm_system, integrator, platform)
        off_context = openmm.Context(off_omm_system, copy.deepcopy(integrator), platform)

        # Get context energies
        amber_energy = get_context_potential_energy(amber_context, positions)
        off_energy = get_context_potential_energy(off_context, positions)

        # Very handy for debugging
        # print(openmm.XmlSerializer.serialize(off_gbsa_force))
        # print(openmm.XmlSerializer.serialize(amber_gbsa_force))


        # Ensure that the GBSA energies (which we put into ForceGroup 1) are identical
        # For Platform=OpenCL, we do get "=="-level identical numbers, but for "Reference", we don't.
        #assert amber_energy[1] == off_energy[1]
        assert abs(amber_energy[1] - off_energy[1]) < 1e-5 * unit.kilojoule/unit.mole

        # Ensure that all system energies are the same
        compare_system_energies(off_omm_system, amber_omm_system, positions, by_force_type=False)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    @pytest.mark.parametrize('zero_charges', [True, False])
    @pytest.mark.parametrize(('gbsa_model'), ['HCT', 'OBC1', 'OBC2'])
    def test_molecule_energy_gb_no_sa(self, zero_charges, gbsa_model):
        """Test creating a GBSA system without a surface energy term, and validate its energy
        against the same system made using OpenMM's AMBER GBSA functionality"""
        from openforcefield.tests.utils import (compare_system_energies, create_system_from_amber,
                                                get_context_potential_energy
                                                )
        import parmed as pmd
        from simtk import openmm
        from simtk.openmm import Platform
        import numpy as np

        # Load an arbitrary molecule from the freesolv set
        molecule = Molecule.from_file(get_data_file_path('molecules/FreeSolv/mol2files_sybyl/mobley_1036761.mol2'))

        molecule.name = 'mobley_1036761' # Name the molecule, otherwise OpenMM will complain
        if zero_charges:
            molecule.partial_charges = np.zeros(molecule.n_atoms) * unit.elementary_charge

        # Give each atom a unique name, otherwise OpenMM will complain
        for idx, atom in enumerate(molecule.atoms):
            atom.name = f'{atom.element.symbol}{idx}'

        positions = np.concatenate((molecule.conformers[0], molecule.conformers[0] + (10 * unit.angstrom)))
        # Create OpenFF System with the current toolkit.

        off_gbsas  = {'HCT': 'test_forcefields/GBSA_HCT-1.0.offxml',
                      'OBC1': 'test_forcefields/GBSA_OBC1-1.0.offxml',
                      'OBC2': 'test_forcefields/GBSA_OBC2-1.0.offxml'
                      }

        ff = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                        off_gbsas[gbsa_model])
        ff.get_parameter_handler('GBSA').sa_model = None
        off_top = Topology.from_molecules([molecule, molecule])
        off_omm_system = ff.create_openmm_system(off_top, charge_from_molecules=[molecule])

        omm_top = off_top.to_openmm()
        pmd_struct = pmd.openmm.load_topology(omm_top, off_omm_system, positions)
        prmtop_file = NamedTemporaryFile(suffix='.prmtop')
        inpcrd_file = NamedTemporaryFile(suffix='.inpcrd')
        pmd_struct.save(prmtop_file.name, overwrite=True)
        pmd_struct.save(inpcrd_file.name, overwrite=True)

        openmm_gbsas = {'HCT': openmm.app.HCT,
                        'OBC1': openmm.app.OBC1,
                        'OBC2': openmm.app.OBC2,
                        }

        (amber_omm_system,
         amber_omm_topology,
         amber_positions) = create_system_from_amber(prmtop_file.name,
                                                     inpcrd_file.name,
                                                     implicitSolvent=openmm_gbsas[gbsa_model],
                                                     nonbondedMethod = openmm.app.forcefield.NoCutoff,
                                                     nonbondedCutoff=None,
                                                     gbsaModel=None,
                                                     implicitSolventKappa=0.,
                                                     )

        # Retrieve the GBSAForce from both the AMBER and OpenForceField systems
        off_gbsa_forces = [force for force in off_omm_system.getForces() if
                              (isinstance(force, openmm.GBSAOBCForce) or
                               isinstance(force, openmm.openmm.CustomGBForce))]
        assert len(off_gbsa_forces) == 1
        off_gbsa_force = off_gbsa_forces[0]
        amber_gbsa_forces = [force for force in amber_omm_system.getForces() if
                                 (isinstance(force, openmm.GBSAOBCForce) or
                                  isinstance(force, openmm.openmm.CustomGBForce))]
        assert len(amber_gbsa_forces) == 1
        amber_gbsa_force = amber_gbsa_forces[0]

        # We get radius and screen values from each model's getStandardParameters method
        if gbsa_model == 'HCT':
            gb_params = openmm.app.internal.customgbforces.GBSAHCTForce.getStandardParameters(omm_top)
        elif gbsa_model == 'OBC1':
            gb_params = openmm.app.internal.customgbforces.GBSAOBC1Force.getStandardParameters(omm_top)
        elif gbsa_model == 'OBC2':
            gb_params = openmm.app.internal.customgbforces.GBSAOBC2Force.getStandardParameters(omm_top)
            # This is only necessary until https://github.com/openmm/openmm/pull/2362 is bundled into a conda release
            amber_gbsa_force.setSurfaceAreaEnergy(0)

        # Use GB params from OpenMM GBSA classes to populate parameters
        for idx, (radius, screen) in enumerate(gb_params):
            # Keep the charge, but throw out the old radius and screen values
            q, old_radius, old_screen = amber_gbsa_force.getParticleParameters(idx)

            if isinstance(amber_gbsa_force, openmm.GBSAOBCForce):
                # Note that in GBSAOBCForce, the per-particle parameters are separate
                # arguments, while in CustomGBForce they're a single iterable
                amber_gbsa_force.setParticleParameters(idx, q, radius, screen)

            elif isinstance(amber_gbsa_force, openmm.CustomGBForce):
                # !!! WARNING: CustomAmberGBForceBase expects different per-particle parameters
                # depending on whether you use addParticle or setParticleParameters. In
                # setParticleParameters, we have to apply the offset and scale BEFORE setting
                # parameters, whereas in addParticle, it is applied afterwards, and the particle
                # parameters are not set until an auxillary finalize() method is called. !!!
                amber_gbsa_force.setParticleParameters(idx, (q, radius - 0.009, screen * (radius - 0.009)))

        # Put the GBSA force into a separate group so we can specifically compare GBSA energies
        amber_gbsa_force.setForceGroup(1)
        off_gbsa_force.setForceGroup(1)


        # Create Contexts
        integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        platform = Platform.getPlatformByName('Reference')
        amber_context = openmm.Context(amber_omm_system, integrator, platform)
        off_context = openmm.Context(off_omm_system, copy.deepcopy(integrator), platform)

        # Get context energies
        amber_energy = get_context_potential_energy(amber_context, positions)
        off_energy = get_context_potential_energy(off_context, positions)

        # Very handy for debugging
        # print(openmm.XmlSerializer.serialize(off_gbsa_force))
        # print(openmm.XmlSerializer.serialize(amber_gbsa_force))

        # Ensure that the GBSA energies (which we put into ForceGroup 1) are identical
        # For Platform=OpenCL, we do get "=="-level identical numbers, but for "Reference", we don't.
        #assert amber_energy[1] == off_energy[1]
        assert abs(amber_energy[1] - off_energy[1]) < 1e-5 * unit.kilojoule/unit.mole

        # If charges are zero, the GB energy component should be 0, so the total GBSA energy should be 0
        if zero_charges:
            assert amber_energy[1] == 0. * unit.kilojoule / unit.mole
        else:
            assert amber_energy[1] != 0. * unit.kilojoule / unit.mole

        # Ensure that all system energies are the same
        compare_system_energies(off_omm_system, amber_omm_system, positions, by_force_type=False)

    @pytest.mark.slow
    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_protein(self, toolkit_registry, registry_description):
        """Test that ForceField assigns parameters correctly for a protein
        """

        mol_path = get_data_file_path('proteins/T4-protein.mol2')
        molecule = Molecule.from_file(mol_path, allow_undefined_stereo=False)
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        topology = Topology.from_molecules(molecule)


        labels = forcefield.label_molecules(topology)[0]
        assert len(labels["Bonds"]) ==            2654
        assert len(labels["Angles"]) ==           4789
        assert len(labels["ProperTorsions"]) ==   6973
        assert len(labels["ImproperTorsions"]) == 528

        fn = forcefield.create_openmm_system
        omm_system = fn(topology,
                        charge_from_molecules=[molecule],
                        toolkit_registry=toolkit_registry,
                        allow_nonintegral_charges=False)

    @pytest.mark.parametrize(
            ('get_molecule', 'k_interpolated', 'central_atoms'),
            [(create_ethanol, 4.953856, (1, 2)), (create_reversed_ethanol, 4.953856, (7, 6))])
    def test_fractional_bondorder(self, get_molecule, k_interpolated, central_atoms):
        """Test the fractional bond orders are used to interpolate k values as we expect"""

        mol = get_molecule()

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                                xml_ff_torsion_bo)
        topology = Topology.from_molecules(mol)

        omm_system = forcefield.create_openmm_system(
                topology,
                charge_from_molecules=[mol],
                partial_bond_orders_from_molecules=[mol])

        off_torsion_force = [force for force in omm_system.getForces() if
                               isinstance(force, openmm.PeriodicTorsionForce)][0]

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if (((atom2 == atom2_mol) and (atom3 == atom3_mol)) or
                    ((atom2 == atom3_mol) and (atom3 == atom2_mol))):
                k = params[-1]
                assert_almost_equal(k/k.unit, k_interpolated)

    def test_fractional_bondorder_multiple_same_mol(self):
        """Check that an error is thrown when essentially the same molecule is entered more than once
        for partial_bond_orders_from_molecules"""

        mol = create_ethanol()
        mol2 = create_reversed_ethanol()

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                                xml_ff_torsion_bo)
        topology = Topology.from_molecules([mol, mol2])

        with pytest.raises(ValueError):
            omm_system = forcefield.create_openmm_system(
                    topology,
                    charge_from_molecules=[mol],
                    partial_bond_orders_from_molecules=[mol, mol2])

    @pytest.mark.parametrize(
            ('get_molecule', 'central_atoms'),
            [(create_ethanol, (1, 2)), (create_reversed_ethanol, (7, 6))])
    def test_fractional_bondorder_superseded_by_standard_torsion(self, get_molecule, central_atoms):
        """Test that matching rules are still respected with fractional bondorder parameters.

        We check here that a standard torsion parameter can still get priority on a match
        if it is further down the list.
        """

        mol = get_molecule()

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                                xml_ff_torsion_bo_standard_supersede)
        topology = Topology.from_molecules([mol])

        omm_system = forcefield.create_openmm_system(
                topology,
                charge_from_molecules=[mol],
                partial_bond_orders_from_molecules=[mol])

        off_torsion_force = [force for force in omm_system.getForces() if
                               isinstance(force, openmm.PeriodicTorsionForce)][0]

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if (((atom2 == atom2_mol) and (atom3 == atom3_mol)) or
                    ((atom2 == atom3_mol) and (atom3 == atom2_mol))):
                k = params[-1]
                assert_almost_equal(k/k.unit, 5.0208)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='Test requires RDKit toolkit')
    @pytest.mark.parametrize(
            ('get_molecule', 'k_interpolated', 'central_atoms'),
            [(create_ethanol, 4.953856, (1, 2)), (create_reversed_ethanol, 4.953856, (7, 6))])
    def test_fractional_bondorder_calculated_rdkit(self, get_molecule, k_interpolated, central_atoms):
        """Test that torsion barrier heights interpolated from fractional bond
        orders calculated with the Amber toolkit are assigned within our
        expectations.

        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper])
        mol = get_molecule()

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                                xml_ff_torsion_bo)
        topology = Topology.from_molecules([mol])

        omm_system, ret_top = forcefield.create_openmm_system(
                topology,
                charge_from_molecules=[mol],
                return_topology=True,
                toolkit_registry=toolkit_registry)

        off_torsion_force = [force for force in omm_system.getForces() if
                               isinstance(force, openmm.PeriodicTorsionForce)][0]

        ret_mol = list(ret_top.reference_molecules)[0]
        bond = ret_mol.get_bond_between(*central_atoms)

        # ambertools appears to yield around 1.0009 for this bond
        assert abs(bond.fractional_bond_order - 1.) > 1.e-6
        assert .95 < bond.fractional_bond_order < 1.05

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if (((atom2 == atom2_mol) and (atom3 == atom3_mol)) or
                    ((atom2 == atom3_mol) and (atom3 == atom2_mol))):
                k = params[-1]

                # do a hand calculation as a sanity check
                slope = (7.5312 - 4.184)/(1.8 - 1.0)
                k_interpolated_ret = slope * (bond.fractional_bond_order - 1.0) + 4.184
                assert_almost_equal(k_interpolated_ret, k/k.unit, 2)

                # check that we *are not* matching the values we'd get if we
                # had offered our molecules to `partial_bond_orders_from_molecules`
                with pytest.raises(AssertionError):
                    assert_almost_equal(k/k.unit, k_interpolated)

    @pytest.mark.skipif( not(OpenEyeToolkitWrapper.is_available()), reason='Test requires OE toolkit')
    @pytest.mark.parametrize(
            ('get_molecule', 'k_interpolated', 'central_atoms'),
            [(create_ethanol, 4.953856, (1, 2)), (create_reversed_ethanol, 4.953856, (7, 6))])
    def test_fractional_bondorder_calculated_openeye(self, get_molecule, k_interpolated, central_atoms):
        """Test that torsion barrier heights interpolated from fractional bond
        orders calculated with the OpenEye toolkit are assigned within our
        expectations.

        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        mol = get_molecule()

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml',
                                xml_ff_torsion_bo)
        topology = Topology.from_molecules([mol])

        omm_system, ret_top = forcefield.create_openmm_system(
                topology,
                charge_from_molecules=[mol],
                return_topology=True,
                toolkit_registry=toolkit_registry)

        off_torsion_force = [force for force in omm_system.getForces() if
                               isinstance(force, openmm.PeriodicTorsionForce)][0]

        ret_mol = list(ret_top.reference_molecules)[0]
        bond = ret_mol.get_bond_between(*central_atoms)

        # openeye toolkit appears to yield around .9945 for this bond
        assert abs(bond.fractional_bond_order - 1.) > 1.e-6
        assert .95 < bond.fractional_bond_order < 1.05

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if (((atom2 == atom2_mol) and (atom3 == atom3_mol)) or
                    ((atom2 == atom3_mol) and (atom3 == atom2_mol))):
                k = params[-1]

                # do a hand calculation as a sanity check
                slope = (7.5312 - 4.184)/(1.8 - 1.0)
                k_interpolated_ret = slope * (bond.fractional_bond_order - 1.0) + 4.184
                assert_almost_equal(k_interpolated_ret, k/k.unit, 2)

                # check that we *are not* matching the values we'd get if we
                # had offered our molecules to `partial_bond_orders_from_molecules`
                with pytest.raises(AssertionError):
                    assert_almost_equal(k/k.unit, k_interpolated)


class TestSmirnoffVersionConverter:

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                        reason='Test requires OE toolkit to read mol2 files')
    @pytest.mark.parametrize(('freesolv_id', 'forcefield_version', 'allow_undefined_stereo'),
                             generate_freesolv_parameters_assignment_cases())
    @pytest.mark.parametrize(('spec'), ['0_1', '0_2', '0_3'])
    def test_read_smirnoff_spec_freesolv(self, freesolv_id, forcefield_version, allow_undefined_stereo, spec):
        """
        Test reading an 0.2 smirnoff spec file, by reading an 0.1 spec representation of a set of parameters,
        and ensuring that it parameterizes molecules identically to the same FF in the most recent spec

        Regression test on parameters assignment based on the FreeSolv set used in the 0.1 paper.

        This, contrarily to the similar AlkEthOH test, checks also constraints
        and improper torsions.

        The AlkEthOH set, however, does not have impropers, which should be
        tested separately. Currently, test_freesolv_parameters_assignment
        does the job.
        """
        from openforcefield.tests.utils import get_freesolv_file_path, compare_system_parameters

        mol2_file_path, xml_file_path = get_freesolv_file_path(freesolv_id, forcefield_version)

        # Load molecules.
        molecule = Molecule.from_file(mol2_file_path, allow_undefined_stereo=allow_undefined_stereo)

        # Create OpenFF System with the current toolkit.
        #forcefield_file_path = 'test_forcefields/old/smirnoff99Frosst_1_0_8_spec_0_2.offxml'
        forcefield_file_path = f'test_forcefields/old/test_ff_{forcefield_version}_spec_{spec}.offxml'
        ff = ForceField(forcefield_file_path, 'test_forcefields/old/hbonds.offxml')
        ff_system = ff.create_openmm_system(molecule.to_topology())

        # Load OpenMM System created with the 0.1 version of the toolkit.
        from simtk import openmm
        with open(xml_file_path, 'r') as f:
            xml_system = openmm.XmlSerializer.deserialize(f.read())

        # Compare parameters. We ignore the improper folds as in 0.0.3 we
        # used a six-fold implementation while we now use a three-fold way.
        # TODO: Reactivate charge comparison once we'll be able to read them from file.
        compare_system_parameters(ff_system, xml_system,
                                  systems_labels=('current OpenFF', 'SMIRNOFF 0.0.4'),
                                  ignore_charges=True, ignore_improper_folds=True)




@pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
def test_electrostatics_options(self):
    """Test parameter assignment using smirnoff99Frosst on laromustine with various long-range electrostatics options.
    """
    molecules_file_path = get_data_file_path('molecules/laromustine_tripos.mol2')
    molecule = openforcefield.topology.Molecule.from_file(molecules_file_path)
    forcefield = ForceField([smirnoff99Frosst_offxml_file_path, charge_increment_offxml_file_path])
    for method in ['PME', 'reaction-field', 'Coulomb']:
        # Change electrostatics method
        forcefield.forces['Electrostatics'].method = method
        f = partial(check_system_creation_from_molecule, forcefield, molecule)
        f.description = 'Testing {} parameter assignment using molecule {}'.format(offxml_file_path, molecule.name)
        #yield f
    # TODO: Implement a similar test, where we compare OpenMM energy evals from an
    #       AMBER-parameterized system to OFF-parameterized systems

@pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
def test_charge_increment(self):
    """Test parameter assignment using smirnoff99Frosst on laromustine with ChargeIncrementModel.
    """
    molecules_file_path = get_data_file_path('molecules/laromustine_tripos.mol2')
    molecule = openforcefield.topology.Molecule.from_file(molecules_file_path)
    forcefield = ForceField(['test_forcefields/smirnoff99Frosst.offxml', 'chargeincrement-test'])
    check_system_creation_from_molecule(forcefield, molecule)
    # TODO: We can't implement a test for chargeincrement yet because we
    #       haven't settled on a SMIRNOFF spec for chargeincrementmodel


@pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
def test_create_system_molecules_parmatfrosst_gbsa(self):
    """Test creation of a System object from small molecules to test parm@frosst forcefield with GBSA support.
    """
    molecules_file_path = get_data_file_path('molecules/AlkEthOH_test_filt1_tripos.mol2')
    check_parameter_assignment(
        offxml_file_path='test_forcefields/Frosst_AlkEthOH_GBSA.offxml', molecules_file_path=molecules_file_path)
    # TODO: Figure out if we just want to check that energy is finite (this is what the original test did,
    #       or compare numerically to a reference system.


# TODO: test_get_new_parameterhandler
# TODO: test_get_existing_parameterhandler
# TODO: test_get_parameter
# TODO: test_add_parameter
# TODO: test_store_cosmetic_attribs
# TODO: test_write_cosmetic_attribs
# TODO: test_store_cosmetic_elements (eg. Author)
# TODO: test_write_cosmetic_elements (eg. Author)
# TODO: add_handler_with_incompatible_kwargs (for example different scale14 vals)
# TODO: test_invalid aromaticity_model
# TODO: test_invalid_file_version
# TODO: test_library_charges
# TODO: test_forcefield_to_dict (ensure that ParameterHandlers serialize without collisions
#     and header-level attribs include handler attribs as well as attached units,
#     note that header attribs are not ordered)
# TODO: test_create_gbsa
