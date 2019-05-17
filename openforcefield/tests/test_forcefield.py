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

import pytest
from tempfile import NamedTemporaryFile

from openforcefield.utils.toolkits import OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry
from openforcefield.utils import get_data_file_path
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField, IncompatibleParameterError, SMIRNOFFSpecError
from openforcefield.typing.engines.smirnoff import XMLParameterIOHandler


#======================================================================
# GLOBAL CONSTANTS
#======================================================================

# File paths.
TIP3P_SDF_FILE_PATH = get_data_file_path(os.path.join('systems', 'monomers', 'water.sdf'))

GENERIC_BOND_LENGTH = 1.09 * unit.angstrom
XML_FF_GENERICS = f"""<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.2" aromaticity_model="OEAroModel_MDL">
  <Bonds length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[*:1]~[*:2]" id="b1" k="680.0" length="{GENERIC_BOND_LENGTH/unit.angstrom}"/>
  </Bonds>
  <Angles angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
    <Angle smirks="[*:1]~[*:2]~[*:3]" angle="109.5" id="a1" k="100.0"/>
  </Angles>
  <ProperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Proper smirks="[*:1]~[*:2]~[*:3]-[*:4]" id="t1" idivf1="1" k1="0.156" periodicity1="3" phase1="0.0"/>
  </ProperTorsions>
  <ImproperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Improper smirks="[*:1]~[*:2](~[*:3])~[*:4]" id="i1" k1="1.1" periodicity1="2" phase1="180."/>
  </ImproperTorsions>
  <vdW potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch_width="1.0" switch_width_unit="angstrom" cutoff="9.0" cutoff_unit="angstrom" method="cutoff">
    <Atom smirks="[*:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
  </vdW>
  <Electrostatics method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0" cutoff_unit="angstrom"/>
  <ToolkitAM1BCC/>
</SMIRNOFF>
"""

simple_xml_ff = str.encode('''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.2" aromaticity_model="OEAroModel_MDL">
  <Bonds length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0" length="1.526"/>
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0" length="1.51"/>
  </Bonds>
  <Angles angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5" id="a1" k="100.0"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5" id="a2" k="70.0"/>
  </Angles>
  <ProperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156" periodicity1="3" phase1="0.0"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180" periodicity1="3" phase1="0.0" periodicity2="2" phase2="180.0" idivf2="1" k2="0.250" periodicity3="1" phase3="180.0" idivf3="1" k3="0.200"/>
  </ProperTorsions>
  <ImproperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1" periodicity1="2" phase1="180."/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5" periodicity1="2" phase1="180."/>
  </ImproperTorsions>
  <vdW potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch_width="1.0" switch_width_unit="angstrom" cutoff="9.0" cutoff_unit="angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157" id="n2" rmin_half="1.4870"/>
  </vdW>
  <Electrostatics method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0" cutoff_unit="angstrom"/>
  <ToolkitAM1BCC/>
</SMIRNOFF>
''')

xml_ff_w_comments = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.2" aromaticity_model="OEAroModel_MDL">
  <!-- SMIRNOFF (SMIRKS Native Open Force Field) template file -->
  <Date>2018-07-14</Date>
  <Author>C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, UC Irvine; D. L. Mobley, UC Irvine</Author>
  <!-- This file is meant for processing via openforcefield.typing.engines.smirnoff -->
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the bond energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Bonds length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0" length="1.526" />
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0" length="1.51"/>
  </Bonds>
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the angle energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Angles angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5" id="a1" k="100.0"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5" id="a2" k="70.0"/>
  </Angles>
  <ProperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156" periodicity1="3" phase1="0.0"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180" periodicity1="3" phase1="0.0" periodicity2="2" phase2="180.0" idivf2="1" k2="0.250" periodicity3="1" phase3="180.0" idivf3="1" k3="0.200"/>
  </ProperTorsions>
  <ImproperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1" periodicity1="2" phase1="180."/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5" periodicity1="2" phase1="180."/>
  </ImproperTorsions>
  <vdW potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch_width_unit="angstrom" cutoff="9.0" cutoff_unit="angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157" id="n2" rmin_half="1.4870"/>
  </vdW>
  <Electrostatics method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0" cutoff_unit="angstrom" pme_tolerance="0.00001"/>
  <ToolkitAM1BCC/>
</SMIRNOFF>
'''

xml_ff_w_cosmetic_elements = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.2" aromaticity_model="OEAroModel_MDL">
  <!-- SMIRNOFF (SMIRKS Native Open Force Field) template file -->
  <Date>2018-07-14</Date>
  <Author>C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, UC Irvine; D. L. Mobley, UC Irvine</Author>
  <!-- This file is meant for processing via openforcefield.typing.engines.smirnoff -->
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the bond energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Bonds length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0" length="1.526" parameters="k, length" parameterize_eval="blah=blah2" />
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0" length="1.51"/>
  </Bonds>
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the angle energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Angles angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2" cosmetic_element="why not?">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5" id="a1" k="100.0"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5" id="a2" k="70.0"/>
  </Angles>
  <ProperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156" periodicity1="3" phase1="0.0"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180" periodicity1="3" phase1="0.0" periodicity2="2" phase2="180.0" idivf2="1" k2="0.250" periodicity3="1" phase3="180.0" idivf3="1" k3="0.200"/>
  </ProperTorsions>
  <ImproperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1" periodicity1="2" phase1="180."/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5" periodicity1="2" phase1="180."/>
  </ImproperTorsions>
  <vdW potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch_width="8.0" switch_width_unit="angstrom" cutoff="9.0" cutoff_unit="angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157" id="n2" rmin_half="1.4870"/>
  </vdW>
  <Electrostatics method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0" cutoff_unit="angstrom" pme_tolerance="0.00001"/>
  <ToolkitAM1BCC/>
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
        chunksp[0] = str('%.2d' % (float(chunksp[0])))
        chunk = '" sig'.join(chunksp)
        xmlsp[index] = chunk
    return ' q="'.join(xmlsp)

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
    ethanol.add_bond(0, 1, 1, False)  # C0 - C1
    ethanol.add_bond(1, 2, 1, False)  # C1 - O2
    ethanol.add_bond(0, 3, 1, False)  # C0 - H3
    ethanol.add_bond(0, 4, 1, False)  # C0 - H4
    ethanol.add_bond(0, 5, 1, False)  # C0 - H5
    ethanol.add_bond(1, 6, 1, False)  # C1 - H6
    ethanol.add_bond(1, 7, 1, False)  # C1 - H7
    ethanol.add_bond(2, 8, 1, False)  # O2 - H8
    charges = unit.Quantity(np.array([-0.4, -0.3, -0.2, -0.1, 0.01, 0.1, 0.2, 0.3, 0.4]), unit.elementary_charge)
    ethanol.partial_charges = charges
    return ethanol

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
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'reaction-field', 'has_periodic_box': False,
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'PME', 'has_periodic_box': True,
     'omm_force': openmm.NonbondedForce.PME, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'cutoff', 'electrostatics_method': 'PME', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},

    {'vdw_method': 'PME', 'electrostatics_method': 'Coulomb', 'has_periodic_box': True,
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'Coulomb', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'reaction-field', 'has_periodic_box': True,
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'reaction-field', 'has_periodic_box': False,
     'omm_force': None, 'exception': IncompatibleParameterError, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'PME', 'has_periodic_box': True,
     'omm_force': openmm.NonbondedForce.LJPME, 'exception': None, 'exception_match': ''},
    {'vdw_method': 'PME', 'electrostatics_method': 'PME', 'has_periodic_box': False,
     'omm_force': openmm.NonbondedForce.NoCutoff, 'exception': None, 'exception_match': ''},
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

    def test_create_forcefield_from_file(self):
        """Test empty constructor"""
        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        assert len(forcefield._parameter_handlers['Bonds']._parameters) == 87
        assert len(forcefield._parameter_handlers['Angles']._parameters) == 38
        assert len(forcefield._parameter_handlers['ProperTorsions']._parameters) == 158
        assert len(forcefield._parameter_handlers['ImproperTorsions']._parameters) == 4
        assert len(forcefield._parameter_handlers['vdW']._parameters) == 35

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_file_list(self):
        # These offxml files are located in package data path, which is automatically installed and searched
        filenames = [smirnoff99Frosst_offxml_filename, tip3p_offxml_filename]
        # Create a forcefield from multiple offxml files
        forcefield = ForceField(filenames)

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_create_forcefield_from_filename_iterator(self):
        # These offxml files are located in package data path, which is automatically installed and searched
        filenames = [smirnoff99Frosst_offxml_filename, tip3p_offxml_filename]
        # A generator should work as well
        forcefield = ForceField(iter(filenames))

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
        assert len(forcefield._parameter_handlers['ProperTorsions']._parameters) == 2
        assert len(forcefield._parameter_handlers['ImproperTorsions']._parameters) == 2
        assert len(forcefield._parameter_handlers['vdW']._parameters) == 2

    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    def test_deep_copy(self):
        forcefield = ForceField(smirnoff99Frosst_offxml_filename)
        # Deep copy
        forcefield2 = copy.deepcopy(cls.forcefield)
        assert_forcefields_equal(cls.forcefield, forcefield2,
                                 "ForceField deep copy does not match original ForceField")


    @pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
    # TODO: This should check the output of forcefield.to_dict
    def test_serialize(self):

        forcefield = ForceField(smirnoff99Frosst_offxml_filename)
        # Serialize/deserialize
        serialized_forcefield = cls.forcefield.__getstate__()
        forcefield2 = ForceField.__setstate__(serialized_forcefield)
        assert_forcefields_equal(cls.forcefield, forcefield2,
                                 "Deserialized serialized ForceField does not match original ForceField")

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
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg {'parameters': 'k, length'} passed") as excinfo:
            forcefield = ForceField(xml_ff_w_cosmetic_elements)

        # Create a forcefield from XML successfully, by explicitly permitting cosmetic attributes
        forcefield_1 = ForceField(xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)

        # Convert the forcefield back to XML
        string_1 = forcefield_1.to_string('XML', discard_cosmetic_attributes=False)

        # Ensure that the new XML string has cosmetic attributes in it
        assert 'cosmetic_element="why not?"' in string_1
        assert 'parameterize_eval="blah=blah2"' in string_1
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg {'parameters': 'k, length'} passed") as excinfo:
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

    @pytest.mark.parametrize('filename_extension', ['xml', 'XML', 'offxml', 'OFFXML'])
    @pytest.mark.parametrize('specified_format', [None, 'xml', 'XML', '.xml', '.XML',
                                                  'offxml', 'OFFXML', '.offxml', '.OFFXML',
                                                  XMLParameterIOHandler()])
    def test_xml_file_roundtrip(self, filename_extension, specified_format):
        """
        Test roundtripping a ForceField to and from an XML file
        """
        # These files will be deleted once garbage collection runs (end of this function)
        iofile1 = NamedTemporaryFile(suffix='.' + filename_extension)
        iofile2 = NamedTemporaryFile(suffix='.' + filename_extension)
        forcefield_1 = ForceField(simple_xml_ff)
        forcefield_1.to_file(iofile1.name, io_format=specified_format)
        forcefield_2 = ForceField(iofile1.name)
        forcefield_2.to_file(iofile2.name, io_format=specified_format)
        assert open(iofile1.name).read() == open(iofile2.name).read()


    @pytest.mark.parametrize('filename_extension', ['xml', 'XML', 'offxml', 'OFFXML'])
    @pytest.mark.parametrize('specified_format', [None, 'xml', 'XML', '.xml', '.XML',
                                                  'offxml', 'OFFXML', '.offxml', '.OFFXML',
                                                  XMLParameterIOHandler()])
    def test_xml_file_roundtrip_keep_cosmetic(self, filename_extension, specified_format):
        """
        Test roundtripping a forcefield to an XML file with and without retaining cosmetic elements
        """
        # These files will be deleted once garbage collection runs (end of this function)
        iofile1 = NamedTemporaryFile(suffix='.' + filename_extension)
        iofile2 = NamedTemporaryFile(suffix='.' + filename_extension)
        iofile3 = NamedTemporaryFile(suffix='.' + filename_extension)

        # Ensure an exception is raised if we try to read the XML string with cosmetic attributes
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg {'parameters': 'k, length'} passed") as excinfo:
            forcefield = ForceField(xml_ff_w_cosmetic_elements)

        # Create a forcefield from XML successfully
        forcefield_1 = ForceField(xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)

        # Convert the forcefield back to XML, keeping cosmetic attributes
        forcefield_1.to_file(iofile1.name, discard_cosmetic_attributes=False, io_format=specified_format)

        # Ensure that the new XML string has cosmetic attributes in it
        assert 'cosmetic_element="why not?"' in open(iofile1.name).read()
        assert 'parameterize_eval="blah=blah2"' in open(iofile1.name).read()
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg {'parameters': 'k, length'} passed") as excinfo:
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



    def test_load_two_sources(self):
        """Test loading data from two SMIRNOFF data sources"""
        ff = ForceField(simple_xml_ff, xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True)
        assert len(ff.get_parameter_handler('Bonds').parameters) == 4

    def test_load_two_sources_incompatible_tags(self):
        """Test loading data from two SMIRNOFF data sources which have incompatible physics"""
        # Make an XML forcefield with a modifiedvdW 1-4 scaling factor
        nonstandard_xml_ff = xml_ff_w_comments.replace('scale14="0.5"', 'scale14="1.0"')
        with pytest.raises(IncompatibleParameterError, match="handler value: 0.5, incompatible value: 1.0") as excinfo:
            ff = ForceField(simple_xml_ff, nonstandard_xml_ff)

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_ethanol(self, toolkit_registry, registry_description):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        molecules = []
        molecules.append(Molecule.from_smiles('CCO'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology, toolkit_registry=toolkit_registry)

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_1_cyclohexane_1_ethanol(self, toolkit_registry, registry_description):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_cyclohexane_1_ethanol.pdb'))
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = []
        molecules.append(Molecule.from_smiles('CCO'))
        molecules.append(Molecule.from_smiles('C1CCCCC1'))
        # molecules = [Molecule.from_file(get_data_file_path(name)) for name in ('molecules/ethanol.mol2',
        #                                                                      'molecules/cyclohexane.mol2')]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology)

    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_1_cyclohexane_1_ethanol_vacuum(self, toolkit_registry, registry_description):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_cyclohexane_1_ethanol.pdb'))
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = []
        molecules.append(Molecule.from_smiles('CCO'))
        molecules.append(Molecule.from_smiles('C1CCCCC1'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        topology.box_vectors = None

        omm_system = forcefield.create_openmm_system(topology)



    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_parameterize_no_matching_reference(self, toolkit_registry, registry_description):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_cyclohexane_1_ethanol.pdb'))
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = []
        molecules.append(Molecule.from_smiles('CC'))
        with pytest.raises(ValueError, match='No match found for molecule'):
            topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

    @pytest.mark.slow
    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    @pytest.mark.parametrize("box", ['ethanol_water.pdb',
                                     'cyclohexane_water.pdb',
                                     'cyclohexane_ethanol_0.4_0.6.pdb',
                                     'propane_methane_butanol_0.2_0.3_0.5.pdb'])
    def test_parameterize_large_system(self, toolkit_registry, registry_description, box):
        from simtk.openmm import app

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')
        box_filename = get_data_file_path(os.path.join('systems', 'packmol_boxes', box))
        pdbfile = app.PDBFile(box_filename)
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
    def test_charges_from_molecule(self, toolkit_registry, registry_description):
        """Test skipping charge generation and instead getting charges from the original Molecule"""
        # Create an ethanol molecule without using a toolkit
        molecules = [create_ethanol()]

        from simtk.openmm import app, NonbondedForce

        filename = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        omm_system = forcefield.create_openmm_system(topology, charge_from_molecules=molecules,
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
        topology2 = Topology.from_openmm(pdbfile2.topology, unique_molecules=molecules)

        omm_system2 = forcefield.create_openmm_system(topology2, charge_from_molecules=molecules,
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
    def test_some_charges_from_molecule(self, toolkit_registry, registry_description):
        """
        Test creating an OpenMM system where some charges come from a Molecule, but others come from toolkit
        calculation
        """
        ethanol = create_ethanol()
        cyclohexane = create_cyclohexane()
        molecules = [ethanol, cyclohexane]

        from simtk.openmm import app, NonbondedForce

        filename = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
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



    @pytest.mark.parametrize("toolkit_registry,registry_description", toolkit_registries)
    def test_pass_invalid_kwarg_to_create_openmm_system(self, toolkit_registry, registry_description):
        """Test to ensure an exception is raised when an unrecognized kwarg is passed """
        from simtk.openmm import app

        filename = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
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
        forcefield.get_parameter_handler('vdW', {})._method = vdw_method
        forcefield.get_parameter_handler('Electrostatics', {})._method = electrostatics_method

        pdbfile = app.PDBFile(get_data_file_path('systems/test_systems/1_ethanol.pdb'))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        if not(has_periodic_box):
            topology.box_vectors = None

        if exception is None:
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
                omm_system = forcefield.create_openmm_system(topology)


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
                                        bond_length=GENERIC_BOND_LENGTH)


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
        forcefield_file_path = 'test_forcefields/old/smirnoff99Frosst_' + forcefield_version + '.offxml'
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
    molecules_filename = get_data_file_path('molecules/laromustine_tripos.mol2')
    molecule = openforcefield.topology.Molecule.from_file(molecules_filename)
    forcefield = ForceField([smirnoff99Frosst_offxml_filename, chargeincrement_offxml_filename])
    for method in ['PME', 'reaction-field', 'Coulomb']:
        # Change electrostatics method
        forcefield.forces['Electrostatics'].method = method
        f = partial(check_system_creation_from_molecule, forcefield, molecule)
        f.description = 'Testing {} parameter assignment using molecule {}'.format(offxml_filename, molecule.name)
        #yield f
    # TODO: Implement a similar test, where we compare OpenMM energy evals from an
    #       AMBER-parameterized system to OFF-parameterized systems

@pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
def test_chargeincrement(self):
    """Test parameter assignment using smirnoff99Frosst on laromustine with ChargeIncrementModel.
    """
    molecules_filename = get_data_file_path('molecules/laromustine_tripos.mol2')
    molecule = openforcefield.topology.Molecule.from_file(molecules_filename)
    forcefield = ForceField(['test_forcefields/smirnoff99Frosst.offxml', 'chargeincrement-test'])
    check_system_creation_from_molecule(forcefield, molecule)
    # TODO: We can't implement a test for chargeincrement yet because we
    #       haven't settled on a SMIRNOFF spec for chargeincrementmodel


@pytest.mark.skip(reason='Needs to be updated for 0.2.0 syntax')
def test_create_system_molecules_parmatfrosst_gbsa(self):
    """Test creation of a System object from small molecules to test parm@frosst forcefield with GBSA support.
    """
    molecules_filename = get_data_file_path('molecules/AlkEthOH_test_filt1_tripos.mol2')
    check_parameter_assignment(
        offxml_filename='Frosst_AlkEthOH_GBSA.offxml', molecules_filename=molecules_filename)
    # TODO: Figure out if we just want to check that energy is finite (this is what the original test did,
    #       or compare numerically to a reference system.

# TODO: test_get_new_parameterhandler
# TODO: test_get_existing_parameterhandler
# TODO: test_get_parameter
# TODO: test_add_parameter
# TODO: test_add_parameter_fractional_bondorder
# TODO: test_create_force_fractional_bondorder
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
