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

from functools import partial
from unittest import TestCase
from openforcefield import utils
from simtk import unit
import numpy as np
from numpy.testing import assert_almost_equal

import pytest
from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry

from openforcefield.utils import get_data_filename

from openforcefield.topology.molecule import Molecule
from openforcefield.typing.engines.smirnoff import ForceField


#=============================================================================================
# Data
#=============================================================================================

simple_xml_ff = str.encode('''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="1.0" aromaticity_model="OEAroModel_MDL">
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
  <vdW potential="Lennard-Jones-12-6" combining_rules="Loentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch="8.0*angstroms" cutoff="9.0*angstroms" long_range_dispersion="isotropic">
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157" id="n2" rmin_half="1.4870"/>
  </vdW>
</SMIRNOFF>
''')

xml_ff_w_comments = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="1.0" aromaticity_model="OEAroModel_MDL">
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
  <vdW potential="Lennard-Jones-12-6" combining_rules="Loentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch="8.0*angstroms" cutoff="9.0*angstroms" long_range_dispersion="isotropic">
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157" id="n2" rmin_half="1.4870"/>
  </vdW>
</SMIRNOFF>
'''

xml_ff_w_cosmetic_elements = '''<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="1.0" aromaticity_model="OEAroModel_MDL">
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
  <vdW potential="Lennard-Jones-12-6" combining_rules="Loentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" rmin_half_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch="8.0*angstroms" cutoff="9.0*angstroms" long_range_dispersion="isotropic">
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n1" rmin_half="0.6000"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157" id="n2" rmin_half="1.4870"/>
  </vdW>
</SMIRNOFF>
'''


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


#=============================================================================================
# TESTS
#=============================================================================================

class TestForceField(TestCase):
    """Test the ForceField class"""



    def test_create_forcefield_from_file(self):
        """Test empty constructor"""
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        assert len(forcefield._parameter_handlers['Bonds']._parameters) == 87
        assert len(forcefield._parameter_handlers['Angles']._parameters) == 38
        assert len(forcefield._parameter_handlers['ProperTorsions']._parameters) == 158
        assert len(forcefield._parameter_handlers['ImproperTorsions']._parameters) == 4
        assert len(forcefield._parameter_handlers['vdW']._parameters) == 35

    def test_create_forcefield_from_xml_string(self):
        forcefield = ForceField(simple_xml_ff)
        assert len(forcefield._parameter_handlers['Bonds']._parameters) == 2
        assert len(forcefield._parameter_handlers['Angles']._parameters) == 2
        assert len(forcefield._parameter_handlers['ProperTorsions']._parameters) == 2
        assert len(forcefield._parameter_handlers['ImproperTorsions']._parameters) == 2
        assert len(forcefield._parameter_handlers['vdW']._parameters) == 2

    # TODO: Support writing out offxml
    @pytest.mark.skip
    def test_xml_string_roundtrip(self):
        forcefield = ForceField(simple_xml_ff)
        data = forcefield._parameter_io_handlers['XML'].to_string()
        raise Exception(data)

    def test_parameterize_ethanol(self):
        from simtk.openmm import app
        from openforcefield.topology import Topology
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        pdbfile = app.PDBFile(get_data_filename('systems/packmol_boxes/1_ethanol.pdb'))
        #toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [ Molecule.from_file(get_data_filename(name)) for name in ('molecules/ethanol.mol2',) ]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology)

    def test_parameterize_1_cyclohexane_1_ethanol(self):
        from simtk.openmm import app
        from openforcefield.topology import Topology
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        # pdbfile = app.PDBFile(get_data_filename('systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb'))
        pdbfile = app.PDBFile(get_data_filename('systems/packmol_boxes/1_cyclohexane_1_ethanol.pdb'))
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [Molecule.from_file(get_data_filename(name)) for name in ('molecules/ethanol.mol2',
                                                                              'molecules/cyclohexane.mol2')]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology)

    # This test takes too long with the initial implementation of the toolkit
    @pytest.mark.skip
    def test_parameterize_large_system(self):
        from simtk.openmm import app
        from openforcefield.topology import Topology
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        pdbfile = app.PDBFile(get_data_filename('systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb'))
        molecules = [Molecule.from_file(get_data_filename(name)) for name in ('molecules/ethanol.mol2',
                                                                              'molecules/cyclohexane.mol2')]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = forcefield.create_openmm_system(topology)

    def test_parameterize_different_reference_ordering(self):
        """
        Test parameterizing the same PDB, using reference mol2s that have different atom orderings.
        The results of both should be identical.
        """
        from simtk.openmm import app
        from openforcefield.topology import Topology
        from simtk.openmm import XmlSerializer
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        pdbfile = app.PDBFile(get_data_filename('systems/packmol_boxes/1_ethanol.pdb'))
        # Load the unique molecules with one atom ordering
        molecules1 = [Molecule.from_file(get_data_filename('molecules/ethanol.mol2'))]
        topology1 = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules1)
        omm_system1 = forcefield.create_openmm_system(topology1)
        # Load the unique molecules with a different atom ordering
        molecules2 = [Molecule.from_file(get_data_filename('molecules/ethanol_reordered.mol2'))]
        topology2 = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules2)
        omm_system2 = forcefield.create_openmm_system(topology2)

        serialized_1 = XmlSerializer.serialize(omm_system1)
        serialized_2 = XmlSerializer.serialize(omm_system2)

        serialized_1 = round_charge(serialized_1)
        serialized_2 = round_charge(serialized_2)

        assert serialized_1 == serialized_2


# from_filename
# from xml_string
# from_xml_bytes
# from_url
# get_new_parameterhandler
# get_existing_parameterhandler
# get_parameter
# add_parameter
# add_parameter_fractional_bondorder
# create_force_fractional_bondorder
# store_cosmetic_attribs
# write_cosmetic_attribs
# store_cosmetic_elements
# write_cosmetic_elements
# add_handler_with_incompatible_kwargs (for example different scale14 vals)
# invalid aromaticity_model
# invalid_file_version
# library_charges
