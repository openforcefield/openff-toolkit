#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================
"""
Tests for the SMIRNOFF ForceField class

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

# TODO: Split this file into many test files, potentially distributing tests within each subpackage next to the classes they test

from functools import partial
from io import StringIO
import tempfile
from tempfile import TemporaryDirectory

from numpy.testing import assert_equal
import parmed
from simtk import unit, openmm
from simtk.openmm import app

# from openforcefield.utils import get_data_file_path  #, generateTopologyFromOEMol, read_molecules
# from openforcefield.utils import check_energy_is_finite, get_energy
# from openforcefield.tests.utils import get_packmol_pdb_file_path, get_monomer_mol2_file_path, compare_system_energies
from openforcefield.tests.utils import *
from openforcefield.typing.engines.smirnoff import *

# TODO: Fix these imports and unskip these tests
import pytest
pytestmark = pytest.mark.skip('tests in test_smirnoff are outdated and should be integrated with test_forcefield.py')

#=============================================================================================
# CONSTANTS
#=============================================================================================

# TODO: Move these to setup?
# These paths should only be found in the test data directories, so we need to use get_data_file_path()
AlkEthOH_offxml_filename = os.path.join('test_forcefields', 'Frosst_AlkEthOH.offxml')
AlkEthOH_molecules_filename = get_data_file_path('molecules/AlkEthOH_test_filt1_tripos.mol2')
MiniDrugBank_molecules_filename = get_data_file_path('molecules/MiniDrugBank_tripos.mol2')
#chargeincrement_offxml_filename = get_data_file_path('chargeincrement-test.offxml')
# TODO: Swtich all paths to use os.path.join
tip3p_molecule_filename = get_data_file_path(os.path.join('systems', 'monomers', 'tip3p_water.mol2'))

# This is the production form of smirnoff99Frosst that should be found in the data directories
smirnoff99Frosst_offxml_filename = os.path.join('test_forcefields', 'smirnoff99Frosst.offxml')
tip3p_offxml_filename = os.path.join('test_forcefields', 'tip3p.offxml')

# TODO: Add tests to compare RDKit and OpenEye derived forcefields to make sure they are the same

# TODO: Move forcefields for testing only to a special test data directory, separate from the data/ paths that are automatically searched and populated with production force fields
#
# SMIRNOFF ForceField XML definitions for testing purposes
#

# This is a test forcefield that is not meant for actual use.
# It just tests various capabilities.
ffxml_standard = u"""\
<!-- Header block (optional) -->
<Date>Date: May-September 2016</Date>
<Author>C. I. Bayly, OpenEye Scientific Software and David Mobley, UCI</Author>

<Bonds potential="harmonic" length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
   <Bond smirks="[#6X4:1]-[#6X4:2]" length="1.526" k="620.0" id="b0001" parent_id="b0001"/> <!-- CT-CT from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#6X4:1]-[#1:2]" length="1.090" k="680.0" id="b0002" parent_id="b0002"/> <!-- CT-H_ from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#8:1]~[#1:2]" length="1.410" k="640.0" id="b0003" parent_id="b0003"/> <!-- DEBUG O-H -->
   <Bond smirks="[#6X4:1]-[O&amp;X2&amp;H1:2]" length="1.410" k="640.0" id="b0004" parent_id="b0004"/> <!-- CT-OH from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#6X4:1]-[O&amp;X2&amp;H0:2]" length="1.370" k="640.0" id="b0005" parent_id="b0005"/> <!--CT-OS from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#8X2:1]-[#1:2]" length="0.960" k="1106.0" id="b0006" parent_id="b0003"/> <!-- OH-HO from frcmod.Frosst_AlkEthOH -->
    <Bond smirks="[#6X3:1]!#[#6X3:2]" id="b5" parent_id="b5" k_bondorder1="820.0" k_bondorder2="1098" length_bondorder1="1.45" length_bondorder2="1.35"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
   <Bond smirks="[#6X3:1]-[#1:2]" id="b27" parent_id="b27" k="734.0" length="1.080"/> <!-- Christopher Bayly from par99, Aug 2016 -->
</Bonds>

<Angles potential="harmonic" angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
   <Angle smirks="[a,A:1]-[#6X4:2]-[a,A:3]" angle="109.50" k="100.0" id="a0001" parent_id="a0001"/> <!-- consensus matches all X-Csp3-X -->
   <Angle smirks="[*:1]~[#8:2]~[*:3]" angle="113." k="100.0" id="a0001b" parent_id="a0001"/> <!-- consensus matches all X-O-X, such as in water -->
   <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.50" k="70.0" id="a0002" parent_id="a0001"/> <!-- H1-CT-H1 from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]" angle="109.50" k="80.0" id="a0003" parent_id="a0001"/> <!-- CT-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#8X2:1]-[#6X4:2]-[#8X2:3]" angle="109.50" k="140.0" id="a0004" parent_id="a0001"/> <!-- O_-CT-O_ from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#6X4:1]-[#8X2:2]-[#1:3]" angle="108.50" k="110.0" id="a0005" parent_id="a0005"/> <!-- CT-OH-HO from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#6X4:1]-[#8X2:2]-[#6X4:3]" angle="109.50" k="120.0" id="a0006" parent_id="a0006"/> <!-- CT-OS-CT from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[*:1]~[#6X3:2]~[*:3]" angle="120." id="a10" parent_id="a10" k="140.0"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
   <Angle smirks="[#1:1]-[#6X3:2]~[*:3]" angle="120." id="a11" parent_id="a11" k="100.0"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
   <Angle smirks="[#1:1]-[#6X3:2]-[#1:3]" angle="120." id="a12" parent_id="a12" k="70.0"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
</Angles>

<ProperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
   <Proper smirks="[a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]" idivf1="9" periodicity1="3" phase1="0.0" k1="1.40" id="t0001" parent_id="t0001"/> <!-- X -CT-CT-X from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[a,A:1]-[#6X4:2]-[#8X2:3]-[#1:4]" idivf1="3" periodicity1="3" phase1="0.0" k1="0.50" id="t0002" parent_id="t0002"/> <!--X -CT-OH-X from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]" idivf1="3" periodicity1="3" phase1="0.0" k1="1.15" id="t0003" parent_id="t0003"/> <!-- X -CT-OS-X from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.15" id="t0004" parent_id="t0001"/> <!-- HC-CT-CT-HC from frcmod.Frosst_AlkEthOH; note discrepancy with parm@frosst which applies this ONLY to HC-CT-CT-HC and other hydrogens get the generic X -CT-CT-X -->
   <Proper smirks="[#1:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.16" id="t0005" parent_id="t0001"/> <!-- HC-CT-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#8X2:3]-[#1:4]" idivf1='1' periodicity1="3" phase1="0.0" k1="0.16" idivf2="1" periodicity2="1" phase2="0.0" k2="0.25" id="t0006" parent_id="t0002"/> <!-- HO-OH-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.18" idivf2="1" periodicity2="2" phase2="180.0" k2="0.25" idivf3="1" periodicity3="1" phase3="180.0" k3="0.20" id="t0007" parent_id="t0001"/> <!-- CT-CT-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.383" idivf2="1" periodicity2="2" phase2="180.0" k2="0.1" id="t0008" parent_id="t0003"/> <!-- CT-CT-OS-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#8X2:2]-[#6X4:3]-[O&amp;X2&amp;H0:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.10" idivf2="1" periodicity2="2" phase2="180.0" k2="0.85" idivf3="1" periodicity3="1" phase3="180.0" k3="1.35" id="t0009" parent_id="t0003"/> <!-- CT-OS-CT-OS from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#8X2:1]-[#6X4:2]-[#6X4:3]-[#8X2:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.144" idivf2="1" periodicity2="2" phase2="0.0" k2="1.175" id="t0010" parent_id="t0001"/> <!-- O_-CT-CT-O_ from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#8X2:1]-[#6X4:2]-[#6X4:3]-[#1:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.0" idivf2="1" periodicity2="1" phase2="0.0" k2="0.25" id="t0011" parent_id="t0001"/> <!-- H_-CT-CT-O_ from frcmod.Frosst_AlkEthOH; discrepancy with parm@frosst with H2,H3-CT-CT-O_ per C Bayly -->
   <Proper smirks="[#1:1]-[#6X4:2]-[#6X4:3]-[OX2:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.0" idivf2="1" periodicity2="1" phase2="0.0" k2="0.25" id="t0012" parent_id="t0001"/> <!-- HC,H1-CT-CT-OH,OS from frcmod.Frosst_AlkEthOH ; note corresponding H2 and H3 params missing so they presumably get generic parameters (another parm@frosst bug) -->
   <Proper smirks="[*:1]~[#6X3:2]-[#6X4:3]~[*:4]" id="t13" idivf1="1" k1="0.000" periodicity1="3" phase1="0.0"/> <!-- parm99 generic, Bayly, Aug 2016 -->
   <Proper smirks="[#1:1]-[#6X4:2]-[#6X3:3]=[#6X3:4]" id="t15" parent_id="t15" idivf1="1" k1="0.380" periodicity1="3" phase1="180.0" phase2="0.0" k2="1.150" periodicity2="1" idivf2="1"/> <!-- parm99, Bayly, Aug 2016 -->
   <Proper smirks="[*:1]~[#6X3:2]-[#6X3:3]~[*:4]" id="t18" parent_id="t18" idivf1="1" k1="1." periodicity1="2" phase1="180.0"/> <!-- parm99 generic, Bayly, Aug 2016 -->
   <Proper smirks="[*:1]~[#6X3:2]:[#6X3:3]~[*:4]" id="t20" parent_id="t20" idivf1="1" k1="3.625" periodicity1="2" phase1="180.0"/> <!-- parm99 generic, Bayly, Aug 2016 -->
   <Proper smirks="[*:1]-[#6X3:2]=[#6X3:3]-[*:4]" id="t21" parent_id="t21" idivf1="1" k1="6." periodicity1="2" phase1="180.0"/> <!-- parm99 generic, Bayly, Aug 2016 -->
</ProperTorsions>

<ImproperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
   <Improper smirks="[a,A:1]~[#6X3:2]([a,A:3])~[OX1:4]" periodicity1="2" phase1="180.0" k1="10.5" id="i0001" parent_id="i0001"/> <!-- X -X -C -O  from frcmod.Frosst_AlkEthOH; none in set but here as format placeholder -->
</ProperTorsions>

<vdW potential="Lennard-Jones-12-6" combining_rules="Loentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch="8.0" cutoff="9.0" long_range_dispersion="isotropic">
   <!-- sigma is in angstroms, epsilon is in kcal/mol -->
   <Atom smirks="[#1:1]" rmin_half="1.4870" epsilon="0.0157" id="n0001" parent_id="n0001"/> <!-- making HC the generic hydrogen -->
   <Atom smirks="[$([#1]-C):1]" rmin_half="1.4870" epsilon="0.0157" id="n0002" parent_id="n0001"/> <!-- HC from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[$([#1]-C-[#7,#8,F,#16,Cl,Br]):1]" rmin_half="1.3870" epsilon="0.0157" id="n0003" parent_id="n0002"/> <!-- H1 from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[$([#1]-C(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]):1]" rmin_half="1.2870" epsilon="0.0157" id="n0004" parent_id="n0003"/> <!--H2 from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[$([#1]-C(-[#7,#8,F,#16,Cl,Br])(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]):1]" rmin_half="1.1870" epsilon="0.0157" id="n0005" parent_id="n0004"/> <!--H3 from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[#1$(*-[#8]):1]" rmin_half="0.0000" epsilon="0.0000" id="n0006" parent_id="n0001"/> <!-- HO from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[#6:1]" rmin_half="1.9080" epsilon="0.1094" id="n0007" parent_id="n0007"/> <!-- making CT the generic carbon -->
   <Atom smirks="[#6X4:1]" rmin_half="1.9080" epsilon="0.1094" id="n0008" parent_id="n0007"/> <!-- CT from frcmod.Frosst_AlkEthOH-->
   <Atom smirks="[#8:1]" rmin_half="1.6837" epsilon="0.1700" id="n0009" parent_id="n0009"/> <!-- making OS the generic oxygen -->
   <Atom smirks="[#8X2:1]" rmin_half="1.6837" epsilon="0.1700" id="n0010" parent_id="n0009"/> <!-- OS from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[#8X2+0$(*-[#1]):1]" rmin_half="1.7210" epsilon="0.2104" id="n0011" parent_id="n0009"/> <!-- OH from frcmod.Frosst_AlkEthOH -->
</vdW>
"""

ffxml_chargeincrement = u"""\
  <ChargeIncrementModel number_of_conformers="10" quantum_chemical_method="AM1" partial_charge_method="CM2" increment_unit="elementary_charge">
    <!-- A fractional charge can be moved along a single bond -->
    <ChargeIncrement smirks="[#6X4:1]-[#6X3a:2]" charge1increment="-0.0073" charge2increment="+0.0073"/>
    <ChargeIncrement smirks="[#6X4:1]-[#6X3a:2]-[#7]" charge1increment="+0.0943" charge2increment="-0.0943"/>
    <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge1increment="-0.0718" charge2increment="+0.0718"/>
    <!--- Alternatively, factional charges can be redistributed among any number of bonded atoms -->
  <ChargeIncrement smirks="[N:1](H:2)(H:3)" charge1increment="+0.02" charge2increment="-0.01" charge3increment="-0.01"/>
</ChargeIncrementModel>
"""

ffxml_constraints = u"""\
<Constraints distance_unit="angstroms">
  <!-- constrain all bonds to hydrogen to their equilibrium bond length -->
  <Constraint smirks="[#1:1]-[*:2]" id="constraint0001"/>
</Constraints>
"""

ffxml_gbsa = u"""\
<GBSA gb_model="OBC1" solvent_dielectric="78.5" solute_dielectric="1" radius_units="nanometers" sa_model="ACE" surface_area_penalty="5.4*calories/mole/angstroms**2" solvent_radius="1.4*angstroms">
  <Atom smirks="[#1:1]" radius="0.12" scale="0.85" id="gb0001"/>
  <Atom smirks="[#6:1]" radius="0.22" scale="0.72" id="gb0002"/>
  <Atom smirks="[#7:1]" radius="0.155" scale="0.79" id="gb0003"/>
  <Atom smirks="[#8:1]" radius="0.15" scale="0.85" id="gb0004"/>
  <Atom smirks="[#9:1]" radius="0.15" scale="0.88" id="gb0005"/>
  <Atom smirks="[#14:1]" radius="0.21" scale="0.8" id="gb0006"/>
  <Atom smirks="[#15:1]" radius="0.185" scale="0.86" id="gb0007"/>
  <Atom smirks="[#16:1]" radius="0.18" scale="0.96" id="gb0008"/>
  <Atom smirks="[#17:1]" radius="0.17" scale="0.8" id="gb0009"/>
</GBSA>
"""

ffxml_contents_gbsa = u"""\
<?xml version="1.0"?>

<SMIRNOFF use_fractional_bondorder="True">

%(ffxml_standard)s
%(ffxml_constraints)s
%(ffxml_gbsa)s

</SMIRNOFF>
""" % globals()

ffxml_contents_chargeincrement = u"""\
<?xml version="1.0"?>

%(ffxml_standard)s
%(ffxml_constraints)s
%(ffxml_chargeincrement)s

</SMIRNOFF>
""" % globals()

ffxml_contents_noconstraints = u"""\
<?xml version="1.0"?>

<SMIRNOFF use_fractional_bondorder="True">

%(ffxml_standard)s

</SMIRNOFF>
""" % globals()

ffxml_contents = u"""\
<?xml version="1.0"?>

<SMIRNOFF use_fractional_bondorder="True">

%(ffxml_standard)s
%(ffxml_constraints)s

</SMIRNOFF>
""" % globals()

# Set up another, super minimal FF to test that alternate aromaticity
# models implemented properly
ffxml_MDL_contents = u"""\
<?xml version="1.0"?>

<SMIRNOFF version="0.1" aromaticity_model="OEAroModel_MDL">

!-- SMIRKS (SMIRKS Force Field) minimal file, not intended for use -->
  <Date>Date: Sept. 7, 2016</Date>
  <Author>D. L. Mobley, UC Irvine</Author>
  <Bonds potential="harmonic" length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[*:1]~[*:2]" id="b1" k="2000.0" length="4.0"/>
    <Bond smirks="[#6X3:1]-[#8X2:2]" id="b16" k="960.0" length="1.240"/>
  </Bonds>
  <Angles potential="harmonic" angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
    <Angle smirks="[*:1]~[*:2]~[*:3]" angle="109.5" id="a1" k="160.0"/>
  </Angles>
  <ProperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Proper smirks="[*:1]~[*:2]~[*:3]~[*:4]" id="t1" idivf1="4" k1="3.50" periodicity1="2" phase1="180.0"/>
  </ProperTorsions>
  <ImproperTorsions potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Improper smirks="[a,A:1]~[#6X3:2]([a,A:3])~[OX1:4]" periodicity1="2" phase1="180.0" k1="10.5" id="i0001" parent_id="i0001"/> <!-- X -X -C -O  from frcmod.Frosst_AlkEthOH; none in set but here as format placeholder -->
  </ImproperTorsions>
  <vdW potential="Lennard-Jones-12-6" combining_rules="Loentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole" switch="8.0" cutoff="9.0" long_range_dispersion="isotropic">
    <Atom smirks="[*:1]" epsilon="0.2100" id="n1" rmin_half="1.6612"/>
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n2" rmin_half="0.6000"/>
  </vdW>

</SMIRNOFF>
"""

#=============================================================================================
# TESTS
#=============================================================================================

#
# Test various components
#


class TestApplyForceField:
    """Test the use of ForceField to parameterize Topology objects.
    """

    def check_system_creation_from_molecule(self, forcefield, molecule):
        """
        Generate a System from the given OEMol and SMIRNOFF forcefield and check that its energy is finite.

        Parameters
        ----------
        forcefield : openforcefield.typing.engines.smirnoff.ForceField
            SMIRNOFF forcefield object
        molecule : openforcefield.topology.Molecule
            Molecule to test (which must have coordinates)

        """
        topology = Topology.from_molecules(molecule)
        system = forcefield.create_system(topology)
        check_energy_is_finite(system, molecule.positions)
        # TODO: Add check_energy_is_finite to test_FF.test_paramerize_ethanol

    def check_system_creation_from_topology(forcefield, topology, positions):
        """
        Generate a System from the given topology, OEMols matching the contents of the topology, and SMIRNOFF forcefield and check that its energy is finite.

        Parameters
        ----------
        forcefield : ForceField
            SMIRNOFF forcefield
        topology : openforcefield.topology.Topology
            Topology for which System should be constructed
        positions : simtk.unit.Quantity with dimension (natoms,3) with units of length
            Positions of all particles

        """
        system = forcefield.create_system(topology)
        check_energy_is_finite(system, positions)

    def check_parameter_assignment(self,
                                   offxml_filename=smirnoff99Frosst_offxml_filename,
                                   molecules_filename=AlkEthOH_molecules_filename):
        """Test parameter assignment using specified forcefield on all molecules from specified source.

        Parameters
        ----------
        offxml_filename : str, optional, default='test_forcefields/smirnoff99Frosst.offxml'
            Filename from which SMIRNOFF .offxml XML file is to be loaded.
        molecules_filename : str, optional, default='molecules/AlkEthOH_test_filt1_tripos.mol2
            Filename from which molecule identities and positions are to be loaded.
        """
        forcefield = ForceField(offxml_filename)
        for molecule in openforcefield.topology.Molecule.from_file(molecules_filename):
            f = partial(check_system_creation_from_molecule, forcefield, molecule)
            f.description = 'Testing {} parameter assignment using molecule {}'.format(offxml_filename, molecule.name)
            yield f


    def test_partial_bondorder(verbose=False):
        """Test energies of a molecule which activates partial bond order code."""
        # Create benzene molecule
        molecule = Molecule.from_iupac('benzene')
        topology = molecule.to_topology()

        # Load forcefield
        ff = ForceField(ffxml_contents_noconstraints)

        # Set up once using AM1BCC charges
        # TODO: Use OECharges_AM1BCCSym charges
        system = ff.create_system(topology)

        # Check that energy is what it ought to be -- the partial bond order
        # for benzene makes the energy a bit higher than it would be without it
        energy = get_energy(system, positions)
        if energy < 7.50 or energy > 7.60:
            raise Exception(
                "Partial bond order code seems to have issues, as energy for benzene is outside of tolerance in tests."
            )

        # Set up once also without asking for charges
        # TODO: Don't create charges
        system = ff.create_system(topology)
        energy = get_energy(system, positions)
        # Energy is lower with user supplied charges (which in this case are zero)
        if energy < 4.00 or energy > 6.0:
            raise Exception(
                "Partial bond order code seems to have issues when run with user-provided charges, as energy for benzene is out of tolerance in tests."
            )

    def test_improper(verbose=False):
        """Test implement of impropers on benzene."""
        # Create benzene molecule
        molecule = Molecule.from_iupac('benzene')
        topology = molecule.to_topology()

        # Load forcefield
        ff = ForceField('test_forcefields/benzene_minimal.offxml')

        # Load AMBER files and compare
        inpcrd = get_data_file_path('molecules/benzene.crd')
        prmtop = get_data_file_path('molecules/benzene.top')
        g0, g1, e0, e1 = compare_amber_smirnoff(prmtop, inpcrd, ff, mol, skip_assert=True)

        # Check that torsional energies the same to 1 in 10^6
        rel_error = np.abs((g0['torsion'] - g1['torsion']) / g0['torsion'])
        if rel_error > 6e-3:  #Note that this will not be tiny because we use three-fold impropers and they use a single improper
            raise Exception(
                "Improper torsion energy for benzene differs too much (relative error %.4g) between AMBER and SMIRNOFF."
                % rel_error)


class TestForceFieldLabeling:
    """Tests for labeling parameter types in molecules
    """

    def test_label_molecules(self):
        """Test labeling/getting stats on labeling molecules"""
        molecules = read_molecules(get_data_file_path('molecules/AlkEthOH_test_filt1_tripos.mol2'), verbose=verbose)
        ffxml = get_data_file_path('test_forcefields/Frosst_AlkEthOH.offxml')
        get_molecule_parameterIDs(molecules, ffxml)

    def test_molecule_labeling(self):
        """Test using labelMolecules to see which parameters applied to a molecule
        """
        from openforcefield.topology.testsystems import SMILESTopology
        topology = SMILESTopology('CCC')
        forcefield = ForceField('test_forcefields/Frosst_AlkEthOH.offxml')
        labels = forcefield.label_molecules(topology)[0]

        # Check that force terms aren't empty
        for forcename in ['Bonds', 'Angles', 'ProperTorsions', 'ImproperTorsions', 'vdW', 'Electrostatics']:
            if not forcename in labels:
                raise Exception("No force term assigned for {}.".format(forcename))


class TestExceptionHandling:
    """Tests for exception handling for incomplete parameterization
    """

    def test_parameter_completeness_check(self):
        """Test that proper exceptions are raised if a force field fails to assign parameters to valence terms in a molecule.
        """
        from openforcefield.topology.testsystems import SMILESTopology
        topology = SMILESTopology('CCC')

        # Test vdW error checking by wiping out required parameter
        parameter_edits = [
            ('vdW', '[#136:1]'),
            ('Bonds', '[#136:1]~[*:2]'),
            ('Angles', '[#136:1]~[*:2]~[*:3]'),
            ('ProperTorsions', '[#136:1]~[*:2]~[*:3]~[*:4]'),
        ]
        for (tag, smirks) in parameter_edits:
            forcefield = ForceField('test_forcefields/Frosst_AlkEthOH.offxml')
            forcefield.forces[tag].parameters[0].smirks = smirks
            with self.assertRaises(Exception):
                system = forcefield.create_system(topology)


#
# TODO: Finish updating the tests below this line
#


# TODO: Is there any way to make this test not depend on OpenEye?
@pytest.mark.skip
#@requires_openeye_licenses
def test_improper_pyramidal(verbose=False):
    """Test implement of impropers on ammonia."""
    from openeye import oeomega
    from openeye import oequacpac
    from oeommtools.utils import oemol_to_openmmTop, openmmTop_to_oemol
    from openforcefield.utils import get_data_file_path, extractPositionsFromOEMol
    # Prepare ammonia
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, 'N')
    oechem.OEAddExplicitHydrogens(mol)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(100)
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(False)
    # Generate charges
    chargeEngine = oequacpac.OEAM1BCCCharges()
    status = omega(mol)
    oequacpac.OEAssignCharges(mol, chargeEngine)
    # Assign atom names
    oechem.OETriposAtomTypes(mol)
    oechem.OETriposAtomNames(mol)
    # Set up minimization
    ff = ForceField('test_forcefields/ammonia_minimal.offxml')
    topology, positions = oemol_to_openmmTop(mol)
    system = ff.createSystem(topology, [mol], verbose=verbose)
    positions = extractPositionsFromOEMol(mol)
    integrator = openmm.VerletIntegrator(2.0 * unit.femtoseconds)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    # Minimize energy
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    newpositions = state.getPositions()
    outmol = openmmTop_to_oemol(topology, state.getPositions())
    # Sum angles around the N atom
    for atom in outmol.GetAtoms(oechem.OEIsInvertibleNitrogen()):
        aidx = atom.GetIdx()
        nbors = list(atom.GetAtoms())
        ang1 = math.degrees(oechem.OEGetAngle(outmol, nbors[0], atom, nbors[1]))
        ang2 = math.degrees(oechem.OEGetAngle(outmol, nbors[1], atom, nbors[2]))
        ang3 = math.degrees(oechem.OEGetAngle(outmol, nbors[2], atom, nbors[0]))
    ang_sum = math.fsum([ang1, ang2, ang3])
    # Check that sum of angles around N is within 1 degree of 353.5
    if abs(ang_sum - 353.5) > 1.0:
        raise Exception(
            "Improper torsion for ammonia differs too much from reference partly pyramidal geometry."
            "The reference case has a sum of H-N-H angles (3x) at 353.5 degrees; the test case has a sum of {}.".
            format(ang_sum))


def test_MDL_aromaticity(verbose=False):
    """Test support for alternate aromaticity models."""
    ffxml = StringIO(ffxml_MDL_contents)
    ff = ForceField(ffxml)
    from openeye import oechem
    mol = oechem.OEMol()
    oechem.OEParseSmiles(mol, 'c12c(cccc1)occ2')
    oechem.OEAddExplicitHydrogens(mol)

    labels = ff.labelMolecules([mol], verbose=True)
    # The bond 6-7 should get the b16 parameter iff the MDL model is working, otherwise it will pick up just the generic
    details = labels[0]['HarmonicBondGenerator']
    found = False
    for (atom_indices, pid, smirks) in details:
        if pid == 'b16' and atom_indices == [6, 7]:
            found = True
    if not found: raise Exception("Didn't find right param.")


def test_change_parameters(verbose=False):
    """Test modification of forcefield parameters."""
    from openeye import oechem
    # Load simple OEMol
    ifs = oechem.oemolistream(get_data_file_path('molecules/AlkEthOH_c100.mol2'))
    mol = oechem.OEMol()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol)
    oechem.OETriposAtomNames(mol)

    # Load forcefield file
    ff = ForceField('test_forcefields/Frosst_AlkEthOH.offxml')

    topology = generateTopologyFromOEMol(mol)
    # Create initial system
    system = ff.createSystem(topology, [mol], verbose=verbose)
    # Get initial energy before parameter modification
    positions = positions_from_oemol(mol)
    old_energy = get_energy(system, positions)

    # Get params for an angle
    params = ff.getParameter(smirks='[a,A:1]-[#6X4:2]-[a,A:3]')
    # Modify params
    params['k'] = '0.0'
    ff.setParameter(params, smirks='[a,A:1]-[#6X4:2]-[a,A:3]')
    # Write params
    ff.writeFile(tempfile.TemporaryFile(suffix='.offxml'))
    # Make sure params changed energy! (Test whether they get rebuilt on system creation)
    system = ff.createSystem(topology, [mol], verbose=verbose)
    energy = get_energy(system, positions)
    if verbose:
        print("New energy/old energy:", energy, old_energy)
    if np.abs(energy - old_energy) < 0.1:
        raise Exception("Error: Parameter modification did not change energy.")


class TestForceFieldExport:
    """Test the ability to export ForceField-parameterized systems to other major simulation packages.
    """

    @staticmethod
    def get_mixture_tesystem(packmol_boxname, monomers):
        """Retrieve test system information for specified packmol mixed-component box and monomers.

        Parameters
        ----------
        packmol_boxname : str
            Prefix of PDB file to be retrieved from systems/packmolboxes in test data directory
        monomers : list of str
            List of monomer names within packmol box

        Returns
        -------
        topology : openforcefield.topology.Topology
            Topology for the packmol box
        positions : simtk.unit.Quantity with dimension (nparticles,3) with units compatible with angstroms
            Positions corresponding to particles in ``topology``
        """
        pdbfile = app.PDBFile(get_packmol_pdb_file_path(packmol_boxname))
        molecules = [Molecule.from_file(get_monomer_mol2_file_path(name)) for name in monomers]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        positions = pdbfile.positions
        return topology, positions

    def test_amber_roundtrip(self):
        """Save a System (a mixture) to AMBER, read back in, verify yields same energy and force terms."""
        forcefield = ForceField(AlkEthOH_offxml_filename)
        topology, positions = get_mixture_testsystem('cyclohexane_ethanol_0.4_0.6', ['ethanol', 'cyclohexane'])
        system = forcefield.create_system(pdbfile.topology, mols)
        openmm_topology = topology.to_openmm()

        # TODO: Can we provide a way to more directly export to AMBER?
        from openforcefield.utils import save_system_to_amber
        with TemporaryDirectory() as tmpdir:
            # Write AMBER prmtop, inpcrd
            prmtop_filename = os.path.join(tmpdir, 'output.prmtop')
            inpcrd_filename = os.path.join(tmpdir, 'output.inpcrd')
            save_system_to_amber(openmm_topology, system, positions, prmtop_filename, inpcrd_filename)

            # Read back in and cross-check energies
            structure = parmed.load_file(prmtop, inpcrd)  # TODO: Is this a structure or prmtop object in ParmEd?
            # TODO: Make sure the SMIRNOFF systems is also created with no cutoff, no constraints, and no implicit solvent
            amber_system = structure.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=None)
            groups0, groups1, energy0, energy1 = compare_system_energies(
                openmm_topology, openmm_topology, amber_system, system, positions, verbose=False)

    def test_gromacs_roundtrip(self):
        """Save a System (a mixture) to gromacs, read back in, verify yields same energy and force terms."""
        forcefield = ForceField(AlkEthOH_offxml_filename)
        topology, positions = get_mixture_testsystem('cyclohexane_ethanol_0.4_0.6', ['ethanol', 'cyclohexane'])
        system = forcefield.create_system(pdbfile.topology, mols)
        openmm_topology = topology.to_openmm()

        # TODO: Can we provide a way to more directly export to gromacs?
        from openforcefield.utils import save_system_to_amber
        with TemporaryDirectory() as tmpdir:
            # Write gromacs .gro and .top files
            gro_filename = os.path.join(tmpdir, 'output.gro')
            top_filename = os.path.join(tmpdir, 'output.top')
            save_system_to_gromacs(openmm_topology, system, positions, top_filename, gro_filename)

            # Read back in and cross-check energies
            top = parmed.load_file(top_filename)
            gro = parmed.load_file(gro_filename)
            # TODO: Make sure the SMIRNOFF systems is also created with no cutoff, no constraints, and no implicit solvent
            gromacs_system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=None)
            groups0, groups1, energy0, energy1 = compare_system_energies(
                openmm_topology, openmm_topology, gromacs_system, system, positions, verbose=False)


class TestSolventSupport:
    """Test support for rigid solvents like TIP3P.
    """

    # TODO: Move this to test_forcefield::TestForceFieldConstraints when it will be possible to specify TIP3P electrostatics.
    def test_tip3p_constraints(self):
        """Test that TIP3P distance costraints are correctly applied."""
        # Specify TIP3P constraint distances
        tip3p_oh_distance = 0.9572  # angstrom
        tip3p_hoh_angle = 104.52  # angle
        tip3p_hh_distance = tip3p_oh_distance * np.sin(np.radians(tip3p_hoh_angle / 2)) * 2
        expected_distances = [tip3p_oh_distance, tip3p_oh_distance, tip3p_hh_distance]

        # Load TIP3P molecule
        molecule = Molecule.from_file(tip3p_molecule_filename)

        # Extract topology and positions
        topology = Topology.from_molecules(molecule)
        positions = molecule.positions

        # Create OpenMM System
        ff = ForceField(tip3p_offxml_filename)
        system = ff.create_system(topology)

        # TODO: This is probably unnecessary and we can simply check that the System has all the Constraints in their position.
        # Run dynamics.
        integrator = openmm.VerletIntegrator(2.0 * unit.femtoseconds)
        context = openmm.Context(system, integrator)
        context.setPositions(positions)
        integrator.step(50)

        # Ensure constrained distances are correct
        state = context.getState(getPositions=True)
        new_positions = state.getPositions(asNumpy=True) / unit.angstroms
        distances = []
        for atom_1, atom_2 in [(0, 1), (0, 2), (1, 2)]:  # pair of atoms O-H1, O-H2, H1-H2
            distances.append(np.linalg.norm(new_positions[atom_1] - new_positions[atom_2]))
        err_msg = 'expected distances [O-H1, O-H2, H1-H2]: {} A, new distances: {} A'
        assert np.allclose(expected_distances, distances), err_msg.format(expected_distances, distances)

    # TODO: Move this to test_forcefield::TestForceFieldParameterAssignment when it will be possible to specify TIP3P electrostatics.
    def test_tip3p_solvated_molecule_energy():
        """Check the energy of a TIP3P solvated molecule is the same with SMIRNOFF and OpenMM.

        This test makes also use of defining the force field with multiple FFXML files,
        and test the energy of mixture systems (i.e. molecule plus water).

        """

        # TODO remove this function when ParmEd#868 is fixed
        def add_water_bonds(system):
            """Hack to make tip3p waters work with ParmEd Structure.createSystem.

            Bond parameters need to be specified even when constrained.
            """
            k_tip3p = 462750.4 * unit.kilojoule_per_mole / unit.nanometers**2
            length_tip3p = 0.09572 * unit.nanometers
            for force in system.getForces():
                if isinstance(force, openmm.HarmonicBondForce):
                    force.addBond(particle1=0, particle2=1, length=length_tip3p, k=k_tip3p)
                    force.addBond(particle1=0, particle2=2, length=length_tip3p, k=k_tip3p)

        monomers_dir = get_data_file_path(os.path.join('systems', 'monomers'))

        # Load TIP3P molecule
        tip3p_molecule = Molecule.from_file(tip3p_molecule_filename)
        tip3p_positions = tip3p_molecule.positions

        # Extract topology and positions
        top3p_topology = Topology.from_molecules(molecule)

        # Create OpenMM System
        tip3p_forcefield = ForceField(tip3p_offxml_filename)
        tip3p_openmm_system = tip3p_forcefield.create_system(tip3p_topology)

        # TODO remove this line when ParmEd#868 is fixed
        # Hack to make tip3p waters work with ParmEd Structure.createSystem.
        add_water_bonds(tip3p_openmm_system)

        tip3p_openmm_structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), tip3p_openmm_system,
                                                                       tip3p_positions)

        smirnoff_forcefield = ForceField(smirnoff99Frosst_offxml_filename)

        for monomer in ['ethanol', 'methane']:
            # Load molecule
            other_molecule_filepath = os.path.join(monomers_dir, monomer + '.mol2')
            other_molecule = Molecule.from_file(other_molecule_filepath)

            # Create ParmEd Structure
            other_molecule_structure = smirnoff_forcefield.create_structure(other_molecule.to_topology(),
                                                                            other_molecule.positions)

            # Merge molecule and OpenMM TIP3P water
            combined_structure = tip3p_openmm_structure + other_molecule_structure
            #structure.positions = np.append(tip3p_molecule.positions, molecule.positions, axis=0) # TODO: This doesn't seem necessary since positions are correctly concatenated already
            system = combined_structure.createSystem()
            energy_openmm = get_energy(system, combined_structure.positions)

            # Test the creation of system with multiple force field files.
            ff = ForceField(smirnoff_offxml_filename, tip3p_offxml_filename)
            combined_topology = Topology.from_molecules([tip3p_molecule, molecule])
            system = ff.create_system(topology)

            # TODO remove this line when ParmEd#868 is fixed
            # Add bonds to match ParmEd hack above.
            add_water_bonds(system)

            energy_smirnoff = get_energy(system, structure.positions)

            # The energies of the systems must be the same.
            assert np.isclose(energy_openmm, energy_smirnoff), 'OpenMM: {}, SMIRNOFF: {}'.format(
                energy_openmm, energy_smirnoff)
