from functools import partial
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, generateTopologyFromOEMol, read_molecules
from openeye import oechem
from simtk.openmm import app
from simtk import unit, openmm
import numpy as np
from io import StringIO
import tempfile
import parmed
import os

# This is a test forcefield that is not meant for actual use.
# It just tests various capabilities.
ffxml_standard = u"""\
<!-- Header block (optional) -->
<Date>Date: May-September 2016</Date>
<Author>C. I. Bayly, OpenEye Scientific Software and David Mobley, UCI</Author>

<HarmonicBondForce length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2" fractional_bondorder="interpolate-linear">
   <Bond smirks="[#6X4:1]-[#6X4:2]" length="1.526" k="620.0" id="b0001" parent_id="b0001"/> <!-- CT-CT from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#6X4:1]-[#1:2]" length="1.090" k="680.0" id="b0002" parent_id="b0002"/> <!-- CT-H_ from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#8:1]~[#1:2]" length="1.410" k="640.0" id="b0003" parent_id="b0003"/> <!-- DEBUG O-H -->
   <Bond smirks="[#6X4:1]-[O&amp;X2&amp;H1:2]" length="1.410" k="640.0" id="b0004" parent_id="b0004"/> <!-- CT-OH from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#6X4:1]-[O&amp;X2&amp;H0:2]" length="1.370" k="640.0" id="b0005" parent_id="b0005"/> <!--CT-OS from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#8X2:1]-[#1:2]" length="0.960" k="1106.0" id="b0006" parent_id="b0003"/> <!-- OH-HO from frcmod.Frosst_AlkEthOH -->
    <Bond smirks="[#6X3:1]!#[#6X3:2]" id="b5" parent_id="b5" k_bondorder1="820.0" k_bondorder2="1098" length_bondorder1="1.45" length_bondorder2="1.35"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
   <Bond smirks="[#6X3:1]-[#1:2]" id="b27" parent_id="b27" k="734.0" length="1.080"/> <!-- Christopher Bayly from par99, Aug 2016 -->
</HarmonicBondForce>

<HarmonicAngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
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
</HarmonicAngleForce>

<PeriodicTorsionForce phase_unit="degrees" k_unit="kilocalories_per_mole">
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
   <Improper smirks="[a,A:1]~[#6X3:2]([a,A:3])~[OX1:4]" periodicity1="2" phase1="180.0" k1="10.5" id="i0001" parent_id="i0001"/> <!-- X -X -C -O  from frcmod.Frosst_AlkEthOH; none in set but here as format placeholder -->
</PeriodicTorsionForce>

<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole">
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
</NonbondedForce>

<BondChargeCorrections method="AM1" increment_unit="elementary_charge">
  <BondChargeCorrection smirks="[#6X4:1]-[#6X3a:2]" increment="+0.0073" id="c0001" parent_id="c0001"/> <!-- tetrahedral carbon bonded to aromatic carbon correction -->
  <BondChargeCorrection smirks="[#6X4:1]-[#6X3a:2]-[#7]" increment="-0.0943" id="c0002" parent_id="c0002"/> <!-- tetrahedral carbon bonded to aromatic carbon (bonded to a nitrogen) -->
  <BondChargeCorrection smirks="[#6X4:1]-[#8:2]" increment="+0.0718" id="c0003" parent_id="c0003"/> <!-- tetrahedral carbon bonded to an oxygen :    13   11   31    1   0.0718 -->
</BondChargeCorrections>
"""

ffxml_constraints = u"""\
<Constraints distance_unit="angstroms">
  <!-- constrain all bonds to hydrogen to their equilibrium bond length -->
  <Constraint smirks="[#1:1]-[*:2]" id="constraint0001"/>
</Constraints>
"""

ffxml_gbsa = u"""\
<GBSAForce gb_model="OBC1" solvent_dielectric="78.5" solute_dielectric="1" radius_units="nanometers" sa_model="ACE" surface_area_penalty="5.4*calories/mole/angstroms**2" solvent_radius="1.4*angstroms">
  <Atom smirks="[#1:1]" radius="0.12" scale="0.85" id="gb0001"/>
  <Atom smirks="[#6:1]" radius="0.22" scale="0.72" id="gb0002"/>
  <Atom smirks="[#7:1]" radius="0.155" scale="0.79" id="gb0003"/>
  <Atom smirks="[#8:1]" radius="0.15" scale="0.85" id="gb0004"/>
  <Atom smirks="[#9:1]" radius="0.15" scale="0.88" id="gb0005"/>
  <Atom smirks="[#14:1]" radius="0.21" scale="0.8" id="gb0006"/>
  <Atom smirks="[#15:1]" radius="0.185" scale="0.86" id="gb0007"/>
  <Atom smirks="[#16:1]" radius="0.18" scale="0.96" id="gb0008"/>
  <Atom smirks="[#17:1]" radius="0.17" scale="0.8" id="gb0009"/>
</GBSAForce>
"""

ffxml_contents_gbsa = u"""\
<?xml version="1.0"?>

<SMIRNOFF use_fractional_bondorder="True">

%(ffxml_standard)s
%(ffxml_constraints)s
%(ffxml_gbsa)s

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
  <HarmonicBondForce length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[*:1]~[*:2]" id="b1" k="2000.0" length="4.0"/>
    <Bond smirks="[#6X3:1]-[#8X2:2]" id="b16" k="960.0" length="1.240"/>
  </HarmonicBondForce>
  <HarmonicAngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
    <Angle smirks="[*:1]~[*:2]~[*:3]" angle="109.5" id="a1" k="160.0"/>
</HarmonicAngleForce>
<PeriodicTorsionForce phase_unit="degrees" k_unit="kilocalories_per_mole">
<Improper smirks="[a,A:1]~[#6X3:2]([a,A:3])~[OX1:4]" periodicity1="2" phase1="180.0" k1="10.5" id="i0001" parent_id="i0001"/> <!-- X -X -C -O  from frcmod.Frosst_AlkEthOH; none in set but here as format placeholder -->
    <Proper smirks="[*:1]~[*:2]~[*:3]~[*:4]" id="t1" idivf1="4" k1="3.50" periodicity1="2" phase1="180.0"/>
  </PeriodicTorsionForce>
  <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole">
    <Atom smirks="[*:1]" epsilon="0.2100" id="n1" rmin_half="1.6612"/>
    <Atom smirks="[#1:1]" epsilon="0.0157" id="n2" rmin_half="0.6000"/>
  </NonbondedForce>
</SMIRNOFF>
"""

def positions_from_oemol(mol):
    """
    Extract OpenMM positions from OEMol.

    Parameters
    ----------
    mol : oechem.openeye.OEMol
        OpenEye molecule from which to extract coordinates.

    Returns
    -------
    positions : simtk.unit.Quantity of dimension (nparticles,3)

    """
    from openeye import oeomega
    if mol.GetDimension() != 3:
        # Assign coordinates
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(1)
        omega.SetIncludeInput(False)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
        status = omega(mol)  # generate conformation

    coordinates = mol.GetCoords()
    natoms = len(coordinates)
    positions = np.zeros([natoms,3], np.float32)
    for index in range(natoms):
        (x,y,z) = coordinates[index]
        positions[index,0] = x
        positions[index,1] = y
        positions[index,2] = z
    positions = unit.Quantity(positions, unit.angstroms)
    return positions

def check_energy_is_finite(system, positions):
    """
    Check the potential energy is not NaN.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use

    """
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    if np.isnan(energy):
        raise Exception('Potential energy is NaN')

def get_energy(system, positions):
    """
    Return the potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    return energy

def test_read_ffxml():
    """Test reading of ffxml files.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))

def test_read_ffxml_gbsa():
    """Test reading of ffxml files with GBSA support.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH_GBSA.ffxml'))

def check_system_creation_from_molecule(forcefield, mol, chargeMethod=None, verbose=False):
    """
    Generate a System from the given OEMol and SMIRNOFF forcefield and check that its energy is finite.

    Parameters
    ----------
    forcefield : ForceField
        SMIRNOFF forcefield
    mol : oechem.OEMol
        Molecule to test (need not have coordinates)
    chargeMethod : str, optional, default=None
        Charge method to use in creating system

    """

    topology = generateTopologyFromOEMol(mol)
    system = forcefield.createSystem(topology, [mol], chargeMethod=chargeMethod, verbose=verbose)
    # Test energy computation.
    positions = positions_from_oemol(mol)
    check_energy_is_finite(system, positions)

def check_system_creation_from_topology(forcefield, topology, mols, positions, chargeMethod=None, verbose=False):
    """
    Generate a System from the given topology, OEMols matching the contents of the topology, and SMIRNOFF forcefield and check that its energy is finite.

    Parameters
    ----------
    forcefield : ForceField
        SMIRNOFF forcefield
    topology : simtk.openmm.app.Topology
        Topology to construct system from
    mols : list of oechem.OEMol
        Reference molecules
    positions : simtk.unit.Quantity with dimension (natoms,3) with units of length
    chargeMethod : str, optional, default=None
        Charge method to use in creating system

    """
    from openforcefield.typing.engines.smirnoff import CutoffPeriodic
    system = forcefield.createSystem(topology, mols, verbose=verbose, chargeMethod=chargeMethod, nonbondedMethod=CutoffPeriodic)
    # Test energy computation.
    check_energy_is_finite(system, positions)

def check_AlkEthOH(forcefield, description="", chargeMethod=None, verbose=False):
    """Test creation of System from AlkEthOH small molecules.
    """
    from openeye import oechem
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEthOH-tripos.mol2.gz'))
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        args = { 'verbose' : verbose, 'chargeMethod' : chargeMethod }
        f = partial(check_system_creation_from_molecule, forcefield, mol, **args)
        f.description ='Testing creation of system object from small molecules (%s) %s' % (mol.GetTitle(), description)
        yield f

def test_create_system_molecules_features(verbose=False):
    """Test creation of a System object from small molecules to test various ffxml features
    """
    ffxml = StringIO(ffxml_contents_gbsa)
    forcefield = ForceField(ffxml)

    for chargeMethod in [None, 'BCC', 'OECharges_AM1BCCSym']:
        for f in check_AlkEthOH(forcefield, description="to test ffxml features with charge method %s" % str(chargeMethod), chargeMethod=chargeMethod, verbose=verbose):
            yield f

def test_create_system_molecules_parmatfrosst(verbose=False):
    """Test creation of a System object from small molecules to test parm@frosst forcefield.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
    for f in check_AlkEthOH(forcefield, "to test Parm@Frosst parameters", verbose=verbose):
        yield f

def test_create_system_molecules_parmatfrosst_gbsa(verbose=False):
    """Test creation of a System object from small molecules to test parm@frosst forcefield with GBSA support.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH_GBSA.ffxml'))
    for f in check_AlkEthOH(forcefield, "to test Parm@Frosst parameters", verbose=verbose):
        yield f

def check_boxes(forcefield, description="", chargeMethod=None, verbose=False):
    """Test creation of System from boxes of mixed solvents.
    """
    # Read monomers
    mols = list()
    monomers = ['water', 'cyclohexane', 'ethanol', 'propane', 'methane', 'butanol']
    from openeye import oechem
    mol = oechem.OEGraphMol()
    for monomer in monomers:
        filename = get_data_filename(os.path.join('systems', 'monomers', monomer + '.sdf'))
        ifs = oechem.oemolistream(filename)
        while oechem.OEReadMolecule(ifs, mol):
            oechem.OETriposAtomNames(mol)
            mols.append( oechem.OEGraphMol(mol) )
    if verbose: print('%d reference molecules loaded' % len(mols))

    # Read systems.
    boxes = ['ethanol_water.pdb',  'cyclohexane_water.pdb',
        'cyclohexane_ethanol_0.4_0.6.pdb', 'propane_methane_butanol_0.2_0.3_0.5.pdb']
    from simtk.openmm.app import PDBFile
    for box in boxes:
        filename = get_data_filename(os.path.join('systems', 'packmol_boxes', box))
        pdbfile = PDBFile(filename)
        f = partial(check_system_creation_from_topology, forcefield, pdbfile.topology, mols, pdbfile.positions, chargeMethod=chargeMethod, verbose=verbose)
        f.description = 'Test creation of System object from %s %s' % (box, description)
        yield f

def test_create_system_boxes_features(verbose=False):
    """Test creation of a System object from some boxes of mixed solvents to test ffxml features.
    """
    ffxml = StringIO(ffxml_contents)
    forcefield = ForceField(ffxml)
    for chargeMethod in [None, 'BCC', 'OECharges_AM1BCCSym']:
        for f in check_boxes(forcefield, description="to test Parm@frosst parameters with charge method %s" % str(chargeMethod), chargeMethod=chargeMethod, verbose=verbose):
            yield f

def test_create_system_boxes_smirnoff99Frosst(verbose=False):
    """Test creation of a System object from some boxes of mixed solvents to test parm@frosst forcefield.
    """
    forcefield = ForceField(get_data_filename('forcefield/smirnoff99Frosst.ffxml'))
    for f in check_boxes(forcefield, description="to test Parm@frosst parameters", verbose=verbose):
        yield f

def test_smirnoff_energies_vs_parmatfrosst(verbose=False):
    """Test evaluation of energies from parm@frosst ffxml files versus energies of equivalent systems."""

    from openeye import oechem
    prefix = 'AlkEthOH_'
    molecules = [ 'r118', 'r12', 'c1161', 'r0', 'c100', 'c38', 'c1266' ]

    # Loop over molecules, load OEMols and prep for comparison/do comparison
    for molnm in molecules:
        f_prefix = os.path.join('molecules', prefix+molnm )
        mol2file = get_data_filename( f_prefix+'.mol2')
        prmtop = get_data_filename( f_prefix+'.top')
        crd = get_data_filename( f_prefix+'.crd')
        # Load special parm@frosst with parm99/parm@frosst bugs re-added for testing
        forcefield = ForceField( get_data_filename('forcefield/Frosst_AlkEthOH_parmAtFrosst.ffxml') )

        # Load OEMol
        mol = oechem.OEGraphMol()
        ifs = oechem.oemolistream(mol2file)
        flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
        ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
        oechem.OEReadMolecule(ifs, mol )
        oechem.OETriposAtomNames(mol)

        # Do comparison
        results = compare_molecule_energies( prmtop, crd, forcefield, mol, verbose = verbose )

def test_label_molecules(verbose=False):
    """Test labeling/getting stats on labeling molecules"""
    molecules = read_molecules(get_data_filename('molecules/AlkEthOH-tripos.mol2.gz'), verbose=verbose)
    ffxml = get_data_filename('forcefield/Frosst_AlkEthOH.ffxml')
    get_molecule_parameterIDs( molecules, ffxml)

def test_molecule_labeling(verbose = False):
    """Test using labelMolecules to see which parameters applied to an oemol."""
    from openeye import oechem
    mol = oechem.OEMol()
    oechem.OEParseSmiles(mol, 'CCC')
    oechem.OEAddExplicitHydrogens(mol)
    ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
    labels = ff.labelMolecules( [mol], verbose = verbose)

    # Check that force terms aren't empty
    print(labels[0].keys())
    if not 'HarmonicBondGenerator' in labels[0].keys():
        raise Exception("No force term assigned for harmonic bonds.")
    if not 'HarmonicAngleGenerator' in labels[0].keys():
        raise Exception("No force term assigned for harmonic angles.")
    if not 'PeriodicTorsionGenerator' in labels[0].keys():
        raise Exception("No force term assigned for periodic torsions.")
    if not 'NonbondedGenerator' in labels[0].keys():
        raise Exception("No nonbonded force term assigned.")

import unittest
class TestExceptionHandling(unittest.TestCase):
    def test_parameter_completeness_check(self):
        """Test that proper exceptions are raised if a force field fails to assign parameters to valence terms in a molecule."""
        from openeye import oechem
        mol = oechem.OEMol()
        oechem.OEParseSmiles(mol, 'CCC')
        oechem.OEAddExplicitHydrogens(mol)
        oechem.OETriposAtomNames(mol)
        ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
        topology = generateTopologyFromOEMol(mol)

        # Test nonbonded error checking by wiping out required LJ parameter
        params = ff.getParameter(paramID='n0001')
        params['smirks']='[#136:1]'
        ff.setParameter(paramID='n0001', params=params)
        ff.setParameter(paramID='n0002', params=params)
        with self.assertRaises(Exception):
            system = ff.createSystem( topology, [mol])
        ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))

        # Test bond error checking by wiping out a required bond parameter
        params = ff.getParameter(paramID='b0001')
        params['smirks'] = '[#136:1]~[*:2]'
        ff.setParameter( paramID='b0001', params=params)
        with self.assertRaises(Exception):
            system = ff.createSystem( topology, [mol])
        ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))

        # Test angle error checking by wiping out a required angle parameter
        params = ff.getParameter(paramID='a0001')
        params['smirks'] = '[#136:1]~[*:2]~[*:3]'
        ff.setParameter( paramID='a0001', params=params)
        with self.assertRaises(Exception):
            system = ff.createSystem( topology, [mol])
        ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))

        # Test torsion error checking by wiping out a required torsion parameter
        params = ff.getParameter(paramID='t0001')
        params['smirks'] = '[#136:1]~[*:2]~[*:3]~[*:4]'
        ff.setParameter( paramID='t0001', params=params)
        ff.setParameter( paramID='t0004', params=params)
        with self.assertRaises(Exception):
            system = ff.createSystem( topology, [mol])
        ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))



def test_partial_bondorder(verbose = False):
    """Test setup of a molecule which activates partial bond order code."""
    from openeye import oechem
    mol = oechem.OEMol()
    from openeye import oeiupac
    oeiupac.OEParseIUPACName(mol, 'benzene')
    positions = positions_from_oemol(mol)
    oechem.OETriposAtomNames(mol)
    topology = generateTopologyFromOEMol(mol)
    # Load forcefield from above
    ffxml = StringIO(ffxml_contents_noconstraints)
    ff = ForceField(ffxml)

    # Set up once using AM1BCC charges
    system = ff.createSystem(topology, [mol], chargeMethod = 'OECharges_AM1BCCSym', verbose = verbose)

    # Check that energy is what it ought to be -- the partial bond order
    # for benzene makes the energy a bit higher than it would be without it
    energy = get_energy(system, positions)
    if energy < 7.50 or energy > 7.60:
        raise Exception("Partial bond order code seems to have issues, as energy for benzene is outside of tolerance in tests.")

    # Set up once also without asking for charges
    system = ff.createSystem(topology, [mol], verbose = verbose)
    energy = get_energy(system, positions)
    # Energy is lower with user supplied charges (which in this case are zero)
    if energy < 4.00 or energy > 6.0:
        raise Exception("Partial bond order code seems to have issues when run with user-provided charges, as energy for benzene is out of tolerance in tests.")

def test_improper(verbose = False):
    """Test implement of impropers on benzene."""
    from openeye import oechem
    # Load benzene
    ifs = oechem.oemolistream(get_data_filename('molecules/benzene.mol2'))
    mol = oechem.OEMol()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol )
    ifs.close()
    # Load forcefield
    ffxml = get_data_filename('forcefield/benzene_minimal.ffxml')
    ff = ForceField(ffxml)

    # Load AMBER files and compare
    crd = get_data_filename('molecules/benzene.crd')
    top = get_data_filename('molecules/benzene.top')
    g0, g1, e0, e1 = compare_molecule_energies( top, crd, ff, mol, skip_assert = True)

    # Check that torsional energies the same to 1 in 10^6
    rel_error = np.abs(( g0['torsion']-g1['torsion'])/ g0['torsion'])
    if rel_error > 2e-5: #Note that this will not be tiny because we use six-fold impropers and they use a single improper
        raise Exception("Improper torsion energy for benzene differs too much (relative error %.4g) between AMBER and SMIRNOFF." % rel_error )


def test_MDL_aromaticity(verbose=False):
    """Test support for alternate aromaticity models."""
    ffxml = StringIO(ffxml_MDL_contents)
    ff = ForceField(ffxml)
    from openeye import oechem
    mol = oechem.OEMol()
    oechem.OEParseSmiles(mol, 'c12c(cccc1)occ2')
    oechem.OEAddExplicitHydrogens(mol)

    labels=ff.labelMolecules( [mol], verbose = True)
    # The bond 6-7 should get the b16 parameter iff the MDL model is working, otherwise it will pick up just the generic
    details = labels[0]['HarmonicBondGenerator']
    found = False
    for (atom_indices, pid, smirks) in details:
        if pid == 'b16' and atom_indices==[6,7]:
            found = True
    if not found: raise Exception("Didn't find right param.")


def test_change_parameters(verbose=False):
    """Test modification of forcefield parameters."""
    from openeye import oechem
    # Load simple OEMol
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEthOH_c100.mol2'))
    mol = oechem.OEMol()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol )
    oechem.OETriposAtomNames(mol)

    # Load forcefield file
    ffxml = get_data_filename('forcefield/Frosst_AlkEthOH.ffxml')
    ff = ForceField(ffxml)

    topology = generateTopologyFromOEMol(mol)
    # Create initial system
    system = ff.createSystem(topology, [mol], verbose=verbose)
    # Get initial energy before parameter modification
    positions = positions_from_oemol(mol)
    old_energy=get_energy(system, positions)

    # Get params for an angle
    params = ff.getParameter(smirks='[a,A:1]-[#6X4:2]-[a,A:3]')
    # Modify params
    params['k']='0.0'
    ff.setParameter(params, smirks='[a,A:1]-[#6X4:2]-[a,A:3]')
    # Write params
    ff.writeFile( tempfile.TemporaryFile(suffix='.ffxml') )
    # Make sure params changed energy! (Test whether they get rebuilt on system creation)
    system=ff.createSystem(topology, [mol], verbose=verbose)
    energy=get_energy(system, positions)
    if verbose:
        print("New energy/old energy:", energy, old_energy)
    if np.abs(energy-old_energy)<0.1:
        raise Exception("Error: Parameter modification did not change energy.")

def test_amber_roundtrip():
    """Save a System (a mixture) to AMBER, read back in, verify yields same energy and force terms."""

    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
    filename = get_data_filename(os.path.join('systems', 'packmol_boxes', 'cyclohexane_ethanol_0.4_0.6.pdb'))
    from simtk.openmm.app import PDBFile
    pdbfile = PDBFile(filename)
    mol2files = [get_data_filename(os.path.join('systems', 'monomers', 'ethanol.mol2')), get_data_filename(os.path.join('systems', 'monomers', 'cyclohexane.mol2'))]

    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    mols = []
    mol = oechem.OEMol()
    for mol2file in mol2files:
        ifs = oechem.oemolistream(mol2file)
        ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            oechem.OETriposAtomNames(mol)
            mols.append(oechem.OEGraphMol(mol))

    # setup system
    system = forcefield.createSystem( pdbfile.topology, mols)

    # Create ParmEd structure, save to AMBER
    a, prmtop = tempfile.mkstemp(suffix='.prmtop')
    a, crd = tempfile.mkstemp(suffix='.crd')
    save_system_to_amber( pdbfile.topology, system, pdbfile.positions, prmtop, crd)

    # Read back in and cross-check energies
    parm = parmed.load_file(prmtop, crd)
    ambersys = parm.createSystem(nonbondedMethod= app.NoCutoff, constraints = None, implicitSolvent = None)
    groups0, groups1, energy0, energy1 = compare_system_energies( pdbfile.topology, pdbfile.topology, ambersys, system, pdbfile.positions, verbose = False)

    # Remove temp files
    os.remove(prmtop)
    os.remove(crd)


def test_gromacs_roundtrip():
    """Save a System (a mixture) to GROMACS, read back in, verify yields same energy and force terms."""

    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
    filename = get_data_filename(os.path.join('systems', 'packmol_boxes', 'cyclohexane_ethanol_0.4_0.6.pdb'))
    from simtk.openmm.app import PDBFile
    pdbfile = PDBFile(filename)
    mol2files = [get_data_filename(os.path.join('systems', 'monomers', 'ethanol.mol2')), get_data_filename(os.path.join('systems', 'monomers', 'cyclohexane.mol2'))]

    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    mols = []
    mol = oechem.OEMol()
    for mol2file in mol2files:
        ifs = oechem.oemolistream(mol2file)
        ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            oechem.OETriposAtomNames(mol)
            mols.append(oechem.OEGraphMol(mol))

    # setup system
    system = forcefield.createSystem( pdbfile.topology, mols)

    # Create ParmEd structure, save to AMBER
    a, topfile = tempfile.mkstemp(suffix='.top')
    a, grofile = tempfile.mkstemp(suffix='.gro')
    save_system_to_gromacs( pdbfile.topology, system, pdbfile.positions, topfile, grofile)

    # Read back in and cross-check energies
    top = parmed.load_file(topfile)
    gro = parmed.load_file(grofile)
    gromacssys = top.createSystem(nonbondedMethod= app.NoCutoff, constraints = None, implicitSolvent = None)

    groups0, groups1, energy0, energy1 = compare_system_energies( pdbfile.topology, pdbfile.topology, gromacssys, system, pdbfile.positions, verbose = False)

    # Remove temp files
    os.remove(topfile)
    os.remove(grofile)


def test_tip3p_constraints():
    """Test that TIP3P distance costraints are correctly applied."""
    # TIP3P constrained distances.
    tip3p_oh_distance = 0.9572  # angstrom
    tip3p_hoh_angle = 104.52  # angle
    tip3p_hh_distance = tip3p_oh_distance * np.sin(np.radians(tip3p_hoh_angle/2)) * 2
    expected_distances = [tip3p_oh_distance, tip3p_oh_distance, tip3p_hh_distance]

    # Load tip3p molecule as OEMol.
    tip3p_mol2_filepath = get_data_filename(os.path.join('systems', 'monomers', 'tip3p_water.mol2'))
    tip3p_oemol = read_molecules(tip3p_mol2_filepath, verbose=False)[0]

    # Extract topology and positions.
    tip3p_topology = generateTopologyFromOEMol(tip3p_oemol)
    tip3p_positions = positions_from_oemol(tip3p_oemol)

    # Create tip3p water system.
    ff = ForceField('forcefield/tip3p.ffxml')
    system = ff.createSystem(tip3p_topology, [tip3p_oemol])

    # Run dynamics.
    integrator = openmm.VerletIntegrator(2.0*unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(tip3p_positions)
    integrator.step(50)

    # Constrained distances are correct.
    state = context.getState(getPositions=True)
    new_positions = state.getPositions(asNumpy=True) / unit.angstroms
    distances = []
    for atom_1, atom_2 in [(0, 1), (0, 2), (1, 2)]:  # pair of atoms O-H1, O-H2, H1-H2
        distances.append(np.linalg.norm(new_positions[atom_1] - new_positions[atom_2]))
    err_msg = 'expected distances [O-H1, O-H2, H1-H2]: {} A, new distances: {} A'
    assert np.allclose(expected_distances, distances), err_msg.format(expected_distances, distances)


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
        k_tip3p = 462750.4*unit.kilojoule_per_mole/unit.nanometers**2
        length_tip3p = 0.09572*unit.nanometers
        for force in system.getForces():
            if isinstance(force, openmm.HarmonicBondForce):
                force.addBond(particle1=0, particle2=1, length=length_tip3p, k=k_tip3p)
                force.addBond(particle1=0, particle2=2, length=length_tip3p, k=k_tip3p)


    monomers_dir = get_data_filename(os.path.join('systems', 'monomers'))

    # Load tip3p molecule as OEMol.
    tip3p_mol2_filepath = os.path.join(monomers_dir, 'tip3p_water.mol2')
    tip3p_oemol = read_molecules(tip3p_mol2_filepath, verbose=False)[0]

    # Create OpenMM-parametrized water ParmEd Structure.
    tip3p_openmm_ff = openmm.app.ForceField('tip3p.xml')
    tip3p_topology = generateTopologyFromOEMol(tip3p_oemol)
    tip3p_positions = positions_from_oemol(tip3p_oemol)
    tip3p_openmm_system = tip3p_openmm_ff.createSystem(tip3p_topology)

    # TODO remove this line when ParmEd#868 is fixed
    # Hack to make tip3p waters work with ParmEd Structure.createSystem.
    add_water_bonds(tip3p_openmm_system)

    tip3p_openmm_structure = parmed.openmm.topsystem.load_topology(
        tip3p_topology, tip3p_openmm_system, tip3p_positions)

    # Test case: monomer_filename
    test_cases = [('smirnoff99Frosst.ffxml', 'ethanol.mol2'),
                  ('smirnoff99Frosst.ffxml', 'methane.mol2')]

    for ff_name, molecule_filename in test_cases:
        # Load molecule as OEMol.
        molecule_filepath = os.path.join(monomers_dir, molecule_filename)
        molecule_oemol = read_molecules(molecule_filepath, verbose=False)[0]

        # Create molecule parametrized ParmEd args.
        molecule_ff = ForceField('forcefield/' + ff_name)
        molecule_args = forcefield_utils.create_system_from_molecule(molecule_ff, molecule_oemol)
        molecule_structure = parmed.openmm.topsystem.load_topology(*molecule_args)
        _, _, molecule_positions = molecule_args

        # Merge molecule and OpenMM TIP3P water.
        structure = tip3p_openmm_structure + molecule_structure
        structure.positions = np.append(tip3p_positions, molecule_positions, axis=0)
        system = structure.createSystem()
        energy_openmm = get_energy(system, structure.positions)

        # Test the creation of system with multiple force field files.
        ff = ForceField('forcefield/' + ff_name, 'forcefield/tip3p.ffxml')
        system = ff.createSystem(structure.topology, [tip3p_oemol, molecule_oemol])

        # TODO remove this line when ParmEd#868 is fixed
        # Add bonds to match ParmEd hack above.
        add_water_bonds(system)

        energy_smirnoff = get_energy(system, structure.positions)

        # The energies of the systems must be the same.
        assert np.isclose(energy_openmm, energy_smirnoff), 'OpenMM: {}, SMIRNOFF: {}'.format(
            energy_openmm, energy_smirnoff)


def test_component_combination():
    """Test that a system still yields the same energy after building it again out of its components."""

    # We've had issues where subsequent instances of a molecule might have zero charges
    # Here we'll try to catch this (and also explicitly check the charges) by re-building
    # a system out of its components

    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
    filename = get_data_filename(os.path.join('systems', 'packmol_boxes', 'cyclohexane_ethanol_0.4_0.6.pdb'))
    from simtk.openmm.app import PDBFile
    pdbfile = PDBFile(filename)
    mol2files = [get_data_filename(os.path.join('systems', 'monomers', 'ethanol.mol2')), get_data_filename(os.path.join('systems', 'monomers', 'cyclohexane.mol2'))]

    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    mols = []
    mol = oechem.OEMol()
    for mol2file in mol2files:
        ifs = oechem.oemolistream(mol2file)
        ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            oechem.OETriposAtomNames(mol)
            mols.append(oechem.OEGraphMol(mol))

    # setup system
    system = forcefield.createSystem( pdbfile.topology, mols, chargeMethod = 'OECharges_AM1BCCSym')

    # Make parmed structure
    structure = parmed.openmm.topsystem.load_topology( pdbfile.topology, system, pdbfile.positions)

    # Split the system, then re-compose it out of its components
    tmp = structure.split()
    strs, nums = [], []
    for s, n in tmp:
        strs.append(s)
        nums.append(n)
    nums = [ len(n) for n in nums]

    # Re-compose system from components
    new_structure = strs[0]*nums[0]
    for idx in range(1,len(nums)):
        new_structure += strs[idx]*nums[idx]
    # Swap in coordinates again
    new_structure.positions = structure.positions

    # Create System
    newsys = new_structure.createSystem(nonbondedMethod= app.NoCutoff, constraints = None, implicitSolvent = None)

    # Cross check energies
    groups0, groups1, energy0, energy1 = compare_system_energies( pdbfile.topology, pdbfile.topology, system, newsys, pdbfile.positions, verbose = False)

    # Also check that that the number of components is equal to the number I expect
    if not len(nums)==2:
        print("Error: Test system has incorrect number of components.")
        raise Exception('Incorrect number of components in cyclohexane/ethanol test system.')

    # Also check that none of residues have zero charge
    for resnr in range(len(structure.residues)):
        abscharges = [ abs(structure.residues[resnr].atoms[idx].charge) for idx in range(len(structure.residues[resnr].atoms))]
        if sum(abscharges)==0:
            raise Exception('Error: Residue %s in cyclohexane-ethanol test system has a charge of zero, which is incorrect.' % resnr)

def test_merge_system():
    """Test merging of a system created from AMBER and another created from SMIRNOFF."""

    #Create System from AMBER
    prefix = os.path.join('systems', 'amber', 'cyclohexane_ethanol_0.4_0.6')
    prmtop = get_data_filename( prefix+'.prmtop')
    incrd = get_data_filename( prefix+'.inpcrd')

    from openforcefield.typing.engines.smirnoff import create_system_from_amber
    topology0, system0, positions0 = create_system_from_amber( prmtop, incrd )

    from openeye import oechem
    # Load simple OEMol
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEthOH_c100.mol2'))
    mol = oechem.OEMol()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol )
    oechem.OETriposAtomNames(mol)

    # Load forcefield file
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
    topology1, system1, positions1 = create_system_from_molecule(forcefield, mol)

    merge_system( topology0, topology1, system0, system1, positions0, positions1, verbose=True )



if __name__ == '__main__':
    #test_smirks()
    test_merge_system()
