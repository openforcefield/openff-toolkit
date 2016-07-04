from functools import partial
from smarty import ForceField
import smarty
import openeye
import os
from smarty.utils import get_data_filename
from simtk.openmm import app
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology
from simtk import unit, openmm
import numpy as np
from StringIO import StringIO

# This is a test forcefield that is not meant for actual use.
# It just tests various capabilities.
ffxml_contents = """\
<?xml version="1.0"?>

<SMARFF>

<!-- Header block (optional) -->
<Date>Date: Tue May  3 2016</Date>
<Author>C. I. Bayly, OpenEye Scientific Software</Author>
<!-- This is a comment -->

<HarmonicBondForce length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
   <Bond smirks="[#6X4:1]-[#6X4:2]" length="1.526" k="620.0"/> <!-- CT-CT from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#6X4:1]-[#1:2]" length="1.090" k="680.0"/> <!-- CT-H_ from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#8:1]~[#1:2]" length="1.410" k="640.0"/> <!-- DEBUG O-H -->
   <Bond smirks="[#6X4:1]-[O&amp;X2&amp;H1:2]" length="1.410" k="640.0"/> <!-- CT-OH from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#6X4:1]-[O&amp;X2&amp;H0:2]" length="1.370" k="640.0"/> <!--CT-OS from frcmod.Frosst_AlkEthOH -->
   <Bond smirks="[#8X2:1]-[#1:2]" length="0.960" k="1106.0"/> <!-- OH-HO from frcmod.Frosst_AlkEthOH -->
</HarmonicBondForce>

<HarmonicAngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
   <Angle smirks="[a,A:1]-[#6X4:2]-[a,A:3]" angle="109.50" k="100.0"/> <!-- consensus matches all X-Csp3-X -->
   <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.50" k="70.0"/> <!-- H1-CT-H1 from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]" angle="109.50" k="80.0"/> <!-- CT-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#8X2:1]-[#6X4:2]-[#8X2:3]" angle="109.50" k="140.0"/> <!-- O_-CT-O_ from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#6X4:1]-[#8X2:2]-[#1:3]" angle="108.50" k="110.0"/> <!-- CT-OH-HO from frcmod.Frosst_AlkEthOH -->
   <Angle smirks="[#6X4:1]-[#8X2:2]-[#6X4:3]" angle="109.50" k="120.0"/> <!-- CT-OS-CT from frcmod.Frosst_AlkEthOH -->
</HarmonicAngleForce>

<PeriodicTorsionForce phase_unit="degrees" k_unit="kilocalories_per_mole">
   <Proper smirks="[a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]" periodicity1="3" phase1="0.0" k1="1.40"/> <!-- X -CT-CT-X from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[a,A:1]-[#6X4:2]-[#8X2:3]-[#1:4]" idivf1="3" periodicity1="3" phase1="0.0" k1="0.50"/> <!--X -CT-OH-X from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="1.15"/> <!-- X -CT-OS-X from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]" periodicity1="3" phase1="0.0" k1="0.15"/> <!-- HC-CT-CT-HC from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#1:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" periodicity1="3" phase1="0.0" k1="0.16"/> <!-- HC-CT-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#8X2:3]-[#1:4]" idivf1="3" periodicity1="3" phase1="0.0" k1="0.16" idivf2="1" periodicity2="1" phase2="0.0" k2="0.25"/> <!-- HO-OH-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" idivf1="3" periodicity1="3" phase1="0.0" k1="0.18" idivf2="2" periodicity2="2" phase2="180.0" k2="0.25" periodicity3="1" phase3="180.0" k3="0.20"/> <!-- CT-CT-CT-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]" idivf1="3" periodicity1="3" phase1="0.0" k1="0.383" periodicity2="2" phase2="180.0" k2="0.1"/> <!-- CT-CT-OS-CT from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#6X4:1]-[#8X4:2]-[#6X4:3]-[O&amp;X2&amp;H0:4]" periodicity1="3" phase1="0.0" k1="0.10" periodicity2="2" phase2="180.0" k2="0.85" periodicity3="1" phase3="180.0" k3="1.35"/> <!-- CT-OS-CT-OS from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#8X2:1]-[#6X4:2]-[#6X4:3]-[#8X2:4]" periodicity1="3" phase1="0.0" k1="0.144" periodicity2="2" phase2="0.0" k2="0.175"/> <!-- O_-CT-CT-O_ from frcmod.Frosst_AlkEthOH -->
   <Proper smirks="[#8X2:1]-[#6X4:2]-[#6X4:3]-[#1:4]" periodicity1="3" phase1="0.0" k1="0.0" periodicity2="1" phase2="0.0" k2="0.25"/> <!-- O_-CT-CT-O_ from frcmod.Frosst_AlkEthOH; discrepancy with parm@frosst with H2,H3-CT-CT-O_ per C Bayly -->
   <Improper smirks="[a,A:1]~[#6X3:2]([a,A:3])~[OX1:4]" periodicity1="2" phase1="180.0" k1="10.5"/> <!-- X -X -C -O  from frcmod.Frosst_AlkEthOH; none in set but here as format placeholder -->
</PeriodicTorsionForce>

<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole">
   <!-- sigma is in angstroms, epsilon is in kcal/mol -->
   <Atom smirks="[#1:1]" sigma="1.4870" epsilon="0.0157"/> <!-- making HC the generic hydrogen -->
   <Atom smirks="[$([#1]-C):1]" rmin_half="1.4870" epsilon="0.0157"/> <!-- HC from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[$([#1]-C-[#7,#8,F,#16,Cl,Br]):1]" rmin_half="1.3870" epsilon="0.0157"/> <!-- H1 from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[$([#1]-C(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]):1]" rmin_half="1.2870" epsilon="0.0157"/> <!--H2 from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[$([#1]-C(-[#7,#8,F,#16,Cl,Br])(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]):1]" rmin_half="1.1870" epsilon="0.0157"/> <!--H3 from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[#1$(*-[#8]):1]" rmin_half="0.0000" epsilon="0.0000"/> <!-- HO from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[#6:1]" rmin_half="1.9080" epsilon="0.1094"/> <!-- making CT the generic carbon -->
   <Atom smirks="[#6X4:1]" rmin_half="1.9080" epsilon="0.1094"/> <!-- CT from frcmod.Frosst_AlkEthOH-->
   <Atom smirks="[#8:1]" rmin_half="1.6837" epsilon="0.1700"/> <!-- making OS the generic oxygen -->
   <Atom smirks="[#8X2:1]" rmin_half="1.6837" epsilon="0.1700"/> <!-- OS from frcmod.Frosst_AlkEthOH -->
   <Atom smirks="[#8X2+0$(*-[#1]):1]" rmin_half="1.7210" epsilon="0.2104"/> <!-- OH from frcmod.Frosst_AlkEthOH -->
</NonbondedForce>

</SMARFF>
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

def test_read_ffxml():
    """Test reading of ffxml files.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

def check_system_creation_from_molecule(forcefield, mol, verbose=False):
    """
    Generate a System from the given OEMol and SMIRFF forcefield and check that its energy is finite.

    Parameters
    ----------
    forcefield : smarty.ForceField
        SMIRFF forcefield
    mol : oechem.OEMol
        Molecule to test (need not have coordinates)

    """

    from smarty.forcefield import generateTopologyFromOEMol
    topology = generateTopologyFromOEMol(mol)
    system = forcefield.createSystem(topology, [mol], verbose=verbose)
    # Test energy computation.
    positions = positions_from_oemol(mol)
    check_energy_is_finite(system, positions)

def check_system_creation_from_topology(forcefield, topology, mols, positions, verbose=False):
    """
    Generate a System from the given topology, OEMols matching the contents of the topology, and SMIRFF forcefield and check that its energy is finite.

    Parameters
    ----------
    forcefield : smarty.ForceField
        SMIRFF forcefield
    topology : simtk.openmm.app.Topology
        Topology to construct system from
    mols : list of oechem.OEMol
        Reference molecules
    positions : simtk.unit.Quantity with dimension (natoms,3) with units of length

    """
    from smarty.forcefield import CutoffPeriodic
    system = forcefield.createSystem(topology, mols, verbose=verbose, nonbondedMethod=CutoffPeriodic)
    # Test energy computation.
    check_energy_is_finite(system, positions)

def check_AlkEtOH(forcefield, description="", verbose=False):
    """Test creation of System from AlkEtOH small molecules.
    """
    from openeye import oechem
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEtOH-tripos.mol2.gz'))
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        args = { 'verbose' : verbose }
        f = partial(check_system_creation_from_molecule, forcefield, mol, **args)
        f.description ='Testing creation of system object from small molecules (%s) %s' % (mol.GetTitle(), description)
        yield f

def test_create_system_molecules_features(verbose=False):
    """Test creation of a System object from small molecules to test various ffxml features
    """
    ffxml = StringIO(ffxml_contents)
    forcefield = ForceField(ffxml)
    check_AlkEtOH(forcefield, description="to test ffxml features", verbose=verbose)

def test_create_system_molecules_parmatfrosst(verbose=False):
    """Test creation of a System object from small molecules to test parm@frosst forcefield.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))
    check_AlkEtOH(forcefield, "to test Parm@Frosst parameters", verbose=verbose)

def check_boxes(forcefield, description="", verbose=False):
    """Test creation of System from boxes of mixed solvents.
    """
    # Read monomers
    mols = list()
    monomers = ['cyclohexane', 'ethanol', 'propane', 'methane', 'butanol']
    from openeye import oechem
    mol = oechem.OEGraphMol()
    for monomer in monomers:
        filename = get_data_filename(os.path.join('systems', 'monomers', monomer + '.sdf'))
        ifs = oechem.oemolistream(filename)
        while oechem.OEReadMolecule(ifs, mol):
            oechem.OETriposAtomNames(mol)
            mols.append( oechem.OEGraphMol(mol) )
    print('%d reference molecules loaded' % len(mols))

    # Read systems.
    boxes = ['cyclohexane_ethanol_0.4_0.6.pdb', 'propane_methane_butanol_0.2_0.3_0.5.pdb']
    from simtk.openmm.app import PDBFile
    for box in boxes:
        filename = get_data_filename(os.path.join('systems', 'packmol_boxes', box))
        pdbfile = PDBFile(filename)
        f = partial(check_system_creation_from_topology, forcefield, pdbfile.topology, mols, pdbfile.positions, verbose=verbose)
        f.description = 'Test creation of System object from %s' % box
        yield f

def test_create_system_boxes_features(verbose=False):
    """Test creation of a System object from some boxes of mixed solvents to test ffxml features.
    """
    ffxml = StringIO(ffxml_contents)
    forcefield = ForceField(ffxml)
    check_boxes(forcefield, description="to test Parm@frosst parameters", verbose=verbose)

def test_create_system_boxes_parmatfrosst(verbose=False):
    """Test creation of a System object from some boxes of mixed solvents to test parm@frosst forcefield.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))
    check_boxes(forcefield, description="to test Parm@frosst parameters", verbose=verbose)

if __name__ == '__main__':
    #test_smirks()
    test_create_system_boxes(verbose=True)
