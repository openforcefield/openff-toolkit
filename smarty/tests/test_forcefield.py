from functools import partial
from smarty import ForceField
import smarty
import openeye
import os
from smarty.utils import get_data_filename
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

def test_read_ffxml():
    """Test reading of ffxml files.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

def test_create_system_molecules(verbose=False):
    """Test creation of a System object from small molecules.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

    from openeye import oechem
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEtOH-tripos.mol2.gz'))
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        from smarty.forcefield import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol(mol)
        system = forcefield.createSystem(topology, [mol], verbose=verbose)

def test_create_system_boxes(verbose=False):
    """Test creation of a System object from some boxes of mixed solvents.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

    # Read monomers
    mols = list()
    monomers = ['cyclohexane', 'ethanol', 'propane', 'methane', 'butanol']
    from openeye import oechem
    mol = oechem.OEGraphMol()
    for monomer in monomers:
        filename = get_data_filename(os.path.join('systems', 'monomers', monomer + '.mol2'))
        ifs = oechem.oemolistream(filename)
        while oechem.OEReadMolecule(ifs, mol):
            mols.append( oechem.OEGraphMol(mol) )
    print('%d reference molecules loaded' % len(mols))

    # Read systems.
    boxes = ['cyclohexane_ethanol_0.4_0.6.pdb', 'propane_methane_butanol_0.2_0.3_0.5.pdb']
    from smarty.forcefield import CutoffPeriodic
    from simtk.openmm.app import PDBFile
    for box in boxes:
        filename = get_data_filename(os.path.join('systems', 'packmol_boxes', box))
        pdbfile = PDBFile(filename)
        system = forcefield.createSystem(pdbfile.topology, mols, verbose=verbose, nonbondedMethod=CutoffPeriodic)

if __name__ == '__main__':
    test_create_system_boxes(verbose=True)
