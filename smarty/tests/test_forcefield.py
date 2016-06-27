from functools import partial
from smarty import ForceField
import smarty
import openeye
from smarty.utils import get_data_filename
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

def test_read_ffxml():
    """Test reading of ffxml files.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

def test_create_system():
    """Test creation of a System object from small molecules.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

    from openeye import oechem
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEtOH-tripos.mol2.gz'))
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        from smarty.forcefield import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol(mol)
        system = forcefield.createSystem(topology, [mol])
