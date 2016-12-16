from simtk.openmm import app
import simtk.openmm as mm
from simtk.openmm.app import Topology
import numpy as np
import openmoltools
from openeye import oechem
from simtk.openmm import unit
from openmoltools import forcefield_generators
from openmoltools.system_checker import *
from smarty import *
import parmed as pmd


def gen_proteinStructure():
    #Create OpeMM Topology from PDB
    pdbfilename = "receptor.pdbfixer.pdb"
    pdbfile = app.PDBFile(pdbfilename)

    #Initialize a forcefield
    gaff_xml_filename = openmoltools.utils.get_data_filename("parameters/gaff.xml")
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml', gaff_xml_filename)
    forcefield.registerTemplateGenerator(forcefield_generators.gaffTemplateGenerator)

    #Setup proteim
    modeller = app.Modeller(pdbfile.topology, pdbfile.positions)
    modeller.addHydrogens(forcefield, pH=7.0)

    ###Solvate and Ionize...?

    #Generate OpenMM system with FF parameters
    topology = modeller.getTopology()
    system = forcefield.createSystem( modeller.getTopology() )

    #Load System to ParmEd Structure
    structure = pmd.openmm.load_topology(topology, system)

    return structure

def gen_ligandStructure():
    #Create OpenMM Topology from OEMol for ligand
    molfilename = "ligand.tripos.mol2"

    # Load OEMol
    mol = OEGraphMol()
    ifs = oemolistream(molfilename)
    flavor = OEIFlavor_Generic_Default | OEIFlavor_MOL2_Default | OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( OEFormat_MOL2, flavor)
    OEReadMolecule(ifs, mol )
    OETriposAtomNames(mol)

    topology = generateTopologyFromOEMol(mol)
    positions = unit.Quantity( np.zeros([mol.NumAtoms(),3], np.float64), unit.angstroms)
    for (index, atom) in enumerate(mol.GetAtoms()):
        [x,y,z] = mol.GetCoords(atom)
        positions[index,0] = x*unit.angstroms
        positions[index,1] = y*unit.angstroms
        positions[index,2] = z*unit.angstroms

    #Initialize SMIRFF parameters
    # using smarty ForceField object
    smirnoff_xml_filename ="smirff99Frosst.ffxml"
    forcefield = ForceField(smirnoff_xml_filename)

    #Create OpenMM System
    system = forcefield.createSystem(topology, [mol])

    #Load System to ParmEd Structure
    structure = pmd.openmm.load_topology(topology, system)

    return structure


#Combine protein-ligand structures
pl_struct = gen_proteinStructure() + gen_ligandStructure()
pl_topology = pl_struct.topology
pl_system = pl_struct.createSystem()

pl_struct.save('pl-test.top', overwrite=True)
