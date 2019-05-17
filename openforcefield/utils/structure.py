#!/usr/bin/env python

"""
Utility subroutines for manipulating parmed Structure objects

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import time

import numpy as np

from simtk import openmm, unit

#=============================================================================================
# PARMED UTILITIES
#=============================================================================================

# TODO: Can we get rid of this, or convert it to use Molecule instead?
def generateSMIRNOFFStructure(oemol):
    """
    Given an OpenEye molecule (oechem.OEMol), create an OpenMM System and use to
    generate a ParmEd structure using the SMIRNOFF forcefield parameters.

    Parameters
    ----------
    oemol : openeye.oechem.OEMol
        OpenEye molecule

    Returns
    -------
    molecule_structure : parmed.Structure
        The resulting Structure

    """
    from openforcefield.topology import Molecule, Topology
    from openforcefield.typing.engines.smirnoff import ForceField

    off_mol = Molecule.from_openeye(oemol)
    off_top = Topology.from_molecules([off_mol])
    mol_ff = ForceField('test_forcefields/smirnoff99Frosst.offxml')

    # Create OpenMM System and Topology.
    omm_top = generateTopologyFromOEMol(oemol)

    # If it's a nonperiodic box, then we can't use default (PME) settings
    if omm_top.getPeriodicBoxVectors() is None:
        mol_ff.get_parameter_handler("Electrostatics", {})._method = 'Coulomb'

    system = mol_ff.create_openmm_system(off_top)

    # Convert to ParmEd structure.
    import parmed
    xyz = extractPositionsFromOEMol(oemol)
    molecule_structure = parmed.openmm.load_topology(omm_top, system, xyz=xyz)

    return molecule_structure

def generateProteinStructure(proteinpdb, protein_forcefield='amber99sbildn.xml', solvent_forcefield='tip3p.xml'):
    """
    Given an OpenMM PDBFile, create the OpenMM System of the protein using OpenMM's ForceField and then generate the parametrized ParmEd Structure of the protein.

    Parameters
    ----------
    proteinpdb : simtk.openmm.app.PDBFile object,
        Loaded PDBFile object of the protein.
    protein_forcefield : xml file, default='amber99sbildn.xml'
        Forcefield parameters for protein
    solvent_forcefield : xml file, default='tip3p.xml'
        Forcefield parameters for solvent

    Returns
    -------
    solv_structure : parmed.structure.Structure
        The parameterized Structure of the protein with solvent molecules. (No ligand).

    """
    # Generate protein Structure object using OpenMM ForceField
    from simtk.openmm import app
    import parmed
    forcefield = openmm.app.ForceField(protein_forcefield, solvent_forcefield)
    protein_system = forcefield.createSystem( proteinpdb.topology )
    protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                    protein_system,
                                                    xyz=proteinpdb.positions)
    return protein_structure

def combinePositions(proteinPositions, molPositions):
    """
    Loops through the positions from the ParmEd structures of the protein and ligand,
    divides by unit.angstroms which will ensure both positions arrays are in the same units.

    Parameters
    ----------
    proteinPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the protein Structure.
    molPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the molecule Structure.

    Returns
    -------
    positions : list of 3-element Quantity tuples.
        ex. unit.Quantity(positions, positions_unit)
        Combined positions of the protein and molecule Structures.
    """

    positions_unit = unit.angstroms
    positions0_dimensionless = np.array(proteinPositions / positions_unit)
    positions1_dimensionless = np.array(molPositions / positions_unit)
    coordinates = np.vstack(
        (positions0_dimensionless, positions1_dimensionless))
    natoms = len(coordinates)
    positions = np.zeros([natoms, 3], np.float32)
    for index in range(natoms):
            (x, y, z) = coordinates[index]
            positions[index, 0] = x
            positions[index, 1] = y
            positions[index, 2] = z
    positions = unit.Quantity(positions, positions_unit)
    return positions

def mergeStructure(proteinStructure, molStructure):
    """
    Combines the parametrized ParmEd structures of the protein and ligand to
    create the Structure for the protein:ligand complex, while retaining the SMIRNOFF
    parameters on the ligand. Preserves positions and box vectors.
    (Not as easily achieved using native OpenMM tools).

    Parameters
    ----------
    proteinStructure : parmed.structure.Structure
        The parametrized structure of the protein.
    moleculeStructure : parmed.structure.Structure
        The parametrized structure of the ligand.

    Returns
    -------
    structure : parmed.structure.Structure
        The parametrized structure of the protein:ligand complex.
    """
    structure = proteinStructure + molStructure
    positions = combinePostions(proteinStructure.positions, molStructure.positions)
    # Concatenate positions arrays (ensures same units)
    structure.positions = positions
    # Restore original box vectors
    structure.box = proteinStructure.box
    return structure

# TODO: Migrate into Molecule
def generateTopologyFromOEMol(molecule):
    """
    Generate an OpenMM Topology object from an OEMol molecule.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule from which a Topology object is to be generated.

    Returns
    -------
    topology : simtk.openmm.app.Topology
        The Topology object generated from `molecule`.

    """
    from openeye import oechem

    # Avoid manipulating the molecule
    mol = oechem.OEMol(molecule)

    # Create a Topology object with one Chain and one Residue.
    from simtk.openmm.app import Topology
    topology = Topology()
    chain = topology.addChain()
    resname = mol.GetTitle()
    residue = topology.addResidue(resname, chain)

    # Make sure the atoms have names, otherwise bonds won't be created properly below
    if any([atom.GetName() =='' for atom in mol.GetAtoms()]):
        oechem.OETriposAtomNames(mol)
    # Check names are unique; non-unique names will also cause a problem
    atomnames = [ atom.GetName() for atom in mol.GetAtoms() ]
    if any( atomnames.count(atom.GetName())>1 for atom in mol.GetAtoms()):
        raise Exception("Error: Reference molecule must have unique atom names in order to create a Topology.")

    # Create atoms in the residue.
    for atom in mol.GetAtoms():
        name = atom.GetName()
        element = openmm.app.element.Element.getByAtomicNumber(atom.GetAtomicNum())
        topology.addAtom(name, element, residue)

    # Create bonds.
    atoms = { atom.name : atom for atom in topology.atoms() }
    for bond in mol.GetBonds():
        aromatic = None
        if bond.IsAromatic(): aromatic = 'Aromatic'
        # Add bond, preserving order assessed by OEChem
        topology.addBond(atoms[bond.GetBgn().GetName()], atoms[bond.GetEnd().GetName()], type=aromatic, order=bond.GetOrder())

    return topology


# TODO: Do we need this?
def normalize_molecules(molecules):
    """
    Normalize all molecules in specified set.

    ARGUMENTS

    molecules (list of OEMol) - molecules to be normalized (in place)

    """
    from openeye import oechem

    # Add explicit hydrogens.
    for molecule in molecules:
        oechem.OEAddExplicitHydrogens(molecule)

    # Build a conformation for all molecules with Omega.
    print("Building conformations for all molecules...")
    from openeye import oeomega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetFromCT(True)
    for molecule in molecules:
        #omega.SetFixMol(molecule)
        omega(molecule)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("%.3f s elapsed" % elapsed_time)

    # Regularize all molecules through writing as mol2.
    print("Regularizing all molecules...")
    ligand_mol2_dirname  = os.path.dirname(mcmcDbName) + '/mol2'
    if( not os.path.exists( ligand_mol2_dirname ) ):
        os.makedirs(ligand_mol2_dirname)
    ligand_mol2_filename = ligand_mol2_dirname + '/temp' + os.path.basename(mcmcDbName) + '.mol2'
    start_time = time.time()
    omolstream = oechem.oemolostream(ligand_mol2_filename)
    for molecule in molecules:
        # Write molecule as mol2, changing molecule through normalization.
        oechem.OEWriteMolecule(omolstream, molecule)
    omolstream.close()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("%.3f s elapsed" % elapsed_time)

    # Assign AM1-BCC charges.
    print("Assigning AM1-BCC charges...")
    start_time = time.time()
    for molecule in molecules:
        # Assign AM1-BCC charges.
        if molecule.NumAtoms() == 1:
            # Use formal charges for ions.
            OEFormalPartialCharges(molecule)
        else:
            # Assign AM1-BCC charges for multiatom molecules.
            OEAssignPartialCharges(molecule, OECharges_AM1BCC, False) # use explicit hydrogens
        # Check to make sure we ended up with partial charges.
        if OEHasPartialCharges(molecule) == False:
            print("No charges on molecule: '%s'" % molecule.GetTitle())
            print("IUPAC name: %s" % OECreateIUPACName(molecule))
            # TODO: Write molecule out
            # Delete themolecule.
            molecules.remove(molecule)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("%.3f s elapsed" % elapsed_time)
    print("%d molecules remaining" % len(molecules))

    return

# TODO: Migrate into Molecule
def read_molecules(file_path, verbose=True):
    """
    Read molecules from an OpenEye-supported file.

    Parameters
    ----------
    file_path : str
        Filename from which molecules are to be read (e.g. mol2, sdf)

    Returns
    -------
    molecules : list of OEMol
        List of molecules read from file

    """
    from openeye import oechem
    from openforcefield.utils import get_data_file_path

    if not os.path.exists(file_path):
        built_in = get_data_file_path('molecules/%s' % file_path)
        if not os.path.exists(built_in):
            raise Exception("File '%s' not found." % file_path)
        file_path = built_in

    if verbose: print("Loading molecules from '%s'..." % file_path)
    start_time = time.time()
    molecules = list()
    input_molstream = oechem.oemolistream(file_path)

    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    input_molstream.SetFlavor(oechem.OEFormat_MOL2, flavor)

    molecule = oechem.OECreateOEGraphMol()
    while oechem.OEReadMolecule(input_molstream, molecule):
        # If molecule has no title, try getting SD 'name' tag
        if molecule.GetTitle() == '':
            name = oechem.OEGetSDData(molecule, 'name').strip()
            molecule.SetTitle(name)
        # Append to list.
        molecule_copy = oechem.OEMol(molecule)
        molecules.append(molecule_copy)
    input_molstream.close()
    if verbose: print("%d molecules read" % len(molecules))
    end_time = time.time()
    elapsed_time = end_time - start_time
    if verbose: print("%.3f s elapsed" % elapsed_time)

    return molecules

# TODO: Migrate into Molecule
def setPositionsInOEMol(oemol, positions):
    """Set the positions in an OEMol using a position array with units from simtk.unit, i.e. from OpenMM. Atoms must have same order.

    Arguments:
    ---------
    oemol : OEMol
        OpenEye molecule
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """
    from openeye import oechem

    if oemol.NumAtoms() != len(positions): raise ValueError("Number of atoms in molecule does not match length of position array.")
    pos_unitless = positions/unit.angstroms

    coordlist = []
    for idx in range(len(pos_unitless)):
        for j in range(3):
            coordlist.append(pos_unitless[idx][j])
    oemol.SetCoords(oechem.OEFloatArray(coordlist))

# TODO: Migrate into Molecule
def extractPositionsFromOEMol(oemol):
    """Get the positions from an OEMol and return in a position array with units via simtk.unit, i.e. foramtted for OpenMM.
    Adapted from choderalab/openmoltools test function extractPositionsFromOEMOL

    Arguments:
    ----------
    oemol : OEMol
        OpenEye molecule

    Returns:
    --------
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """
    positions = unit.Quantity(np.zeros([oemol.NumAtoms(), 3], np.float32), unit.angstroms)
    coords = oemol.GetCoords()
    for index in range(oemol.NumAtoms()):
        positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
    return positions

def read_typelist(file_path):
    """
    Read a parameter type or decorator list from a file.
    Lines in these files have the format
    "SMARTS/SMIRKS  shorthand"
    lines beginning with '%' are ignored

    Parameters
    ----------
    file_path : str
        Path and name of file to be read
        Could be file in openforcefield/data/

    Returns
    -------
    typelist : list of tuples
        Typelist[i] is element i of the typelist in format (smarts, shorthand)
    """
    from openforcefield.utils import get_data_file_path

    if file_path is None:
        return None

    if not os.path.exists(file_path):
        built_in = get_data_file_path(file_path)
        if not os.path.exists(built_in):
            raise Exception("File '%s' not found." % file_path)
        file_path = built_in

    typelist = list()
    ifs = open(file_path)
    lines = ifs.readlines()
    used_typenames = list()

    for line in lines:
        # Strip trailing comments
        index = line.find('%')
        if index != -1:
            line = line[0:index]

        # Split into tokens.
        tokens = line.split()
        # Process if we have enough tokens
        if len(tokens) >= 2:
            smarts = tokens[0]
            typename = ' '.join(tokens[1:])
            if typename not in used_typenames:
                typelist.append([smarts,typename])
                used_typenames.append(typename)
            else:
                raise Exception("Error in file '%s' -- each entry must "
                        "have a unique name." % file_path )

    ifs.close()

    return typelist

# TODO: This only works with the OpenEye toolkit installed; replace with Molecule API
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
    from openeye import oechem, oeomega
    if mol.GetDimension() != 3:
        # Assign coordinates
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(1)
        omega.SetIncludeInput(False)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
        omega(mol)  # generate conformation

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

#=============================================================================
# OPENMM MERGING AND EXPORTING UTILITY FUNCTIONS
#=============================================================================

# TODO: Reorganize this file, moving exporters to openforcefield.exporters

def get_molecule_parameterIDs(molecules, forcefield):
    """Process a list of molecules with a specified SMIRNOFF ffxml file and determine which parameters are used by
    which molecules, returning collated results.

    Parameters
    ----------
    molecules : list of openforcefield.topology.Molecule
        List of molecules (with explicit hydrogens) to parse
    forcefield : openforcefield.typing.engines.smirnoff.ForceField
        The ForceField to apply

    Returns
    -------
    parameters_by_molecule : dict
        Parameter IDs used in each molecule, keyed by isomeric SMILES
        generated from provided OEMols. Each entry in the dict is a list
        which does not necessarily have unique entries; i.e. parameter IDs
        which are used more than once will occur multiple times.

    parameters_by_ID : dict
        Molecules in which each parameter ID occur, keyed by parameter ID.
        Each entry in the dict is a set of isomeric SMILES for molecules
        in which that parameter occurs. No frequency information is stored.

    """
    from openforcefield.topology import Topology
    # Create storage
    parameters_by_molecule = dict()
    parameters_by_ID = dict()

    # Generate isomeric SMILES for each molecule, ensuring all molecules are unique
    isosmiles = [ molecule.to_smiles() for molecule in molecules ]
    already_seen = set()
    duplicates = set(smiles for smiles in isosmiles if smiles in already_seen or already_seen.add(smiles))
    if len(duplicates) > 0:
        raise ValueError("Error: get_molecule_parameterIDs has been provided a list of oemols which contains some duplicates: {}".format(duplicates))

    # Assemble molecules into a Topology
    topology = Topology()
    for molecule in molecules:
        topology.add_molecule(molecule)

    # Label molecules
    labels = forcefield.label_molecules(topology)

    # Organize labels into output dictionary by looping over all molecules/smiles
    for idx in range(len(isosmiles)):
        # Pull smiles, initialize storage
        smi = isosmiles[idx]
        parameters_by_molecule[smi] = []

        # Organize data for this molecule
        data = labels[idx]
        for force_type in data.keys():
            for atom_indices, parameter_type in data[force_type].items():

                pid = parameter_type.id
                # Store pid to molecule
                parameters_by_molecule[smi].append(pid)

                # Store which molecule this pid occurred in
                if pid not in parameters_by_ID:
                    parameters_by_ID[pid] = set()
                    parameters_by_ID[pid].add(smi)
                else:
                    parameters_by_ID[pid].add(smi)

    return parameters_by_molecule, parameters_by_ID

def getMolParamIDToAtomIndex(molecule, forcefield):
    """Take a Molecule and a SMIRNOFF forcefield object and return a dictionary, keyed by parameter ID, where each entry is a tuple of ( smirks, [[atom1, ... atomN], [atom1, ... atomN]) giving the SMIRKS corresponding to that parameter ID and a list of the atom groups in that molecule that parameter is applied to.

    Parameters
    ----------
    molecule : openforcefield.topology.Molecule
        Molecule to investigate
    forcefield : ForceField
        SMIRNOFF ForceField object (obtained from an ffxml via ForceField(ffxml)) containing FF of interest.

    Returns
    -------
    param_usage : dictionary
        Dictionary, keyed by parameter ID, where each entry is a tuple of ( smirks, [[atom1, ... atomN], [atom1, ... atomN]) giving the SMIRKS corresponding to that parameter ID and a list of the atom groups in that molecule that parameter is applied to.

    """

    topology = Topology()
    topology.add_molecule(molecule)
    labels = ff.labal_molecules(topology)

    param_usage = {}
    for mol_entry in range(len(labels)):
        for force in labels[mol_entry].keys():
            for (atom_indices, pid, smirks) in labels[mol_entry][force]:
                if not pid in param_usage:
                    param_usage[pid] = (smirks, [atom_indices])
                else:
                    param_usage[pid][1].append( atom_indices )

    return param_usage

def merge_system(topology0, topology1, system0, system1, positions0, positions1, label0="AMBER system", label1="SMIRNOFF system", verbose=True):
    """Merge two given OpenMM systems. Returns the merged OpenMM System.

    Parameters
    ----------
    topology0 : OpenMM Topology
        Topology of first system (i.e. a protein)
    topology1 : OpenMM Topology
        Topology of second system (i.e. a ligand)
    system0 : OpenMM System
        First system for merging (usually from AMBER)
    system1 : OpenMM System
        Second system for merging (usually from SMIRNOFF)
    positions0 : simtk.unit.Quantity wrapped
        Positions to use for energy evaluation comparison
    positions1 (optional) : simtk.unit.Quantity wrapped (optional)
        Positions to use for second OpenMM system
    label0 (optional) : str
        String labeling system0 for output. Default, "AMBER system"
    label1 (optional) : str
        String labeling system1 for output. Default, "SMIRNOFF system"
    verbose (optional) : bool
        Print out info on topologies, True/False (default True)

    Returns
    ----------
    topology : OpenMM Topology
    system : OpenMM System
    positions: unit.Quantity position array
    """

    #Load OpenMM Systems to ParmEd Structures
    import parmed
    structure0 = parmed.openmm.load_topology( topology0, system0 )
    structure1 = parmed.openmm.load_topology( topology1, system1 )

    #Merge parameterized Structure
    structure = structure0 + structure1
    topology = structure.topology

    #Concatenate positions arrays
    positions_unit = unit.angstroms
    positions0_dimensionless = np.array( positions0 / positions_unit )
    positions1_dimensionless = np.array( positions1 / positions_unit )

    coordinates = np.vstack((positions0_dimensionless,positions1_dimensionless))
    natoms = len(coordinates)
    positions = np.zeros([natoms,3], np.float32)
    for index in range(natoms):
        (x,y,z) = coordinates[index]
        positions[index,0] = x
        positions[index,1] = y
        positions[index,2] = z
    positions = unit.Quantity(positions, positions_unit)

    #Generate merged OpenMM system
    system = structure.createSystem()

    if verbose:
        print("Generating ParmEd Structures...\n \t{}: {}\n \t{}: {}\n".format(label0, structure0, label1, structure1))
        print("Merged ParmEd Structure: {}".format( structure ))

    return topology, system, positions

def save_system_to_amber(openmm_topology, system, positions, prmtop, inpcrd):
    """Save an OpenMM System, with provided topology and positions, to AMBER prmtop and coordinate files.

    Parameters
    ----------
    openmm_topology : OpenMM Topology
        Topology of the system to be saved, perhaps as loaded from a PDB file or similar.
    system : OpenMM System
        Parameterized System to be saved, containing components represented by Topology
    positions : unit.Quantity position array
        Position array containing positions of atoms in topology/system
    prmtop : filename
        AMBER parameter file name to write
    inpcrd : filename
        AMBER coordinate file name (ASCII crd format) to write

    """

    import parmed
    structure = parmed.openmm.topsystem.load_topology(openmm_topology, system, positions)
    structure.save(prmtop, overwrite = True, format="amber")
    structure.save(inpcrd, format='rst7', overwrite = True)

def save_system_to_gromacs(openmm_topology, system, positions, top, gro):
    """Save an OpenMM System, with provided topology and positions, to AMBER prmtop and coordinate files.

    Parameters
    ----------
    openmm_topology : OpenMM Topology
        Topology of the system to be saved, perhaps as loaded from a PDB file or similar.
    system : OpenMM System
        Parameterized System to be saved, containing components represented by Topology
    positions : unit.Quantity position array
        Position array containing positions of atoms in topology/system
    top : filename
        GROMACS topology file name to write
    gro : filename
        GROMACS coordinate file name (.gro format) to write

    """

    import parmed
    structure = parmed.openmm.topsystem.load_topology(openmm_topology, system, positions )
    structure.save(top, overwrite = True, format="gromacs")
    structure.save(gro, overwrite = True, format="gro")
