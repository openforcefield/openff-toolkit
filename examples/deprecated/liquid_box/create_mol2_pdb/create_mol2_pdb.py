#!/usr/bin/env python

from forcebalance.molecule import Molecule
from forcebalance.nifty import which
from openeye import oechem
import openmoltools 
import os, sys, time, argparse, subprocess

def CalculateMolecularWeight(mol):
    """
    Calculate the molecular weight for an OpenEye molecule.
    
    Parameters
    ----------
    mol : OEGraphMol
    
    Returns
    -------
    float
        Molecular weight in g/mol
    """
    result = 0.0
    for atom in mol.GetAtoms():
        elem = atom.GetAtomicNum()
        mass = atom.GetIsotope()
        if (elem != 0 and mass != 0):
            result += oechem.OEGetIsotopicWeight(elem,mass)
        else:
            result += oechem.OEGetAverageWeight(elem)
    return result

def CalculateBoxSize(nmol, molwt, density):
    """
    Calculate the size of a solvent box.
    
    Parameters
    ----------
    nmol : int
        Number of molecules desired for the box
    molwt : float
        Molecular weight in g/mol
    density : float
        Estimated density in kg/m3 (this should be about 40-50% lower than the real liquid density)
    
    Returns
    -------
    float
        Length of a cubic solvent box in nm.
    """
    # Calculate total mass of the box in kg
    mass = nmol * molwt / 1000 / 6.022e23
    volume = mass / density
    length = volume**(1./3)/1e-9
    return length

def GenerateBox(pdbin, pdbout, box, nmol, tries):
    """
    Call genbox. (Confirmed working with Gromacs version 4.6.7 and 5.1.4).
    Mainly checks whether genbox ran correctly.
    
    Parameters
    ----------
    pdbin : str
        Name of input PDB file containing a single molecule.
    pdbout : str
        Name of output PDB file containing solvent box.
    box : float
        Solvent box size, should be determined previously.
    nmol : int
        Number of molecules to go into the solvent box
    tries : int
        Parameter for genbox to try inserting each molecule (tries) times

    Returns
    -------
    None
        If successful, produces "pdbout" containing solvent box.
    """
    if which('gmx'):
        gmxcmd='gmx insert-molecules'
    elif which('genbox'):
        gmxcmd='genbox'
    else:
        raise RuntimeError('gmx and/or genbox not in PATH. Please source Gromacs environment variables.')

    fout=open('genbox.out', 'w')
    ferr=open('genbox.err', 'w')
    subprocess.Popen('%s -ci %s -o genbox.pdb -box %.3f %.3f %.3f -nmol %i -try %i'
                     % (gmxcmd, pdbin, box, box, box, nmol, tries),
                     shell=True, stdout=fout, stderr=ferr)
    fout.close()
    ferr.close()
    t0 = time.time()
    print("Running %s to create a solvent box..." % gmxcmd)
    print("Time elapsed: % .3f seconds" % (time.time() - t0))
    nmol_out = 0
    for line in open('genbox.err').readlines():
        if 'Output configuration contains' in line:
            nmol_out = int(line.split()[-2])
    if nmol_out == 0:
        raise RuntimeError('genbox failed to produce an output configuration')
    elif nmol_out != nmol:
        raise RuntimeError('genbox failed to create a box with %i molecules (actual %i); '
                           'please retry with increased box size or number of tries'  % (nmol, nmol_out))
    else:
        # genbox throws away the CONECT records in the PDB, this operation adds them back.
        M1 = Molecule(pdbin, build_topology=False)
        M = Molecule('genbox.pdb', build_topology=False)
        solventbox_bonds = []
        # Loop over the number of molecules in the solvent box
        for i in range(nmol):
            for j in M1.bonds:
                # Add the bonds for molecule number "i" in the solvent box
                solventbox_bonds.append((j[0]+i*M1.na, j[1]+i*M1.na))
        M.bonds = solventbox_bonds
        M.write(pdbout)
        print("-=# Output #=- Created %s containing solvent box with %i molecules and length %.3f" % (pdbout, nmol, box))

def run_create_mol2_pdb(**kwargs):

    nmol = kwargs['nmol']
    input_txt = kwargs['input']
    resname = kwargs['resname']
    density = kwargs['density']
    tries = kwargs['tries']

    # Disable Gromacs backup file creation
    os.environ['GMX_MAXBACKUP']="-1"
    
    smiles_string = open(input_txt).readlines()[0].strip()
    print("The following SMILES string will be converted: %s" % smiles_string)
    
    # create a new molecule
    oemol = oechem.OEGraphMol()
    # convert the SMILES string into a molecule
    if oechem.OESmilesToMol(oemol, smiles_string):
        # do something interesting with mol
        pass
    else:
        print("SMILES string was invalid!")
    
    # Add explicit
    oechem.OEAddExplicitHydrogens(oemol)
    
    # Generate a single conformer
    oemol = openmoltools.openeye.generate_conformers(oemol, max_confs=1)
    
    # Modify residue names
    oechem.OEPerceiveResidues(oemol, oechem.OEPreserveResInfo_All)
    for atom in oemol.GetAtoms():
        thisRes = oechem.OEAtomGetResidue(atom)
        thisRes.SetName(resname)
    
    # Write output files
    ofs = oechem.oemolostream()
    output_fnms = ['%s.mol2' % resname, '%s.pdb' % resname]
    for output_fnm in output_fnms:
        if not ofs.open(output_fnm):
            oechem.OEThrow.Fatal("Unable to create %s" % output_fnm)
        oechem.OEWriteMolecule(ofs, oemol)
        print("-=# Output #=- Created %s containing single molecule" % output_fnm)
    
    grams_per_mole = CalculateMolecularWeight(oemol)
    
    boxlen = CalculateBoxSize(nmol, grams_per_mole, density)
    GenerateBox('%s.pdb' % resname, '%s-box.pdb' % resname, boxlen, nmol, tries)

def main():
    """
    Provide a text file containing a single SMILES string and three-letter residue name. 
    Receive (res).pdb and (res).mol2 files containing a single molecule with conformation.
    Receive (res)-box.pdb containing a box with specified number

    Dependencies:
    OpenEye tools (for creating molecule from SMILES)
    openmoltools (for calling OpenEye to generate conformer)
    Gromacs 4.6.7 or 5.1.4 (for calling genbox to create solvent box)
    ForceBalance 1.5.x (for putting information back that was thrown away by genbox)
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--density', type=int, default=600, help='Specify target density of the solvent box; should be somewhat smaller than true liquid density due to imperfect packing.')
    parser.add_argument('--nmol', type=int, default=256, help='Specify desired number of molecules in the solvent box.')
    parser.add_argument('--tries', type=int, default=10, help='Pass number of tries per molecule to be passed to genbox. Higher = longer runtime but may achieve higher density.')
    parser.add_argument('input', type=str, help='Input file containing a single SMILES string')
    parser.add_argument('resname', type=str, help='Specify a custom residue name for the molecule.')
    print('%s called with the following command line:' % __file__)
    print(' '.join(sys.argv))
    args = parser.parse_args(sys.argv[1:])
    # Create the desired files (.mol2 file containing a single conformation and .pdb file containing solvent box).
    run_create_mol2_pdb(**vars(args))
    
if __name__ == "__main__":
    main()
    
