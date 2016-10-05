from openeye import oechem
from openeye import oeomega
import smarty

def check_valence(mol):
    """
    Given an OEMol it returns True if no small (atomic number < 10)
    has a valence greater than 4
    """
    for atom in mol.GetAtoms():
        atomNum = atom.GetAtomicNum()
        valence = atom.GetValence()
        if atomNum <= 10:
            if valence > 4:
                print("Found a #%i atom with valence %i in molecule %s" % (atomNum, valence, oechem.OECreateIsoSmiString(mol)))
                return False
    return True

def generate_conformer(mol):
    """
    Generates a 3D conformer of the OEMol
    """
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(False)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(True)
    status = omega(mol)
    if not status:
        raise(Exception("omega returned error code %d" % status))
    return mol
def check_atomtype(mol, check_type):
    for atom in mol.GetAtoms():
        if atom.GetType() == check_type:
            print("Found type %s atom in molecule %s" % (atom.GetType(), oechem.OECreateIsoSmiString(mol)))
            return False
    return True

def read_Elements(elements):
    f = open(elements, 'r')
    elements_list = f.read().splitlines()
    f.close()
    return elements_list

def keep_molecule(mol, max_heavy_atoms = 100,
        remove_smirks = list(), max_metals = 0, check_type = None):
    if oechem.OECount(mol, oechem.OEIsMetal()) > max_metals:
        return False
    if oechem.OECount(mol, oechem.OEIsHeavy()) > max_heavy_atoms:
        return False
    # 
    if oechem.OECount(mol, oechem.OEIsHeavy()) < 5:
        return False
    for smirks in remove_smirks:
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smirks):
            continue
        ss = oechem.OESubSearch(qmol)
        matches = [match for match in ss.Match(mol, False)]
        if len(matches) > 0:
            return False
    if elements != None:
        elements_list = read_Elements(elements)
        if not check_element(mol, elements_list):
            return False
    if check_type != None:
        if not check_atomtype(mol, check_type):
            return False
    return check_valence(mol)


def filter_molecules(input_molstream, output_molstream, allow_repeats = False,
        allow_warnings = False, max_heavy_atoms = 100, remove_smirks = list(),
        max_metals = 0, explicitHs = True, check_type = None,
        number = 10000, conformer = False):
    """
    Takes input file and removes molecules using given criteria then
    writes a new output file
    """
    errs = oechem.oeosstream()
    oechem.OEThrow.SetOutputStream(errs)

    molecule = oechem.OECreateOEGraphMol()
    smiles = list()

    count = 0
    warnings = 0
    smile_count = 0
    saved = 0

    while oechem.OEReadMolecule(input_molstream, molecule):
        count +=1
        if ("warning" in errs.str().lower()) and not allow_warnings:
            warnings += 1
            errs.clear()
            continue

        smi = oechem.OECreateIsoSmiString(molecule)
        mol_copy = oechem.OEMol(molecule)
        if explicitHs:
            oechem.OEAddExplicitHydrogens(mol_copy)
        new_smile = smi not in smiles
        if not new_smile:
            smile_count += 1

        if new_smile or allow_repeats:
            keep = keep_molecule(mol_copy, max_heavy_atoms, remove_smirks, max_metals, check_type)
            if keep:
                if conformer:
                    print("generating conformer for molecule %i" % count)
                    mol_copy = generate_conformer(mol_copy)
                smiles.append(smi)
                oechem.OEWriteMolecule(output_molstream, mol_copy)
                saved += 1
                if saved >= number:
                    break
        errs.clear()

    print("%i molecules in input file" % (count))
    print("%i molecules resulted in warnings when parsing" % warnings)
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to output" % (saved))

if __name__=="__main__":
    from optparse import OptionParser
    usage_string = """\
            Given a set of molecules filter for
            number of heavy atoms per molecule,
            number of metals atoms per molecule,
            inappropriate valence (carbons or nitrogens with 5+ bonds),
            SMIRKS patterns you do not wish to have included

    usage:
    filter_molecule_sets.py --input DrugBank.sdf --output updated_DrugBank.mol2.gz\
            --repeats False --warnings False --heavy 100 --SMIRKS remove_smikrs_example.smarts\
            --metals 0 --hydrogens True
    """
    parser = OptionParser(usage = usage_string)

    parser.add_option('-f', '--input', type = 'string', dest = 'input_file', default = None, action = 'store',
            help = "Name of initial molecule file")
    parser.add_option('-o', '--output', type ='string', dest = 'output_file', default = 'output_molecules.mol2.gz', action = 'store',
            help = "Name of file to save filtered Molecules")
    parser.add_option('-r', '--repeats', type = 'choice', dest = 'repeats', default = 'False',
            choices = ['True', 'False'],
            help = "If True allows for multiple molecules with the same isomeric SMILES string, default = False OPTIONAL")
    parser.add_option('-w', '--warnings', type = 'choice', dest = 'warnings', default = 'False',
            choices = ['True', 'False'],
            help = "If True keeps molecules that result in warning while loading, default = False, OPTIONAL")
    parser.add_option('-n', '--heavy', type = 'int', dest = 'heavy', default = 100,
            help = "Maximum number of heavy atoms per molecule, default = 100, OPTIONAL")
    parser.add_option('-s', '--SMIRKS', type = 'string', dest = 'smirks_file', default = None,
            help = "If not None, the file of SMIRKS are parsed and molecules that match that pattern are removed. The file should have a single column with SMIRKS/SMARTS strings where commented lines begin with a '%', default = None, OPTIONAL")
    parser.add_option('-m', '--metals', type = 'int', dest = 'metals', default = 0,
            help = "Maximum number of metals per molecule, default = 0, OPTIONAL")
    parser.add_option('-H', '--hydrogens', type = 'choice', dest = 'hydrogens', default = 'True',
            choices = ['True', 'False'],
            help = "If True, hydrogens are added to the output molecules")
    parser.add_option('-a', '--atoms', type = 'string', dest = 'atoms', default = None, action = 'store',
            help = "File name which contains the atomic number of the elements that you do not want to be in your set of molecules.")
    parser.add_option('-t', '--type', type = 'string', dest = 'atomtype', default = None,
            help = "Atom type that you do not want in your set of molecules")
    parser.add_option('-y', '--flavor', type = 'choice', dest = 'flavor', default = 'tripos',
            choices = ['tripos', 'ff'],
            help = "If you want to include different flavor - options: tripos or ff")

    parser.add_option('-z', '--number', type = 'int', dest = 'number', default = 10000, help = "Limit for the number of molecules in the final set, returns the first z molecules that meet the filtering requirements")

    parser.add_option('-C', '--conformers', type = 'choice', dest = 'conformers',
            default = 'False', choices = ['True', 'False'],
            help = "If True, generates a single conformer of each molecule before writing to the output file")

    (opt, args) = parser.parse_args()

    # Check input files
    if (opt.input_file is None) or (opt.output_file is None):
        parser.print_help()
        parser.error("Input molecules file cannot be None")

    # Load and check input file
    ifs = oechem.oemolistream(opt.input_file)
    # TODO: get input file format/flavor based on input file name?
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    if not ifs.IsValid():
        parser.print_help()
        parser.error("Error: input_file (%s) was not valid" % opt.input_file)

    # Load and check output file
    ofs = oechem.oemolostream(opt.output_file)
    if opt.flavor == 'ff':
        flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
        ofs.SetFlavor(oechem.OEFormat_MOL2, flavor)
    if not ofs.IsValid():
        parser.print_help()
        parser.error("Error: output_file (%s) was not valid" % opt.output_file)

    if opt.smirks_file is None:
        smirks = []
    else:
        smirks = smarty.utils.parse_odds_file(opt.smirks_file, False)
        smirks = smirks[0]

    # if atoms is not None add to smirks papers
    if opt.atoms is not None:
        try:
            elements = read_Elements(opt.atoms)
            elements = [ '[#%i]' % e for e in elements ]
        except:
            elementReadError = Exception("problem reading element file %s" % opt.atoms)
            raise elementReadError
        smirks += elements

    repeats = opt.repeats == 'True'
    warnings = opt.warnings == 'True'
    hydrogens = opt.hydrogens == 'True'
    conformer = opt.conformers == 'True'

    filter_molecules(ifs, ofs, repeats, warnings, opt.heavy, smirks, opt.metals, hydrogens, opt.atomtype, opt.number, conformer)

    ifs.close()
    ofs.close()
