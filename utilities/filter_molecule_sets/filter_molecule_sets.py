from openeye import oechem
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

def keep_molecule(mol, max_heavy_atoms = 100,
        remove_smirks = list(), max_metals = 0):
    if oechem.OECount(mol, oechem.OEIsMetal()) > max_metals:
        return False
    if oechem.OECount(mol, oechem.OEIsHeavy()) > max_heavy_atoms:
        return False
    for smirks in remove_smirks:
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smirks):
            continue
        ss = oechem.OESubSearch(qmol)
        matches = [match for match in ss.Match(mol, False)]
        if len(matches) > 0:
            return False
    return check_valence(mol)

def filter_molecules(input_molstream, output_molstream, allow_repeats = False,
        allow_warnings = False, max_heavy_atoms = 100, remove_smirks = list(),
        max_metals = 0, explicitHs = True):
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
            keep = keep_molecule(mol_copy, max_heavy_atoms, remove_smirks, max_metals)
            if keep:
                smiles.append(smi)
                oechem.OEWriteMolecule(output_molstream, mol_copy)
                saved += 1

        errs.clear()

    print("%i molecules in %s" % (count, input_file))
    print("%i molecules resulted in warnings when parsing" % warnings)
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to %s" % (saved, output_file))

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

    (opt, args) = parser.parse_args()

    # Check input files
    if (opt.input_file is None) or (opt.output_file is None):
        parser.print_help()
        parser.error("Input molecules file cannot be None")

    # Load and check input file
    ifs = oechem.oemolistream(opt.input_file)
    if not ifs.IsValid():
        parser.print_help()
        parser.error("Error: input_file (%s) was not valid" % opt.input_file)

    # Load and check output file
    ofs = oechem.oemolostream(opt.output_file)
    if not ofs.IsValid():
        parser.print_help()
        parser.error("Error: output_file (%s) was not valid" % opt.output_file)

    smirks = smarty.utils.parse_odds_file(opt.smirks_file, False)
    smirks = smirks[0]

    repeats = opt.repeats == 'True'
    warnings = opt.warnings == 'True'
    hydrogens = opt.hydogens == 'True'

    filter_molecules(ifs, ofs, repeats, warnings, opt.heavy, smirks, opt.metals, hydrogens)

    ifs.close()
    ofs.close()
