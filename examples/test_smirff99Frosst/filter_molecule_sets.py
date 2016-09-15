from openeye import oechem

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

def filter_molecules(input_file, input_format, input_flavor, output_file,
        output_format, output_flavor, allow_repeats = False, allow_warnings = False,
        max_heavy_atoms = 100, remove_smirks = list(), max_metals = 0, explicitHs = True):
    """
    Takes input file and removes molecules using given criteria then
    writes a new output file
    """
    errs = oechem.oeosstream()
    oechem.OEThrow.SetOutputStream(errs)

    input_molstream = oechem.oemolistream(input_file)
    input_molstream.SetFlavor(input_format, input_flavor)
    output_molstream = oechem.oemolostream(output_file)
    output_molstream.SetFlavor(output_format, output_flavor)

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
    input_molstream.close()
    output_molstream.close()

    print("%i molecules in %s" % (count, input_file))
    print("%i molecules resulted in warnings when parsing" % warnings)
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to %s" % (saved, output_file))

if __name__=="__main__":
    from optparse import OptionParser
    usage_string = """\
    TODO: add usage string
    """
    parser = OptionParser(usage = usage_string)

    parser.add_option('-f', '--input', type = 'string', dest = 'input_file', default = None, action = 'store', help = "Name of initial molecule file")
    parser.add_option('-o', '--output', type ='string', dest = 'output_file', default = 'output_molecules.mol2.gz', action = 'store', help = "Name of file to save filtered Molecules")

    (opt, args) = parser.parse_args()
    if (opt.input_file is None) or (opt.output_file is None):
        parser.print_help()
        parser.error("Input molecules file cannot be None")
    filter_molecules(opt.input_file, oechem.OEFormat_SDF, oechem.OEIFlavor_SDF_Default, opt.output_file, oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_DEFAULT)

