from openeye import oechem

errs = oechem.oeosstream()
oechem.OEThrow.SetOutputStream(errs)
input_file = "DrugBank.sdf"

smiles = set()
count = 0
input_molstream = oechem.oemolistream(input_file)
input_molstream.SetFlavor(oechem.OEFormat_SDF, oechem.OEIFlavor_SDF_Default)
molecule = oechem.OECreateOEGraphMol()

errs = oechem.oeosstream()
oechem.OEThrow.SetOutputStream(errs)
molecules = list()

while oechem.OEReadMolecule(input_molstream, molecule):
    errs.clear()
    name = oechem.OEGetSDData(molecule, 'name').strip()
    molecule.SetTitle(name)
    if "Warning" in errs.str():
        count += 1
        print("%i Molecules not parsed" % count)
        continue 

    mol_copy = oechem.OEMol(molecule)
    heavy_atoms = oechem.OECount(mol_copy, oechem.OEIsHeavy())
    metal_atoms = oechem.OECount(mol_copy, oechem.OEIsMetal())
    oechem.OESuppressHydrogens(mol_copy)
    smi = oechem.OECreateIsoSmiString(mol_copy)
    if not smi in smiles:
        if heavy_atoms < 101:
            if metal_atoms == 0:
                smiles.add(smi)

for smile in smiles:
    mol = oechem.OEGraphMol()
    oechem.OEParseSmiles(mol, smile)
    molecules.append(mol)
    
output_molstream = oechem.oemolostream("updated_DrugBank.mol2.gz")
output_molstream.SetFlavor(oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_DEFAULT)

for mol in molecules:
    oechem.OEWriteMolecule(output_molstream, mol)

output_molstream.close()
