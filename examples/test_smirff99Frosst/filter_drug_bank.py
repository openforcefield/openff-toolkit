from openeye import oechem

errs = oechem.oeosstream()
oechem.OEThrow.SetOutputStream(errs)
input_file = "DrugBank.sdf"

smiles = list()
count = 0
input_molstream = oechem.oemolistream(input_file)
input_molstream.SetFlavor(oechem.OEFormat_SDF, oechem.OEIFlavor_SDF_Default)
molecule = oechem.OECreateOEGraphMol()

molecules = list()
count = 0
warnings = 0
smile_count = 0
while oechem.OEReadMolecule(input_molstream, molecule):
    count +=1
    if "warning" in errs.str().lower():
        print("Molecule %i note parsed\n\t%s"% (count, errs.str()))
        warnings += 1
        errs.clear()
        continue

    smi = oechem.OECreateIsoSmiString(molecule)
    mol_copy = oechem.OEMol(molecule)
    heavy_atoms = oechem.OECount(mol_copy, oechem.OEIsHeavy())
    metal_atoms = oechem.OECount(mol_copy, oechem.OEIsMetal())
    oechem.OEAddExplicitHydrogens(mol_copy)
    if smi in smiles:
        print("Cannot add more than one molecule with the same isomeric SMILES string")
        smile_count += 1
    else:
        if heavy_atoms < 101:
            if metal_atoms == 0:
                smiles.append(smi)
                molecules.append(oechem.OEMol(mol_copy))

print("%i molecules found total" % count)
print("%i removed for warnings" % warnings)
print("%i removed for repeated smiles" % smile_count)
print("%i molecules saved" % len(molecules))
#for smile in smiles:
#    mol = oechem.OEGraphMol()
#    oechem.OEParseSmiles(mol, smile)
#    molecules.append(mol)

output_molstream = oechem.oemolostream("updated_DrugBank.mol2.gz")
output_molstream.SetFlavor(oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_DEFAULT)

isosmiles = list()
for mol in molecules:
    smi = oechem.OECreateIsoSmiString(mol)
    if not smi in isosmiles:
        oechem.OEWriteMolecule(output_molstream, mol)
    else:
        print("list of molecules contains more than molecule with the isomeric smiles %s" % smi)

output_molstream.close()


