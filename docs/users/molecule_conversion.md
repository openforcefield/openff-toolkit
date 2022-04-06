# Molecule conversion to other packages


Molecule conversion spec

<!---
Add these sections later

## Chemical data

## Representation differences
* Local vs. global stereo
* No conformers or all-zero coords?

## Sources of "trust"
* Graph stereo is trusted over 3D stereo

## Expectations for imports
* RDKit and OpenEye
    * Stereo perceived
    * Explicit Hs
    * 
--->



## Hierarchy data (chains and residues)

Note that the representations of hierarchy data (namely, chains and residues) in 
different software packages have different expectations. 
For example, in OpenMM, the atoms of a single residue must be contiguous. 
In most packages, it is expected that any atom with a residue name defined will also have a residue number. 

The OpenFF toolkit does not have these restrictions, and records hierarchy 
metadata almost entirely for interoperability and user convenience reasons. 
No code paths in the OpenFF Toolkit consider hierarchy metadata during parameter assignment. 
While users should expect hierarchy metadata to be correctly handled in
simple loading operations and export to other packages, _modifying_ 
hierarchy metadata in the OpenFF Toolkit may lead to unexpected incompatibilities with other packages/representations.  

Another consequence of this difference in representations is that hierarchy iterators (like 
`Molecule.residues` and `Topology.chains`) are not accessed during conversion to other packages. 
Only the underlying hierarchy metadata from the atoms is transferred, and the OpenFF Toolkit makes
no attempt to match the behavior of iterators in other packages.

In cases where only _some_ common metadata fields are set (but not others), the following calls 
happen during conversion TO other packages

* RDKit - We always run `rdatom.SetPDBMetadata`, meaning that undefined values in an OFFMol are set to be the default values in RDKit
* OpenEye - We always run `oechem.OEAtomSetResidue(oe_atom, res)`, meaning that undefined values in an OFFMol will be set by OpenEye.
* OpenMM - ???

| Toolkit                      | residue_name                                                                           | residue_number                                                                                  | chain_id                                                                                  |
|------------------------------|----------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| PDB file ATOM/HETATM columns | Columns 18-20 (resName)                                                                | Columns 23-26 (resSeq)                                                                          | Columns 22 (chainID)                                                                      |
| PDBx/MMCIF fields            | label_comp_id                                                                          | label_seq_id                                                                                    | label_asym_id                                                                             |
| OpenFF getter (defined)      | atom.metadata[‘residue_name’]                                                          | atom.metadata[‘residue_number’]                                                                 | atom.metadata[‘chain_id’]                                                                 |
| OpenFF getter (undefined)    | ‘residue_name’ not in atom.metadata                                                    | ‘residue_number’ not in atom.metadata                                                           | ‘chain_id’ not in atom.metadata                                                           |
| OpenFF setter (defined)      | Defined: atom.metadata[‘residue_name’] = X                                             | atom.metadata[‘residue_number’] = X                                                             | atom.metadata[‘chain_id’] = X                                                             |
| OpenFF setter (undefined)    | del atom.metadata[‘residue_name’]                                                      | del atom.metadata[‘residue_number’]                                                             |  del atom.metadata[‘chain_id’]                                                            |
| OpenMM getter (defined)      | omm_atom.residue.name                                                                  | omm_atom.residue.id                                                                             | omm_atom.residue.chain.id                                                                 |
| OpenMM getter (undefined)    | Not possible                                                                           | Not possible                                                                                    | Not possible                                                                              |
| OpenMM setter (defined)      | omm_atom.residue.name = X                                                              | omm_atom.residue.id = X                                                                         | omm_atom.residue.chain.id = X                                                             |
| OpenMM setter (undefined)    | omm_atom.residue.name = “UNL”                                                          | omm_atom.residue.id = 0                                                                         | omm_atom.residue.chain.id = “X”                                                           |
| RDKit getter (defined)       | rda.GetPDBResidueInfo().GetResidueName()                                               | rda.GetPDBResidueInfo().GetResidueNumber()                                                      | rda.GetPDBResidueInfo().GetChainId()                                                      |
| RDKit getter (undefined)     | rda.GetPDBResidueInfo() is None                                                        | rda.GetPDBResidueInfo() is None                                                                 | rda.GetPDBResidueInfo() is None                                                           |
| RDKit setter (defined)       | rda.GetPDBResidueInfo().SetResidueName(X)                                              | rda.GetPDBResidueInfo().SetResidueNumber(X)                                                     | rda.GetPDBResidueInfo().SetChainId(X)                                                     |
| RDKit setter (undefined)     | No action                                                                              | No action                                                                                       | No action                                                                                 |
| OpenEye getter (defined)     | oechem.OEAtomGetResidue(atom).GetName()                                                | oechem.OEAtomGetResidue(atom).GetResidueNumber()                                                | oechem.OEAtomGetResidue(atom).GetChainID()                                                |
| OpenEye getter (undefined)   | oechem.OEHasResidues(oemol) == False                                                   | oechem.OEHasResidues(oemol) == False:                                                           | oechem.OEHasResidues(oemol) == False                                                      |
| OpenEye setter (defined)     | res = oechem.OEAtomGetResidue(atom) res.SetName(X) oechem.OEAtomSetResidue(atom, res)  | res = oechem.OEAtomGetResidue(atom) res.SetResidueNumber(X) oechem.OEAtomSetResidue(atom, res)  | res = oechem.OEAtomGetResidue(atom) res.SetChainID(X) oechem.OEAtomSetResidue(atom, res)  |
| OpenEye setter (undefined)   | No action                                                                              | No action                                                                                       | No action                                                                                 |