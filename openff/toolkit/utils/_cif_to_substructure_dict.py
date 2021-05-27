#!/usr/bin/env python

"""
Tools for generating library of substructures from the Chemical Component
Dictionary (CCD).
"""

from openff.toolkit.topology import Molecule
from CifFile import ReadCif
from mendeleev import element
import rdkit
import openff
# from amber-ff-porting.ConvertResidueParameters import (
#     remove_charge_and_bond_order_from_imidazole,
#     remove_charge_and_bond_order_from_guanidinium
# )


class CifSubstructures:
    """
    Substructure from CIF formatted files.
    """

    # bond order dictionary translation from mmcif
    bond_order_dict = {
        'SING': 1,
        'DOUB': 2,
        'TRIP': 3,
        'QUAD': 4
    }

    def __init__(self):
        """
        Create object with substructures data from CIF files.
        """
        self.data = {}  # Dictionary where to store substructures data
        self._cif_multi_entry_object = None

    def from_file(self, cif_file, discard_keyword='FRAGMENT'):
        """
        Reads cif formatted file and fills substructure information.

        Parameters
        __________
        cif_file : str or file-like object
            File path string, URL or file-like object.
        """
        # read cif file
        self._cif_multi_entry_object = ReadCif(cif_file)
        # fill data dictionary
        self._fill_substructure_data(discard_keyword=discard_keyword)

    def to_json(self, indent=None):
        """
        Return a JSON serialized representation.

        Parameters
        __________
        indent : int, optional, default=None
            If not None, will pretty-print with specified number of spaces for indentation.

        Returns
        _______
        serialized : str
            A JSON serialized representation of the object.
        """
        return NotImplementedError

    @staticmethod
    def _generate_atom_symbol_number_dictionary(cif_entry):
        """
        Generate symbol to atomic number dictionary, given entry information from cif file.

        Parameters
        __________
        cif_entry : CifFile.StarFile.StarBlock
            Data block entry from cif file.

        Returns
        _______
        symbol_to_num_dict : dict
            Atom symbol to atomic number (key, value) dictionary.
        """
        atom_symbol_list = cif_entry['_chem_comp_atom.type_symbol']
        elements = element(atom_symbol_list)
        atomic_nums_list = []
        for element_ in elements:
            atomic_nums_list.append(element_.atomic_number)
        # Create dictionary
        symbol_to_num_dict = {symbol: number
                              for symbol, number
                              in zip(atom_symbol_list, atomic_nums_list)}
        return symbol_to_num_dict

    @staticmethod
    def _get_atom_by_name(molecule, name: str):
        """
        Gets first occurrence of atom object from a molecule by its name.

        Parameters
        __________
        molecule : Molecule
            Molecule object where to look.
        name : str
            Atom name string to match against molecule atoms.

        Returns
        _______
        atom_match : Atom or None
            First occurrence of Atom with the given name, in molecule. None if not found.
        """
        for atom in molecule.atoms:
            if atom.name == name:
                match = atom
                break
        else:
            match = None
        return match

    def _gather_atoms_information(self, cif_entry_data):
        """
        Gather minimum required atom information for creating toolkit Atoms.

        Returns
        _______
        atom_information_array : tuple
            Multidimensional tuple with atoms information.
        """
        symbol_to_num_dict = self._generate_atom_symbol_number_dictionary(cif_entry_data)
        atom_names = cif_entry_data['_chem_comp_atom.atom_id']
        atomic_numbers = [symbol_to_num_dict[x]
                          for x in cif_entry_data['_chem_comp_atom.type_symbol']]
        formal_charges = [int(x) for x in cif_entry_data['_chem_comp_atom.charge']]
        is_aromatic = [False if x == 'N' else True
                       for x in cif_entry_data['_chem_comp_atom.pdbx_aromatic_flag']]
        stereochemistry = [None if x == 'N' else x
                           for x in cif_entry_data['_chem_comp_atom.pdbx_stereo_config']]
        leaving_atoms_list = cif_entry_data['_chem_comp_atom.pdbx_leaving_atom_flag']
        return atom_names, \
            atomic_numbers,\
            formal_charges,\
            is_aromatic, \
            stereochemistry,\
            leaving_atoms_list

    def _gather_bonds_information(self, cif_entry_data):
        """
        Gather minimum required bond information for creating toolkit Bonds.

        Returns
        _______
        bond_information_array : tuple
            Multidimensional tuple with bonds information.
        """
        atom1_name_list = cif_entry_data['_chem_comp_bond.atom_id_1']
        atom2_name_list = cif_entry_data['_chem_comp_bond.atom_id_2']
        bond_order_list = [self.bond_order_dict[x]
                           for x in cif_entry_data['_chem_comp_bond.value_order']]
        is_aromatic_bond_list = [False if x == 'N' else True
                                 for x in cif_entry_data['_chem_comp_bond.pdbx_aromatic_flag']]
        stereochemistry_bond_list = [None if x == 'N' else x
                                     for x in cif_entry_data['_chem_comp_bond.pdbx_stereo_config']]
        return atom1_name_list, \
            atom2_name_list,\
            bond_order_list,\
            is_aromatic_bond_list, \
            stereochemistry_bond_list

    @staticmethod
    def _get_subrdmol(rdmol: rdkit.Chem.rdchem.Mol,
                      indices: list = [],
                      sanitize: bool = False):
        """
        Create new sub-molecule from selected atom indices

        Parameters
        ----------
        rdmol: rdkit.Chem.rdchem.Mol
            Input molecule
        indices: iterable of ints
            atom indices to include from input molecule, indexed from 0
        sanitize: bool
            whether to sanitize the molecule (recommend: no)

        Returns
        -------
        rdkit.Chem.rdchem.Mol: subset of molecule
        """
        submol = rdkit.Chem.RWMol(rdmol)
        ix = sorted([at.GetIdx() for at in rdmol.GetAtoms()
                     if at.GetIdx() not in indices])
        for i in ix[::-1]:
            submol.RemoveAtom(int(i))
        if sanitize:
            rdkit.Chem.SanitizeMol(submol)

        for atom in submol.GetAtoms():
            # print(dir(atom))
            atom.SetNoImplicit(True)

        # remove_charge_and_bond_order_from_imidazole(submol)
        # remove_charge_and_bond_order_from_guanidinium(submol)
        return submol

    def _get_smiles(self,
                   mol: openff.toolkit.topology.molecule.Molecule,
                   indices: list = [], label_indices: list = []):
        """
        Get SMARTS of selected atoms in molecule

        Parameters
        ----------
        mol: openff.toolkit.topology.molecule.Molecule
            Input molecule
        indices: iterable of ints
            atom indices to include in SMARTS, indexed from 0
        label_indices: iterable of ints
            atom indices to label, indexed from 0. The atoms
            will be labelled in the order specified. Labels
            begin from 1.

        Returns
        -------
        str: SMARTS string
        """
        rdmol = mol.to_rdkit()
        for i, lix in enumerate(label_indices, 1):
            at = rdmol.GetAtomWithIdx(int(lix))
            at.SetAtomMapNum(i)
        indices = sorted(set(indices) | set(label_indices))
        submol = self._get_subrdmol(rdmol, indices)
        smiles = rdkit.Chem.MolToSmiles(submol, allHsExplicit=True)
        smiles = smiles.replace('$', '~')
        return smiles

    def _fill_substructure_data(self, discard_keyword='FRAGMENT', include_leaving=False):
        """
        Fills data dictionary with the substructure information.

        Parameters
        __________
        discard_keyword : str, optional, default='FRAGMENT'
            Keyword in _chem_comp.name for filtering out entries. Default is 'FRAGMENT'.
        include_leaving : bool, optional, default=False
            Whether to include leaving atoms marked in _chem_comp_atom.pdbx_leaving_atom_flag
        """
        for entry, entry_data in self._cif_multi_entry_object.items():
            # create empty molecule to fill with data
            offmol = Molecule()
            # Skip entries matching discard filtering keyword in comp name
            comp_name = entry_data['_chem_comp.name']
            if discard_keyword in comp_name:
                continue
            # Gather necessary data for creating/adding atoms
            atoms_information = self._gather_atoms_information(entry_data)
            atom_names = atoms_information[0]
            atomic_numbers = atoms_information[1]
            formal_charges = atoms_information[2]
            are_aromatic = atoms_information[3]
            stereochemistry = atoms_information[4]
            leaving_atoms_list = atoms_information[5]
            # add atoms
            for atom_idx in range(len(atom_names)):
                if include_leaving:
                    # Add all atoms
                    offmol.add_atom(atomic_numbers[atom_idx],
                                    formal_charges[atom_idx],
                                    are_aromatic[atom_idx],
                                    stereochemistry=stereochemistry[atom_idx],
                                    name=atom_names[atom_idx]
                                    )
                else:
                    # Add only not leaving atoms
                    if leaving_atoms_list[atom_idx] == 'N':
                        offmol.add_atom(atomic_numbers[atom_idx],
                                        formal_charges[atom_idx],
                                        are_aromatic[atom_idx],
                                        stereochemistry=stereochemistry[atom_idx],
                                        name=atom_names[atom_idx]
                                        )
            # Gather information for bonds
            bonds_information = self._gather_bonds_information(entry_data)
            atom1_name_list = bonds_information[0]
            atom2_name_list = bonds_information[1]
            bond_order_list = bonds_information[2]
            is_aromatic_bond_list = bonds_information[3]
            stereochemistry_bond_list = bonds_information[4]
            # add bonds
            # TODO: What about fractional bond order?
            for bond_idx in range(len(atom1_name_list)):
                try:
                    offmol.add_bond(self._get_atom_by_name(offmol, atom1_name_list[bond_idx]),
                                    self._get_atom_by_name(offmol, atom2_name_list[bond_idx]),
                                    bond_order_list[bond_idx],
                                    is_aromatic_bond_list[bond_idx],
                                    stereochemistry=stereochemistry_bond_list[bond_idx]
                                    )
                # TODO: Broad exception. We should raise a specific error.
                except Exception:  # happens when atom is not found, i.e. no leaving atoms.
                    continue
            # Get smiles
            atom_indices = list(range(offmol.n_atoms))
            smiles = self._get_smiles(offmol,
                                      indices=atom_indices,
                                      label_indices=atom_indices)
            self.data[entry] = smiles
