#!/usr/bin/env python

"""
Tools for generating library of substructures from the Chemical Component
Dictionary (CCD).
"""

import copy
from collections import defaultdict
from openff.toolkit.topology import Molecule
from CifFile import ReadCif
import json
from mendeleev import element
import rdkit
import openff


# Amber ff porting functions
def remove_charge_and_bond_order_from_guanidinium(rdmol):
    """
    To correct for chemical perception issues with possible resonance states of arginine,
    remove all charge from the guanidinium group, and set all bond orders to 4. This will
    mark the resonant bonds with a unique "$" character in the SMARTS, which we can later
    replace.
    """
    for atom in rdmol.GetAtoms():
        # print(dir(atom))
        if atom.GetAtomicNum() != 6:  # element.symbol != "C":
            continue
        nitrogen_neighbors = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:
                nitrogen_neighbors += 1
        if nitrogen_neighbors != 3:
            continue
        atom.SetFormalCharge(0)
        for neighbor in atom.GetNeighbors():
            neighbor.SetFormalCharge(0)
        for bond in atom.GetBonds():
            # Set bond order 4, which will produce a "$" character. We later replace this with "~".
            # print(dir(bond))
            bond.SetBondType(rdkit.Chem.BondType.QUADRUPLE)
            # bond.SetBondType(Chem.BondType.AROMATIC)
            # bond.SetBondType(Chem.BondType.OTHER)


def remove_charge_and_bond_order_from_imidazole(rdmol):
    """
    To correct for chemical perception issues with possible resonance states of histidine,
    remove all charge from the imidazole group, and set all bond orders to 4. This will
    mark the resonant bonds with a unique "$" character in the SMARTS, which we can later
    replace.
    """
    # matches = offmol.chemical_environment_matches('[C:1]1~[C:2]~[N:3]~[C:4]~[N:5]1')
    imidazole_substructure = rdkit.Chem.MolFromSmarts('[C:1]1~[C:2]~[N:3]~[C:4]~[N:5]1')
    matches = rdmol.GetSubstructMatches(imidazole_substructure)
    all_imidazole_atoms = set()
    for match in matches:
        for idx in match:
            all_imidazole_atoms.add(idx)

    for atom in rdmol.GetAtoms():
        if atom.GetIdx() in all_imidazole_atoms:
            atom.SetFormalCharge(0)

    for bond in rdmol.GetBonds():
        # print(dir(bond))
        if ((bond.GetBeginAtomIdx() in all_imidazole_atoms) and
                (bond.GetEndAtomIdx() in all_imidazole_atoms)):
            bond.SetBondType(rdkit.Chem.BondType.QUADRUPLE)


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
        self.data = defaultdict(defaultdict)  # Dictionary where to store substructures data
        self._cif_multi_entry_object = None

    def from_file(self,
                  cif_file,
                  include_leaving=False,
                  discard_keyword='FRAGMENT',
                  replace_quadruple_bond_with_any: bool = True,
                  remove_charge_bond_order_resonant: bool = True
                  ):
        """
        Reads cif formatted file and fills substructure information.

        Sequential runs combine information from different files.

        Parameters
        __________
        cif_file : str or file-like object
            File path string, URL or file-like object.
        include_leaving : bool, optional, default=False
            Whether to include leaving atoms marked in _chem_comp_atom.pdbx_leaving_atom_flag
        discard_keyword : str, optional, default='FRAGMENT'
            Keyword in _chem_comp.name for filtering out entries. Default is 'FRAGMENT'.
        remove_charge_bond_order_resonant: bool, optional, default=True
            Removes bond order and charge information from guanidium and imidazole resonant groups.
        replace_quadruple_bond_with_any: bool, optional, default=True:
            Replace expected unique quadruple bond symbol ($) with any order symbol (~) in smarts patterns.
        """
        # read cif file
        self._cif_multi_entry_object = ReadCif(cif_file)
        # fill data dictionary
        self._fill_substructure_data(include_leaving=include_leaving,
                                     discard_keyword=discard_keyword,
                                     replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
                                     remove_charge_bond_order_resonant=remove_charge_bond_order_resonant
                                     )

    def to_json_file(self, output_file, indent=4):
        """
        Return a JSON serialized representation.

        Parameters
        __________
        output_file : str
            Path string for output file.
        indent : int, optional, default=4
            If not None, will pretty-print with specified number of spaces for indentation.

        Returns
        _______
        serialized : str
            A JSON serialized representation of the object.
        """
        with open(output_file, 'w') as library_file:
            json.dump(self.data, library_file, indent=indent)

    def clear_data(self):
        """
        Clears current data structures read from files.

        Useful if you want to start reading cif files from scratch.
        """
        self.data = defaultdict(defaultdict)
        self._cif_multi_entry_object = None

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
                      sanitize: bool = False,
                      remove_charge_bond_order_resonant: bool = True
                      ):
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

        if remove_charge_bond_order_resonant:
            remove_charge_and_bond_order_from_imidazole(submol)
            remove_charge_and_bond_order_from_guanidinium(submol)
        return submol

    def _get_smiles(self,
                    mol: openff.toolkit.topology.molecule.Molecule,
                    indices: list = [],
                    label_indices: list = [],
                    replace_quadruple_bond_with_any: bool = True,
                    remove_charge_bond_order_resonant: bool = True,
                    ):
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
        submol = self._get_subrdmol(rdmol,
                                    indices,
                                    remove_charge_bond_order_resonant=remove_charge_bond_order_resonant
                                    )
        smiles = rdkit.Chem.MolToSmiles(submol, allHsExplicit=True)
        # TODO: This doesn't seem to be having an effect. It is replacing bond orders regardless.
        if replace_quadruple_bond_with_any:
            smiles = smiles.replace('$', '~')
        return smiles

    def _add_substructure_data_entry(self,
                                     entry_data,
                                     include_leaving,
                                     discard_keyword,
                                     additional_leaving=None,
                                     replace_quadruple_bond_with_any: bool = True,
                                     remove_charge_bond_order_resonant: bool = True
                                     ):
        """
        Read substructure from entry data in cif file and fill data object with smiles.

        Parameters
        __________
        entry_data : CifFile.StarFile.StarBlock
            Cif file entry data object from pycifrw module.
        include_leaving : bool, optional, default=False
            Whether to include leaving atoms marked in _chem_comp_atom.pdbx_leaving_atom_flag
        discard_keyword : str, optional, default='FRAGMENT'
            Keyword in _chem_comp.name for filtering out entries. Default is 'FRAGMENT'.
        additional_leaving : None or List, default=None
            Additional list of atoms in entry to treat as leaving atoms.
        """
        offmol = Molecule()
        # Skip entries matching discard filtering keyword in comp name
        comp_name = entry_data['_chem_comp.name']
        if discard_keyword in comp_name:
            return
        # Gather necessary data for creating/adding atoms
        atoms_information = self._gather_atoms_information(entry_data)
        atom_names_orig = atoms_information[0]
        n_atoms = len(atom_names_orig)  # number of atoms in entry
        atomic_numbers = atoms_information[1]
        formal_charges = atoms_information[2]
        are_aromatic = atoms_information[3]
        stereochemistry = atoms_information[4]
        leaving_atoms_list = atoms_information[5]
        if additional_leaving:
            for atom_idx, atom_name in enumerate(atom_names_orig):
                if atom_name in additional_leaving:
                    leaving_atoms_list[atom_idx] = 'Y'
        # copy atom names to object that can be modified if not leaving atoms
        atom_names = copy.copy(atom_names_orig)
        # add atoms
        for atom_idx in range(len(atom_names_orig)):
            if include_leaving:
                # Add all atoms
                offmol.add_atom(atomic_numbers[atom_idx],
                                formal_charges[atom_idx],
                                are_aromatic[atom_idx],
                                stereochemistry=stereochemistry[atom_idx],
                                name=atom_names_orig[atom_idx]
                                )
            else:
                # Add only not leaving atoms
                if leaving_atoms_list[atom_idx] == 'N':
                    offmol.add_atom(atomic_numbers[atom_idx],
                                    formal_charges[atom_idx],
                                    are_aromatic[atom_idx],
                                    stereochemistry=stereochemistry[atom_idx],
                                    name=atom_names_orig[atom_idx]
                                    )
                else:
                    # Remove original leaving atom name from current atom names
                    atom_names.remove(atom_names_orig[atom_idx])
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
                                  label_indices=atom_indices,
                                  replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
                                  remove_charge_bond_order_resonant=remove_charge_bond_order_resonant
                                  )
        # Only take three letter code for key -- Note it uses the new atom_names list
        entry_code = entry_data['_chem_comp.three_letter_code']
        if smiles not in self.data[entry_code]:
            self.data[entry_code][smiles] = atom_names

    def _fill_substructure_data(self,
                                include_leaving=False,
                                discard_keyword='FRAGMENT',
                                replace_quadruple_bond_with_any: bool = True,
                                remove_charge_bond_order_resonant: bool = True
                                ):
        """
        Fills data dictionary with the substructure information.

        Parameters
        __________
        discard_keyword : str, optional, default='FRAGMENT'
            Keyword in _chem_comp.name for filtering out entries. Default is 'FRAGMENT'.
        include_leaving : bool, optional, default=False
            Whether to include leaving atoms marked in _chem_comp_atom.pdbx_leaving_atom_flag
        """
        override_dict = {
            #LYN
            'lys_lfoh_dhz3': ['H2'],
            #GLU
            'glu_lfoh_dhe2': ['H2'],
            #HIE
            'his_lfoh_dhd1': ['H2'],
            #ASP
            'asp_lfoh_dhd2': ['H2'],
            #HID
            'his_lfoh_dhe2': ['H2'],
            #CYX
            'cys': ['HG']
        }
        for entry, entry_data in self._cif_multi_entry_object.items():
            # create empty molecule to fill with data
            # if 'LYS' in entry:
            #     continue
            self._add_substructure_data_entry(entry_data,
                                              include_leaving,
                                              discard_keyword,
                                              replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
                                              remove_charge_bond_order_resonant=remove_charge_bond_order_resonant
                                              )
            if entry in override_dict:
                self._add_substructure_data_entry(entry_data,
                                                  include_leaving,
                                                  discard_keyword,
                                                  override_dict[entry],
                                                  replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
                                                  remove_charge_bond_order_resonant=remove_charge_bond_order_resonant
                                                  )

    def _patch_known_problems(self):
        """
        Monkey-patching known problems with current aa-variants-v1.cif from CCD.

        .. warning: Needed as of Oct-06-2021
        """
        # Fix PRO smarts substructure
        old_smarts = "[N-:1]1[C@@:2]([C-:3]=[O:4])([H:8])[C:5]([H:9])([H:10])[C:6]([H:11])([H:12])[C:7]1([H:13])[H:14]"
        #new_smarts = "[N:1]1[C@@:2]([C:3]=[O:4])([H:8])[C:5]([H:9])([H:10])[C:6]([H:11])([H:12])[C:7]1([H:13])[H:14]"
        #self.data['PRO'][new_smarts] = self.data['PRO'][old_smarts]
        self.data['PRO'].pop(old_smarts)

        # Add common caps
        self.data['ACE'] = {
            "[C:1](=[O:2])[C:3]([H:4])([H:5])[H:6]": [
                "C",
                "O",
                "CH3",
                "H1",
                "H2",
                "H3"
            ]
        }
        self.data["NME"] = {
            "[N:1]([C:2]([H:4])([H:5])[H:6])[H:3]": [
                "N",
                "C",
                "HN2",
                "H1",
                "H2",
                "H3"
            ]
        }
        self.data["NH2"] = {
            "[N:1]([H:2])[H:3]": [
                "N",
                "HN1",
                "HN2"
            ]
        }
