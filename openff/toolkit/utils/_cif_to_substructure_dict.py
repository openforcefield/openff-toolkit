# flake8: noqa
"""
Tools for generating library of substructures from the Chemical Component
Dictionary (CCD).
"""

import copy
import json
from collections import defaultdict

import rdkit
from CifFile import ReadCif
from openff.units.elements import SYMBOLS

import openff
from openff.toolkit.topology import Molecule

_SYMBOL_TO_ATOMIC_NUMBER = {v: k for k, v in SYMBOLS.items()}

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
    imidazole_substructure = rdkit.Chem.MolFromSmarts("[C:1]1~[C:2]~[N:3]~[C:4]~[N:5]1")
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
        if (bond.GetBeginAtomIdx() in all_imidazole_atoms) and (
            bond.GetEndAtomIdx() in all_imidazole_atoms
        ):
            bond.SetBondType(rdkit.Chem.BondType.QUADRUPLE)


class CifSubstructures:
    """
    Substructure from CIF formatted files.
    """

    # bond order dictionary translation from mmcif
    bond_order_dict = {"SING": 1, "DOUB": 2, "TRIP": 3, "QUAD": 4}

    def __init__(self):
        """
        Create object with substructures data from CIF files.
        """
        self.data = defaultdict(
            defaultdict
        )  # Dictionary where to store substructures data
        self._cif_multi_entry_object = None

    def from_file(
        self,
        cif_file,
        discard_keyword="FRAGMENT",
        replace_quadruple_bond_with_any: bool = True,
        remove_charge_bond_order_resonant: bool = True,
    ):
        """
        Reads cif formatted file and fills substructure information.

        Sequential runs combine information from different files.

        Parameters
        __________
        cif_file : str or file-like object
            File path string, URL or file-like object.
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
        self._fill_substructure_data(
            discard_keyword=discard_keyword,
            replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
            remove_charge_bond_order_resonant=remove_charge_bond_order_resonant,
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
        with open(output_file, "w") as library_file:
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
        # Is this de-duplicated? i.e. can it be ["C", "C", "H"] or only ["C", "H"] ?
        atom_symbol_list: List[str] = cif_entry["_chem_comp_atom.type_symbol"]

        # Is this dict a de-duplicated mapping between element symbols and atomic numbers
        # for the elements in this cif file? If so, can maybe do it in one step
        """
        symbol_to_num_dict: Dict[str, int] = {
            symbol: atomic_number for atomic_number, symbol in SYMBOLS.items() if atomic_number in
            cif_entry["_chem_comp_atom.type_symbol"]
        }
        """

        symbol_to_num_dict: Dict[str, int] = {
            symbol: _SYMBOL_TO_ATOMIC_NUMBER[symbol] for symbol in atom_symbol_list
        }

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
        symbol_to_num_dict = self._generate_atom_symbol_number_dictionary(
            cif_entry_data
        )
        atom_names = cif_entry_data["_chem_comp_atom.atom_id"]
        atomic_numbers = [
            symbol_to_num_dict[x] for x in cif_entry_data["_chem_comp_atom.type_symbol"]
        ]
        formal_charges = [int(x) for x in cif_entry_data["_chem_comp_atom.charge"]]
        is_aromatic = [
            False if x == "N" else True
            for x in cif_entry_data["_chem_comp_atom.pdbx_aromatic_flag"]
        ]
        stereochemistry = [
            None if x == "N" else x
            for x in cif_entry_data["_chem_comp_atom.pdbx_stereo_config"]
        ]
        leaving_atoms_list = cif_entry_data["_chem_comp_atom.pdbx_leaving_atom_flag"]
        return (
            atom_names,
            atomic_numbers,
            formal_charges,
            is_aromatic,
            stereochemistry,
            leaving_atoms_list,
        )

    def _gather_bonds_information(self, cif_entry_data):
        """
        Gather minimum required bond information for creating toolkit Bonds.

        Returns
        _______
        bond_information_array : tuple
            Multidimensional tuple with bonds information.
        """
        atom1_name_list = cif_entry_data["_chem_comp_bond.atom_id_1"]
        atom2_name_list = cif_entry_data["_chem_comp_bond.atom_id_2"]
        bond_order_list = [
            self.bond_order_dict[x]
            for x in cif_entry_data["_chem_comp_bond.value_order"]
        ]
        is_aromatic_bond_list = [
            False if x == "N" else True
            for x in cif_entry_data["_chem_comp_bond.pdbx_aromatic_flag"]
        ]
        stereochemistry_bond_list = [
            None if x == "N" else x
            for x in cif_entry_data["_chem_comp_bond.pdbx_stereo_config"]
        ]
        return (
            atom1_name_list,
            atom2_name_list,
            bond_order_list,
            is_aromatic_bond_list,
            stereochemistry_bond_list,
        )

    @staticmethod
    def _get_subrdmol(
        rdmol: rdkit.Chem.rdchem.Mol,
        indices: list = [],
        sanitize: bool = False,
        remove_charge_bond_order_resonant: bool = True,
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
        ix = sorted(
            [at.GetIdx() for at in rdmol.GetAtoms() if at.GetIdx() not in indices]
        )
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

    def _get_smiles(
        self,
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
        submol = self._get_subrdmol(
            rdmol,
            indices,
            remove_charge_bond_order_resonant=remove_charge_bond_order_resonant,
        )
        smiles = rdkit.Chem.MolToSmiles(submol, allHsExplicit=True)
        # TODO: This doesn't seem to be having an effect. It is replacing bond orders regardless.
        if replace_quadruple_bond_with_any:
            smiles = smiles.replace("$", "~")
        return smiles

    def _add_substructure_data_entry(
        self,
        atoms_information,
        bonds_information,
        entry_code,
        replace_quadruple_bond_with_any: bool = True,
        remove_charge_bond_order_resonant: bool = True,
    ):
        """
        Read substructure from entry data in cif file and fill data object with smiles.

        Parameters
        __________
        atoms_information : list[tuple[str,int]]
            Iterable of tuple of atom info.
            Currently [atom_name, element, formal_charge, is_aromatic, stereo, is_leaving]
        bonds_information : list[tuple[str,int]]
            Iterable of tuple of bond info.
            Currently [atom1_name, atom2_name, bond_order, is_aromatic, stereo]
        """
        offmol = Molecule()

        atom_name_to_idx = dict()
        # add atoms
        for (
            atom_name,
            atomic_number,
            formal_charge,
            is_aromatic,
            stereochemistry,
            is_leaving,
        ) in atoms_information:
            # if include_leaving:
            # Add all atoms
            atom_idx = offmol.add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                stereochemistry=stereochemistry,
                name=atom_name,
            )
            atom_name_to_idx[atom_name] = atom_idx
        for (
            atom1_name,
            atom2_name,
            bond_order,
            is_aromatic,
            stereochemistry,
        ) in bonds_information:
            # try:
            atom1_idx = atom_name_to_idx[atom1_name]
            atom2_idx = atom_name_to_idx[atom2_name]
            offmol.add_bond(
                atom1_idx,
                atom2_idx,
                bond_order,
                is_aromatic=is_aromatic,
                stereochemistry=stereochemistry,
            )
            # # TODO: Broad exception. We should raise a specific error.
            # except Exception:  # happens when atom is not found, i.e. no leaving atoms.
            #     continue
        # Get smiles
        atom_indices = list(range(offmol.n_atoms))
        smiles = self._get_smiles(
            offmol,
            indices=atom_indices,
            label_indices=atom_indices,
            replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
            remove_charge_bond_order_resonant=remove_charge_bond_order_resonant,
        )
        if smiles not in self.data[entry_code]:
            atom_names = [atom_info[0] for atom_info in atoms_information]
            self.data[entry_code][smiles] = atom_names

    def _fill_substructure_data(
        self,
        discard_keyword="FRAGMENT",
        replace_quadruple_bond_with_any: bool = True,
        remove_charge_bond_order_resonant: bool = True,
    ):
        """
        Fills data dictionary with the substructure information.

        Parameters
        __________
        discard_keyword : str, optional, default='FRAGMENT'
            Keyword in _chem_comp.name for filtering out entries. Default is 'FRAGMENT'.
        """
        # override_dict contains atoms that we DO want to treat as leaving, even
        # though the actual entries don't have them flagged as leaving. The "H2" entries
        # in this list are hydrogens attached to the backbone N, which are (for some reason)
        # recorded as NH2, with neither H marked as "leaving", which causes the
        # substructures to be unable to match main-chain appearances of the amino acids.
        override_dict = {
            # LYN
            "lys_lfoh_dhz3": ["H2"],
            # GLU
            "glu_lfoh_dhe2": ["H2"],
            # HIE
            "his_lfoh_dhd1": ["H2"],
            # ASP
            "asp_lfoh_dhd2": ["H2"],
            # HID
            "his_lfoh_dhe2": ["H2"],
            # CYX
            "cys": ["HG"],
        }

        for entry, entry_data in self._cif_multi_entry_object.items():
            # Skip entries matching discard filtering keyword in comp name
            comp_name = entry_data["_chem_comp.name"]
            if discard_keyword in comp_name:
                continue
            # Gather necessary data for creating/adding atoms
            atoms_information = [*zip(*self._gather_atoms_information(entry_data))]
            # Gather information for bonds
            bonds_information = [*zip(*self._gather_bonds_information(entry_data))]

            for override_list in override_dict.get(entry, []):
                for idx, atom in enumerate(atoms_information):
                    if atom[0] in override_list:
                        new_atom_tuple = (
                            atom[0],
                            atom[1],
                            atom[2],
                            atom[3],
                            atom[4],
                            "Y",
                        )
                        atoms_information[idx] = new_atom_tuple
            # Only take three letter code for key -- Note it uses the new atom_names list
            entry_code = entry_data["_chem_comp.three_letter_code"]
            # Generate structures with permutations of leaving atoms. For example, if a COH
            # motif has -OH group marked as leaving, this will add
            # 1) the O- intermediate where only the H has left, as well as
            # 2) the structure with the entire OH removed.
            # Only atoms marked as "leaving" may have their formal charge changed by
            # this function, so the remaining C when the OH is removed still has a
            # neutral charge.
            info_tuples = self._recursive_prepare_atom_bond_info(
                atoms_information, bonds_information, return_list=list()
            )

            # Also add a structure with ALL leaving atoms removed
            atoms_info_copy = copy.deepcopy(atoms_information)
            bonds_info_copy = copy.deepcopy(bonds_information)
            atom_idxs_to_remove = [
                idx for idx, at in enumerate(atoms_info_copy) if at[5] == "Y"
            ]
            atom_names_to_remove = [
                atoms_info_copy[idx][0] for idx in atom_idxs_to_remove
            ]
            bond_idxs_to_remove = list()
            for bond_idx, bond in enumerate(bonds_info_copy):
                if (bond[0] in atom_names_to_remove) or (
                    bond[1] in atom_names_to_remove
                ):
                    bond_idxs_to_remove.append(bond_idx)

            for atom_idx_to_remove in atom_idxs_to_remove[::-1]:
                atoms_info_copy.pop(atom_idx_to_remove)
            for bond_idx_to_remove in bond_idxs_to_remove[::-1]:
                bonds_info_copy.pop(bond_idx_to_remove)
            info_tuples.append((atoms_info_copy, bonds_info_copy))

            # Generate substructure entries for all prepared substructures
            for atoms_information, bonds_information in info_tuples:
                self._add_substructure_data_entry(
                    atoms_information,
                    bonds_information,
                    entry_code,
                    replace_quadruple_bond_with_any=replace_quadruple_bond_with_any,
                    remove_charge_bond_order_resonant=remove_charge_bond_order_resonant,
                )

    def _recursive_prepare_atom_bond_info(
        self, atoms_info, bonds_info, return_list=list()
    ):
        """
        This method makes permutations of the substructure with leaving atoms left and removed.
        If some leaving atoms are bound to other leaving atoms, this method produces all permutations
        of the leaving atom chains that are possible without creating a disconnected molecular graph.

        If multiple leaving atoms form a chain, then this method attempts to keep the graph chemically
        valid by subtracting one from the formal charge of each leaving atom when its neighbor is removed.
        However non-leaving atoms will never have their formal charge changed.

        Parameters
        ----------
        atoms_info : list[tuple[str,int]]
            Iterable of tuple of atom info.
            Currently [atom_name, element, formal_charge, is_aromatic, stereo, is_leaving]
        bonds_info : list[tuple[str,int]]
            Iterable of tuple of bond info.
            Currently [atom1_name, atom2_name, bond_order, is_aromatic, stereo]
        """
        return_list.append((atoms_info, bonds_info))
        for leaving_atom_idx, leaving_atom_info_tuple in enumerate(atoms_info):
            # If the atom isn't leaving then we don't need to do anything
            if leaving_atom_info_tuple[5] == "N":
                continue
            # If the atom IS leaving, then find any bonds it's involved in
            leaving_atom_name = leaving_atom_info_tuple[0]
            involved_bond_indices = list()
            for bond_idx, bond_info_tuple in enumerate(bonds_info):
                if leaving_atom_name in bond_info_tuple[:2]:
                    involved_bond_indices.append(bond_idx)
            # If this leaving atom is involved in more than one bond, then
            # don't remove it, since doing so would disconnect the molecule
            if len(involved_bond_indices) != 1:
                continue
            atoms_info_copy = copy.deepcopy(atoms_info)
            bonds_info_copy = copy.deepcopy(bonds_info)

            # Determine whether the neighbor of the leaving atom is also leaving
            bond_info_tuple = bonds_info[involved_bond_indices[0]]
            neighbor_atom_name = [
                name for name in bond_info_tuple[:2] if name != leaving_atom_name
            ][0]
            for neighbor_atom_idx, neighbor_atom_info_tuple in enumerate(atoms_info):
                if neighbor_atom_info_tuple[0] != neighbor_atom_name:
                    continue
                neighbor_atom_is_also_leaving = neighbor_atom_info_tuple[5] == "Y"
                # If the neighbor atom is ALSO leaving, then deduct 1 from its formal charge
                if neighbor_atom_is_also_leaving:
                    new_atom_tuple = (
                        neighbor_atom_info_tuple[0],
                        neighbor_atom_info_tuple[1],
                        neighbor_atom_info_tuple[2] - 1,
                        neighbor_atom_info_tuple[3],
                        neighbor_atom_info_tuple[4],
                        neighbor_atom_info_tuple[5],
                    )
                    atoms_info_copy[neighbor_atom_idx] = new_atom_tuple
            atoms_info_copy.pop(leaving_atom_idx)
            bonds_info_copy.pop(involved_bond_indices[0])
            return_list = self._recursive_prepare_atom_bond_info(
                atoms_info_copy, bonds_info_copy, return_list=return_list
            )

        return return_list

    def _patch_known_problems(self):
        """
        Monkey-patching known problems with current aa-variants-v1.cif from CCD.

        .. warning: Needed as of Oct-21-2021
        """

        substructures_to_fix = {
            "PRO": [
                (
                    "[N-:1]1[C@@:2]([C-:3]=[O:4])([H:8])[C:5]([H:9])([H:10])[C:6]([H:11])([H:12])[C:7]1([H:13])[H:14]",
                    "[N:1]1[C@@:2]([C:3]=[O:4])([H:8])[C:5]([H:9])([H:10])[C:6]([H:11])([H:12])[C:7]1([H:13])[H:14]",
                ),
            ],
            "HIS": [
                (
                    "[N:1]([C@:2]([C:3](=[O:4])[O:11][H:18])([C:5]([C:6]1=[C:8]([H:16])[N-:10][C:9]([H:17])=[N+:7]1[H:15])([H:13])[H:14])[H:12])([H:19])[H:20]",
                    "[N:1]([C@:2]([C:3](=[O:4])[O:11][H:18])([C:5]([C:6]1=[C:8]([H:16])[N:10]=[C:9]([H:17])[N:7]1[H:15])([H:13])[H:14])[H:12])([H:19])[H:20]",
                ),
                (
                    "[N:1]([C@:2]([C:3](=[O:4])[O-:11])([C:5]([C:6]1=[C:8]([H:16])[N-:10][C:9]([H:17])=[N+:7]1[H:15])([H:13])[H:14])[H:12])([H:18])[H:19]",
                    "[N:1]([C@:2]([C:3](=[O:4])[O-:11])([C:5]([C:6]1=[C:8]([H:16])[N:10]=[C:9]([H:17])[N:7]1[H:15])([H:13])[H:14])[H:12])([H:18])[H:19]",
                ),
                (
                    "[N:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:8]([H:15])[N-:10][C:9]([H:16])=[N+:7]1[H:14])([H:12])[H:13])[H:11])([H:17])[H:18]",
                    "[N:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:8]([H:15])[N:10]=[C:9]([H:16])[N:7]1[H:14])([H:12])[H:13])[H:11])([H:17])[H:18]",
                ),
                (
                    "[N:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:8]([H:15])[N-:10][C:9]([H:16])=[N+:7]1[H:14])([H:12])[H:13])[H:11])[H:17]",
                    "[N:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:8]([H:15])[N:10]=[C:9]([H:16])[N:7]1[H:14])([H:12])[H:13])[H:11])[H:17]",
                ),
                (
                    "[N:1]([C@:2]([C:3](=[O:4])[O-:11])([C:5]([C:6]1=[C:8]([H:16])[N-:10][C:9]([H:17])=[N+:7]1[H:15])([H:13])[H:14])[H:12])[H:18]",
                    "[N:1]([C@:2]([C:3](=[O:4])[O-:11])([C:5]([C:6]1=[C:8]([H:16])[N:10]=[C:9]([H:17])[N:7]1[H:15])([H:13])[H:14])[H:12])[H:18]",
                ),
                (
                    "[N:1]([C@:2]([C:3](=[O:4])[O:11][H:18])([C:5]([C:6]1=[C:8]([H:16])[N-:10][C:9]([H:17])=[N+:7]1[H:15])([H:13])[H:14])[H:12])[H:19]",
                    "[N:1]([C@:2]([C:3](=[O:4])[O:11][H:18])([C:5]([C:6]1=[C:8]([H:16])[N:10]=[C:9]([H:17])[N:7]1[H:15])([H:13])[H:14])[H:12])[H:19]",
                ),
                (
                    "[N+:1]([C@:2]([C:3](=[O:4])[O-:11])([C:5]([C:6]1=[C:8]([H:16])[N-:10][C:9]([H:17])=[N+:7]1[H:15])([H:13])[H:14])[H:12])([H:18])([H:19])[H:20]",
                    "[N+:1]([C@:2]([C:3](=[O:4])[O-:11])([C:5]([C:6]1=[C:8]([H:16])[N:10]=[C:9]([H:17])[N:7]1[H:15])([H:13])[H:14])[H:12])([H:18])([H:19])[H:20]",
                ),
                (
                    "[N+:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:8]([H:15])[N-:10][C:9]([H:16])=[N+:7]1[H:14])([H:12])[H:13])[H:11])([H:17])([H:18])[H:19]",
                    "[N+:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:8]([H:15])[N:10]=[C:9]([H:16])[N:7]1[H:14])([H:12])[H:13])[H:11])([H:17])([H:18])[H:19]",
                ),
            ],
            "TRP": [
                (
                    "[N:1]([C@:2]([C:3](=[O:4])[O:15][H:24])([C:5]([C:6]1=[C:7]([H:19])[N-:9][c:10]2[c:8]1[c:11]([H:20])[c:13]([H:22])[c:14]([H:23])[c:12]2[H:21])([H:17])[H:18])[H:16])([H:25])[H:26]",
                    "[N:1]([C@:2]([C:3](=[O:4])[O:15][H:24])([C:5]([C:6]1=[C:7]([H:19])[N:9][c:10]2[c:8]1[c:11]([H:20])[c:13]([H:22])[c:14]([H:23])[c:12]2[H:21])([H:17])[H:18])[H:16])([H:25])[H:26]",
                ),
                (
                    "[N:1]([C@:2]([C:3](=[O:4])[O-:15])([C:5]([C:6]1=[C:7]([H:19])[N-:9][c:10]2[c:8]1[c:11]([H:20])[c:13]([H:22])[c:14]([H:23])[c:12]2[H:21])([H:17])[H:18])[H:16])([H:24])[H:25]",
                    "[N:1]([C@:2]([C:3](=[O:4])[O-:15])([C:5]([C:6]1=[C:7]([H:19])[N:9][c:10]2[c:8]1[c:11]([H:20])[c:13]([H:22])[c:14]([H:23])[c:12]2[H:21])([H:17])[H:18])[H:16])([H:24])[H:25]",
                ),
                (
                    "[N:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:7]([H:18])[N-:9][c:10]2[c:8]1[c:11]([H:19])[c:13]([H:21])[c:14]([H:22])[c:12]2[H:20])([H:16])[H:17])[H:15])([H:23])[H:24]",
                    "[N:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:7]([H:18])[N:9][c:10]2[c:8]1[c:11]([H:19])[c:13]([H:21])[c:14]([H:22])[c:12]2[H:20])([H:16])[H:17])[H:15])([H:23])[H:24]",
                ),
                (
                    "[N+:1]([C@:2]([C:3](=[O:4])[O-:15])([C:5]([C:6]1=[C:7]([H:19])[N-:9][c:10]2[c:8]1[c:11]([H:20])[c:13]([H:22])[c:14]([H:23])[c:12]2[H:21])([H:17])[H:18])[H:16])([H:24])([H:25])[H:26]",
                    "[N+:1]([C@:2]([C:3](=[O:4])[O-:15])([C:5]([C:6]1=[C:7]([H:19])[N:9][c:10]2[c:8]1[c:11]([H:20])[c:13]([H:22])[c:14]([H:23])[c:12]2[H:21])([H:17])[H:18])[H:16])([H:24])([H:25])[H:26]",
                ),
                (
                    "[N+:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:7]([H:18])[N-:9][c:10]2[c:8]1[c:11]([H:19])[c:13]([H:21])[c:14]([H:22])[c:12]2[H:20])([H:16])[H:17])[H:15])([H:23])([H:24])[H:25]",
                    "[N+:1]([C@:2]([C:3]=[O:4])([C:5]([C:6]1=[C:7]([H:18])[N:9][c:10]2[c:8]1[c:11]([H:19])[c:13]([H:21])[c:14]([H:22])[c:12]2[H:20])([H:16])[H:17])[H:15])([H:23])([H:24])[H:25]",
                ),
            ],
        }
        # Fix PRO smarts substructure
        for aa_name, replacement_list in substructures_to_fix.items():
            for (old_smarts, new_smarts) in replacement_list:
                self.data[aa_name][new_smarts] = self.data[aa_name][old_smarts]
                self.data[aa_name].pop(old_smarts)

    def _add_common_substructures(self):
        """
        .. warning: Needed as of Oct-21-2021
        """

        # Add common caps
        self.data["ACE"] = {
            "[C:1](=[O:2])[C:3]([H:4])([H:5])[H:6]": ["C", "O", "CH3", "H1", "H2", "H3"]
        }
        self.data["NME"] = {
            "[N:1]([C:2]([H:4])([H:5])[H:6])[H:3]": ["N", "C", "H", "H1", "H2", "H3"]
        }
        self.data["NH2"] = {"[N:1]([H:2])[H:3]": ["N", "HN1", "HN2"]}

    def _add_common_linkages(self):
        """
        .. warning: Needed as of Oct-21-2021
        """

        self.data["PEPTIDE_BOND"] = {
            "[C:1](=[O:2])[N:3]([C:4])": [
                "C",
                "O",
                "N",
                "CA",
            ]
        }
        self.data["DISULFIDE"] = {
            "[S:1][S:2]": [
                "SG",
                "SG",
            ]
        }
