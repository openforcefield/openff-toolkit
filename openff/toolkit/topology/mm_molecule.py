"""
TypedMolecule TODOs
* Testing
* Detemine whether atomic numbers are required, or whether a TypedMol could have a mix of TypedAtom and TypedParticles, where the latter doesn't need atomic nubmers
* Determine whether typedmols should have isomorphism operations
* Should be able to go to hill formula (?), networkx, openmm
* Topology serialization will have trouble here - Won't know whether it's trying to deserialize a Molecule or a TypedMolecule.

"""


class _TypedMolecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.hierarchy_schemes = dict()
        self.conformers = None

    def add_atom(self, **kwargs):
        atom = _TypedAtom(self, **kwargs)
        self.atoms.append(atom)

    def add_bond(self, atom1, atom2, **kwargs):
        if isinstance(atom1, int) and isinstance(atom2, int):
            atom1_atom = self.atoms[atom1]
            atom2_atom = self.atoms[atom2]
        elif isinstance(atom1, Atom) and isinstance(atom2, Atom):
            atom1_atom = atom1
            atom2_atom = atom2
        else:
            raise Exception(
                "Invalid inputs to molecule._add_bond. Expected ints or Atoms. "
                "Received {} (type {}) and {} (type {}) ".format(
                    atom1, type(atom1), atom2, type(atom2)
                )
            )
        bond = _TypedBond(atom1_atom, atom2_atom, **kwargs)
        self.bonds.append(bond)

    def add_conformer(self, conformer):
        if self.conformers is None:
            self.conformers = list()
        self.conformers.append(conformer)

    @property
    def n_atoms(self):
        return len(self.atoms)

    @property
    def n_particles(self):
        return len(self.atoms)

    @property
    def n_bonds(self):
        return len(self.bonds)

    @property
    def n_conformers(self):
        if self.conformers is None:
            return 0
        return len(self.conformers)

    def atom(self, index):
        return self.atoms[index]

    def atom_index(self, atom):
        return self.atoms.index(atom)

    def bond(self, index):
        return self.bonds[index]

    def get_bond_between(self, atom1_index, atom2_index):
        for bond in self.atom(atom1_index).bonds:
            for atom in bond.atoms:
                if atom.molecule_atom_index == atom2_index:
                    return bond

    def to_dict(self):
        molecule_dict = dict()
        special_serialization_logic = [
            "atoms",
            "bonds",
            "conformers",
            "hierarchy_schemes",
            "molecule",
        ]
        for attr_name, attr_val in self.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in special_serialization_logic:
                continue
            molecule_dict[attr_name] = attr_val

        atom_list = list()
        for atom in self.atoms:
            atom_list.append(atom.to_dict())
        molecule_dict["atoms"] = atom_list

        bond_list = list()
        for bond in self.bonds:
            bond_list.append(bond.to_dict())
        molecule_dict["bonds"] = bond_list

        if self.conformers is None:
            molecule_dict["conformers"] = None
        else:
            molecule_dict["conformers"] = []
            molecule_dict[
                "conformers_unit"
            ] = "angstrom"  # Have this defined as a class variable?
            for conf in self.conformers:
                conf_unitless = conf.m_as(unit.angstrom)
                conf_serialized, conf_shape = serialize_numpy((conf_unitless))
                molecule_dict["conformers"].append(conf_serialized)

        molecule_dict["hierarchy_schemes"] = dict()
        for iter_name, hier_scheme in self.hierarchy_schemes.items():
            molecule_dict["hierarchy_schemes"][iter_name] = hier_scheme.to_dict()

        return molecule_dict

    @classmethod
    def from_dict(cls, molecule_dict):
        molecule = cls()
        molecule.initialize_from_dict(molecule_dict)
        return molecule

    def initialize_from_dict(self, molecule_dict):
        atom_dicts = molecule_dict.pop("atoms")
        for atom_dict in atom_dicts:
            self.add_atom(**atom_dict)

        bond_dicts = molecule_dict.pop("bonds")
        print(bond_dicts)
        for bond_dict in bond_dicts:
            bond_dict["atom1"] = int(bond_dict["atom1"])
            bond_dict["atom2"] = int(bond_dict["atom2"])
            self.add_bond(**bond_dict)

        conformers = molecule_dict.pop("conformers")
        if conformers is None:
            self.conformers = None
        else:
            self.conformers = list()
            for ser_conf in molecule_dict["conformers"]:
                # TODO: Update to use string_to_quantity
                conformers_shape = (self.n_atoms, 3)
                conformer_unitless = deserialize_numpy(ser_conf, conformers_shape)
                c_unit = getattr(unit, molecule_dict["conformers_unit"])
                conformer = unit.Quantity(conformer_unitless, c_unit)
                self.conformers.append(conformer)

        hier_scheme_dicts = molecule_dict.pop("hierarchy_schemes")
        for iter_name, hierarchy_scheme_dict in hier_scheme_dicts.items():
            new_hier_scheme = self.add_hierarchy_scheme(
                hierarchy_scheme_dict["uniqueness_criteria"],
                iter_name,
            )
            for element_dict in hierarchy_scheme_dict["hierarchy_elements"]:
                new_hier_scheme.add_hierarchy_element(
                    element_dict["identifier"], element_dict["particle_indices"]
                )
            self._expose_hierarchy_scheme(iter_name)

        for key, val in molecule_dict:
            setattr(self, key, val)

    def particle(self, index):
        return self.atom(index)

    @property
    def particles(self):
        return self.atoms

    def particle_index(self, particle):
        return self.atom_index(particle)

    def to_hill_formula(self) -> str:
        """
        Generate the Hill formula of this molecule.

        Returns
        ----------
        formula : the Hill formula of the molecule

        Raises
        -----------
        NotImplementedError : if the molecule is not of one of the specified types.
        """
        atom_nums = [atom.atomic_number for atom in self.atoms]

        return _atom_nums_to_hill_formula(atom_nums)


class _TypedAtom:
    def __init__(self, molecule, metadata=None, **kwargs):
        if metadata is None:
            self.metadata = AtomMetadataDict()
        else:
            self.metadata = AtomMetadataDict(metadata)
        self.molecule = molecule
        self.bonds = []
        for key, val in kwargs.items():
            setattr(self, key, val)

    def add_bond(self, bond):
        self.bonds.append(bond)

    @property
    def bonded_atoms(self):
        bonded_atoms = []
        for bond in self.bonds:
            for atom in bond:
                if atom is not self:
                    bonded_atoms.append(atom)
        return bonded_atoms

    @property
    def molecule_atom_index(self):
        return self.molecule.atoms.index(self)

    @property
    def molecule_particle_index(self):
        return self.molecule.atoms.index(self)

    def to_dict(self):
        atom_dict = dict()
        atom_dict["metadata"] = dict(self.metadata)
        special_serialization_logic = ["metadata", "molecule", "bonds"]

        for attr_name, attr_val in self.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in special_serialization_logic:
                continue
            atom_dict[attr_name] = attr_val
        return atom_dict


class _TypedBond:
    def __init__(self, atom1, atom2, **kwargs):
        self.molecule = atom1.molecule
        self.atom1 = atom1
        self.atom2 = atom2

        atom1.add_bond(self)
        atom2.add_bond(self)

        for key, val in kwargs.items():
            setattr(self, key, val)

    @property
    def atoms(self):
        return [self.atom1, self.atom2]

    @property
    def atom1_index(self):
        return self.atom1.molecule_atom_index

    @property
    def atom2_index(self):
        return self.atom2.molecule_atom_index

    def to_dict(self):
        bond_dict = dict()
        bond_dict["atom1"] = self.atom1.molecule_atom_index
        bond_dict["atom2"] = self.atom2.molecule_atom_index
        special_serialization_logic = ["atom1", "atom2", "molecule"]

        for attr_name, attr_val in self.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in special_serialization_logic:
                continue
            bond_dict[attr_name] = attr_val
        return bond_dict
