"""
Functions which create a topological molecule directly, without a toolkit.

These are common to several test modules.

"""

import numpy as np
from openff.units import unit

from openff.toolkit.tests.utils import requires_openeye
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.utils import get_data_file_path


def create_cis_1_2_dichloroethene():
    """
    Creates an openff.toolkit.topology.Molecule representation of cis-1,2-dichloroethene
    without the use of a cheminformatics toolkit.
    """

    cis_dichloroethene = Molecule()
    cis_dichloroethene.add_atom(17, 0, False)
    cis_dichloroethene.add_atom(6, 0, False)
    cis_dichloroethene.add_atom(6, 0, False)
    cis_dichloroethene.add_atom(17, 0, False)
    cis_dichloroethene.add_atom(1, 0, False)
    cis_dichloroethene.add_atom(1, 0, False)
    cis_dichloroethene.add_bond(0, 1, 1, False)
    cis_dichloroethene.add_bond(1, 2, 2, False, "Z")
    cis_dichloroethene.add_bond(2, 3, 1, False)
    cis_dichloroethene.add_bond(1, 4, 1, False)
    cis_dichloroethene.add_bond(2, 5, 1, False)
    return cis_dichloroethene


def create_ethanol():
    """
    Creates an openff.toolkit.topology.Molecule representation of
    ethanol without the use of a cheminformatics toolkit
    """
    # Create an ethanol molecule without using a toolkit
    ethanol = Molecule()
    ethanol.add_atom(6, 0, False)  # C0
    ethanol.add_atom(6, 0, False)  # C1
    ethanol.add_atom(8, 0, False)  # O2
    ethanol.add_atom(1, 0, False)  # H3
    ethanol.add_atom(1, 0, False)  # H4
    ethanol.add_atom(1, 0, False)  # H5
    ethanol.add_atom(1, 0, False)  # H6
    ethanol.add_atom(1, 0, False)  # H7
    ethanol.add_atom(1, 0, False)  # H8
    ethanol.add_bond(0, 1, 1, False, fractional_bond_order=1.33)  # C0 - C1
    ethanol.add_bond(1, 2, 1, False, fractional_bond_order=1.23)  # C1 - O2
    ethanol.add_bond(0, 3, 1, False, fractional_bond_order=1)  # C0 - H3
    ethanol.add_bond(0, 4, 1, False, fractional_bond_order=1)  # C0 - H4
    ethanol.add_bond(0, 5, 1, False, fractional_bond_order=1)  # C0 - H5
    ethanol.add_bond(1, 6, 1, False, fractional_bond_order=1)  # C1 - H6
    ethanol.add_bond(1, 7, 1, False, fractional_bond_order=1)  # C1 - H7
    ethanol.add_bond(2, 8, 1, False, fractional_bond_order=1)  # O2 - H8
    charges = unit.Quantity(
        np.array([-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4]),
        unit.elementary_charge,
    )
    ethanol.partial_charges = charges

    return ethanol


def create_reversed_ethanol():
    """
    Creates an openff.toolkit.topology.Molecule representation of
    ethanol without the use of a cheminformatics toolkit. This function
    reverses the atom indexing of create_ethanol
    """
    # Create an ethanol molecule without using a toolkit
    ethanol = Molecule()
    ethanol.add_atom(1, 0, False)  # H0
    ethanol.add_atom(1, 0, False)  # H1
    ethanol.add_atom(1, 0, False)  # H2
    ethanol.add_atom(1, 0, False)  # H3
    ethanol.add_atom(1, 0, False)  # H4
    ethanol.add_atom(1, 0, False)  # H5
    ethanol.add_atom(8, 0, False)  # O6
    ethanol.add_atom(6, 0, False)  # C7
    ethanol.add_atom(6, 0, False)  # C8
    ethanol.add_bond(8, 7, 1, False, fractional_bond_order=1.33)  # C8 - C7
    ethanol.add_bond(7, 6, 1, False, fractional_bond_order=1.23)  # C7 - O6
    ethanol.add_bond(8, 5, 1, False, fractional_bond_order=1)  # C8 - H5
    ethanol.add_bond(8, 4, 1, False, fractional_bond_order=1)  # C8 - H4
    ethanol.add_bond(8, 3, 1, False, fractional_bond_order=1)  # C8 - H3
    ethanol.add_bond(7, 2, 1, False, fractional_bond_order=1)  # C7 - H2
    ethanol.add_bond(7, 1, 1, False, fractional_bond_order=1)  # C7 - H1
    ethanol.add_bond(6, 0, 1, False, fractional_bond_order=1)  # O6 - H0
    charges = unit.Quantity(
        np.array([0.4, 0.3, 0.2, 0.1, 0.00001, -0.1, -0.2, -0.3, -0.4]),
        unit.elementary_charge,
    )
    ethanol.partial_charges = charges
    return ethanol


def create_benzene_no_aromatic():
    """
    Creates an openff.toolkit.topology.Molecule representation of benzene through the API with aromatic bonds
    not defied, used to test the levels of isomorphic matching.
    """
    benzene = Molecule()
    benzene.add_atom(6, 0, False)  # C0
    benzene.add_atom(6, 0, False)  # C1
    benzene.add_atom(6, 0, False)  # C2
    benzene.add_atom(6, 0, False)  # C3
    benzene.add_atom(6, 0, False)  # C4
    benzene.add_atom(6, 0, False)  # C5
    benzene.add_atom(1, 0, False)  # H6
    benzene.add_atom(1, 0, False)  # H7
    benzene.add_atom(1, 0, False)  # H8
    benzene.add_atom(1, 0, False)  # H9
    benzene.add_atom(1, 0, False)  # H10
    benzene.add_atom(1, 0, False)  # H11
    benzene.add_bond(0, 5, 1, False)  # C0 - C5
    benzene.add_bond(0, 1, 1, False)  # C0 - C1
    benzene.add_bond(1, 2, 1, False)  # C1 - C2
    benzene.add_bond(2, 3, 1, False)  # C2 - C3
    benzene.add_bond(3, 4, 1, False)  # C3 - C4
    benzene.add_bond(4, 5, 1, False)  # C4 - C5
    benzene.add_bond(0, 6, 1, False)  # C0 - H6
    benzene.add_bond(1, 7, 1, False)  # C1 - C7
    benzene.add_bond(2, 8, 1, False)  # C2 - C8
    benzene.add_bond(3, 9, 1, False)  # C3 - C9
    benzene.add_bond(4, 10, 1, False)  # C4 - C10
    benzene.add_bond(5, 11, 1, False)  # C5 - C11
    return benzene


def create_acetaldehyde():
    """
    Creates an openff.toolkit.topology.Molecule representation of acetaldehyde through the API
    """
    acetaldehyde = Molecule()
    acetaldehyde.add_atom(6, 0, False)  # C0
    acetaldehyde.add_atom(6, 0, False)  # C1
    acetaldehyde.add_atom(8, 0, False)  # O2
    acetaldehyde.add_atom(1, 0, False)  # H3
    acetaldehyde.add_atom(1, 0, False)  # H4
    acetaldehyde.add_atom(1, 0, False)  # H5
    acetaldehyde.add_atom(1, 0, False)  # H6
    acetaldehyde.add_bond(0, 1, 1, False)  # C0 - C1
    acetaldehyde.add_bond(1, 2, 2, False)  # C1 = O2
    acetaldehyde.add_bond(0, 3, 1, False)  # C0 - H3
    acetaldehyde.add_bond(0, 4, 1, False)  # C0 - H4
    acetaldehyde.add_bond(0, 5, 1, False)  # C0 - H5
    acetaldehyde.add_bond(1, 6, 1, False)  # C1 - H6
    charges = unit.Quantity(
        np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), unit.elementary_charge
    )
    acetaldehyde.partial_charges = charges
    return acetaldehyde


def create_water():
    """
    Creates an openff.toolkit.topology.Molecule representation of water through the API
    """
    mol = Molecule()
    mol.add_atom(1, 0, False)  # H1
    mol.add_atom(8, 0, False)  # O
    mol.add_atom(1, 0, False)  # H2
    mol.add_bond(0, 1, 1, False)  # H1 - O
    mol.add_bond(1, 2, 1, False)  # O - H2
    charges = unit.Quantity(np.array([0.0, 0.0, 0.0]), unit.elementary_charge)
    mol.partial_charges = charges
    return mol


def create_ammonia():
    """
    Creates an openff.toolkit.topology.Molecule representation of ammonia through the API
    """
    mol = Molecule()
    mol.add_atom(1, 0, False)  # H1
    mol.add_atom(7, 0, False)  # N
    mol.add_atom(1, 0, False)  # H2
    mol.add_atom(1, 0, False)  # H3
    mol.add_bond(0, 1, 1, False)  # H1 - N
    mol.add_bond(1, 2, 1, False)  # N - H2
    mol.add_bond(1, 3, 1, False)  # N - H3
    charges = unit.Quantity(np.array([0.0, 0.0, 0.0, 0.0]), unit.elementary_charge)
    mol.partial_charges = charges
    return mol


def create_cyclic_n3h3():
    """
    Creates an openff.toolkit.topology.Molecule representation of N3H3 through the API
    """
    mol = Molecule()
    mol.add_atom(7, 0, False)  # N1
    mol.add_atom(7, 0, False)  # N2
    mol.add_atom(7, 0, False)  # N3
    mol.add_atom(1, 0, False)  # H1
    mol.add_atom(1, 0, False)  # H2
    mol.add_atom(1, 0, False)  # H3
    mol.add_bond(0, 1, 1, False)  # N1 - N2
    mol.add_bond(1, 2, 1, False)  # N2 - N3
    mol.add_bond(2, 0, 1, False)  # N1 - N3
    mol.add_bond(0, 3, 1, False)  # N1 - H1
    mol.add_bond(1, 4, 1, False)  # N2 - H2
    mol.add_bond(2, 5, 1, False)  # N3 - H3
    return mol


def create_acetate():
    """
    Creates an openff.toolkit.topology.Molecule representation of
    acetate without the use of a cheminformatics toolkit
    """
    # Create an acetate molecule without using a toolkit
    acetate = Molecule()
    acetate.add_atom(6, 0, False)  # C0
    acetate.add_atom(6, 0, False)  # C1
    acetate.add_atom(8, 0, False)  # O2
    acetate.add_atom(8, -1, False)  # O3
    acetate.add_atom(1, 0, False)  # H4
    acetate.add_atom(1, 0, False)  # H5
    acetate.add_atom(1, 0, False)  # H6
    acetate.add_bond(0, 1, 1, False)  # C0 - C1
    acetate.add_bond(1, 2, 2, False)  # C1 = O2
    acetate.add_bond(1, 3, 1, False)  # C1 - O3[-1]
    acetate.add_bond(0, 4, 1, False)  # C0 - H4
    acetate.add_bond(0, 5, 1, False)  # C0 - H5
    acetate.add_bond(0, 6, 1, False)  # C0 - H6
    return acetate


def create_cyclohexane():
    """
    Creates an openff.toolkit.topology.Molecule representation of
    cyclohexane without the use of a cheminformatics toolkit
    """
    cyclohexane = Molecule()
    cyclohexane.add_atom(6, 0, False)  # C0
    cyclohexane.add_atom(6, 0, False)  # C1
    cyclohexane.add_atom(6, 0, False)  # C2
    cyclohexane.add_atom(6, 0, False)  # C3
    cyclohexane.add_atom(6, 0, False)  # C4
    cyclohexane.add_atom(6, 0, False)  # C5
    cyclohexane.add_atom(1, 0, False)  # H6
    cyclohexane.add_atom(1, 0, False)  # H7
    cyclohexane.add_atom(1, 0, False)  # H8
    cyclohexane.add_atom(1, 0, False)  # H9
    cyclohexane.add_atom(1, 0, False)  # H10
    cyclohexane.add_atom(1, 0, False)  # H11
    cyclohexane.add_atom(1, 0, False)  # H12
    cyclohexane.add_atom(1, 0, False)  # H13
    cyclohexane.add_atom(1, 0, False)  # H14
    cyclohexane.add_atom(1, 0, False)  # H15
    cyclohexane.add_atom(1, 0, False)  # H16
    cyclohexane.add_atom(1, 0, False)  # H17
    cyclohexane.add_bond(0, 1, 1, False)  # C0 - C1
    cyclohexane.add_bond(1, 2, 1, False)  # C1 - C2
    cyclohexane.add_bond(2, 3, 1, False)  # C2 - C3
    cyclohexane.add_bond(3, 4, 1, False)  # C3 - C4
    cyclohexane.add_bond(4, 5, 1, False)  # C4 - C5
    cyclohexane.add_bond(5, 0, 1, False)  # C5 - C0
    cyclohexane.add_bond(0, 6, 1, False)  # C0 - H6
    cyclohexane.add_bond(0, 7, 1, False)  # C0 - H7
    cyclohexane.add_bond(1, 8, 1, False)  # C1 - H8
    cyclohexane.add_bond(1, 9, 1, False)  # C1 - H9
    cyclohexane.add_bond(2, 10, 1, False)  # C2 - H10
    cyclohexane.add_bond(2, 11, 1, False)  # C2 - H11
    cyclohexane.add_bond(3, 12, 1, False)  # C3 - H12
    cyclohexane.add_bond(3, 13, 1, False)  # C3 - H13
    cyclohexane.add_bond(4, 14, 1, False)  # C4 - H14
    cyclohexane.add_bond(4, 15, 1, False)  # C4 - H15
    cyclohexane.add_bond(5, 16, 1, False)  # C5 - H16
    cyclohexane.add_bond(5, 17, 1, False)  # C5 - H17
    return cyclohexane


def create_dioxygen():
    """
    Creates an openff.toolkit.topology.Molecule representation of
    dioxygen without the use of a cheminformatics toolkit
    """
    dioxygen = Molecule()
    dioxygen.add_atom(8, 0, False)  # O0
    dioxygen.add_atom(8, 0, False)  # O1
    dioxygen.add_bond(0, 1, 2, False)  # O0 # O1
    charges = unit.Quantity(np.array([0.0, 0.0]), unit.elementary_charge)
    dioxygen.partial_charges = charges

    return dioxygen


def create_dinitrogen():
    """
    Creates an openff.toolkit.topology.Molecule representation of
    dinitrogen without the use of a cheminformatics toolkit
    """
    dinitrogen = Molecule()
    dinitrogen.add_atom(7, 0, False)  # N0
    dinitrogen.add_atom(7, 0, False)  # N1
    dinitrogen.add_bond(0, 1, 3, False)  # N0 - N1
    charges = unit.Quantity(np.array([0.0, 0.0]), unit.elementary_charge)
    dinitrogen.partial_charges = charges
    return dinitrogen


def dipeptide():
    dipeptide = Molecule.from_file(get_data_file_path("proteins/CTerminal_ALA.sdf"))
    return dipeptide


def dipeptide_residues_perceived():
    dipeptide_residues_perceived = Molecule(dipeptide())
    dipeptide_residues_perceived.perceive_residues()
    return dipeptide_residues_perceived


def dipeptide_hierarchy_added():
    dipeptide_hierarchy_perceived = Molecule(dipeptide_residues_perceived())
    dipeptide_hierarchy_perceived.add_default_hierarchy_schemes()
    return dipeptide_hierarchy_perceived


def cyx():
    cyx = Molecule.from_file(get_data_file_path("proteins/MainChain_CYX.sdf"))
    return cyx


def cyx_residues_perceived():
    cyx_residues_perceived = Molecule(cyx())
    cyx_residues_perceived.perceive_residues()
    return cyx_residues_perceived


def cyx_hierarchy_added():
    cyx_hierarchy_perceived = Molecule(cyx_residues_perceived())
    cyx_hierarchy_perceived.add_default_hierarchy_schemes()
    return cyx_hierarchy_perceived


def empty_molecule():
    return Molecule()


def ethane_from_smiles():
    return Molecule.from_smiles("CC")


def ethene_from_smiles():
    return Molecule.from_smiles("C=C")


def propane_from_smiles():
    return Molecule.from_smiles("CCC")


def toluene_from_sdf():
    filename = get_data_file_path("molecules/toluene.sdf")
    return Molecule.from_file(filename)


@requires_openeye
def toluene_from_charged_mol2():
    filename = get_data_file_path("molecules/toluene_charged.mol2")
    # TODO: This will require openeye to load
    return Molecule.from_file(filename)


def charged_methylamine_from_smiles():
    return Molecule.from_smiles("[H]C([H])([H])[N+]([H])([H])[H]")


def ethane_from_smiles_w_vsites():
    molecule = Molecule.from_smiles("CC")
    carbons = [atom for atom in molecule.atoms if atom.atomic_number == 6]
    c0_hydrogens = [atom for atom in carbons[0].bonded_atoms if atom.atomic_number == 1]
    molecule.add_bond_charge_virtual_site(
        (carbons[0], carbons[1]),
        0.1 * unit.angstrom,
        charge_increments=[0.1, 0.05] * unit.elementary_charge,
    )
    molecule.add_monovalent_lone_pair_virtual_site(
        (c0_hydrogens[0], carbons[0], carbons[1]),
        0.2 * unit.angstrom,
        20 * unit.degree,
        25 * unit.degree,
        charge_increments=[0.01, 0.02, 0.03] * unit.elementary_charge,
    )
    return molecule


def propane_from_smiles_w_vsites():
    # Make a propane with virtual sites
    molecule = Molecule.from_smiles("CCC")
    carbons = [atom for atom in molecule.atoms if atom.atomic_number == 6]
    c0_hydrogens = [atom for atom in carbons[0].bonded_atoms if atom.atomic_number == 1]
    # This will add *two* particles (symmetric=True), *one* virtual site
    molecule.add_monovalent_lone_pair_virtual_site(
        (c0_hydrogens[0], carbons[0], carbons[1]),
        0.2 * unit.angstrom,
        20 * unit.degree,
        25 * unit.degree,
        charge_increments=[0.01, 0.02, 0.03] * unit.elementary_charge,
        symmetric=True,
    )
    # This will add *one* particle (symmetric=False), *one* virtual site
    molecule.add_bond_charge_virtual_site(
        (carbons[0], carbons[1]),
        0.1 * unit.angstrom,
        charge_increments=[0.1, 0.05] * unit.elementary_charge,
        symmetric=False,
    )
    return molecule


def tip5_water():
    # Make a TIP5 water
    molecule = Molecule.from_smiles("[H][O][H]")
    O1 = [atom for atom in molecule.atoms if atom.atomic_number == 8][0]
    H1, H2 = [atom for atom in O1.bonded_atoms if atom.atomic_number == 1]
    molecule.add_divalent_lone_pair_virtual_site(
        (H1, O1, H2),
        0.7 * unit.angstrom,
        54.71384225 * unit.degree,
        charge_increments=[0.1205, 0.00, 0.1205] * unit.elementary_charge,
        symmetric=True,
    )
    return molecule


def topology_with_metadata():
    n2 = Molecule.from_smiles("N#N")
    ammonia = create_ammonia()
    chloride = Molecule.from_smiles("[Cl-]")
    water = create_water()
    mols = [
        n2,
        ammonia,
        chloride,
        chloride,
        ammonia,
        chloride,
        chloride,
        n2,
        chloride,
        chloride,
        chloride,
        water,
    ]

    top = Topology.from_molecules(mols)

    atom_hier_info = [
        ("AAA", 1, " ", "A"),  # mol[0]
        ("AAA", 1, " ", "A"),  # mol[0]
        ("AAA", 1, " ", "A"),  # mol[1]
        ("BBB", 1, " ", "A"),  # mol[1]
        ("BBB", 2, " ", "A"),  # mol[1]
        ("BBB", 2, " ", "B"),  # mol[1]
        ("CCC", 2, " ", "B"),  # mol[2]
        ("CCC", 3, " ", "B"),  # mol[3]
        ("CCC", 3, " ", "C"),  # mol[4]
        ("DDD", 4, " ", "C"),  # mol[4]
        ("EEE", 4, " ", "D"),  # mol[4]
        ("EEE", 5, " ", "E"),  # mol[4]
        ("FFF", 6, " ", "E"),  # mol[5]
        ("GGG", 6, " ", "F"),  # mol[6]
        ("GGG", 7, " ", "G"),  # mol[7]
        ("HHH", 8, " ", "H"),  # mol[7]
        ("YZ", 9, " ", "I"),  # mol[8]
        ("AAA", 1, " ", "A"),  # mol[9]
        ("AAA", 1, "X", "A"),  # mol[10]
        (None, None, " ", None),  # mol[11]
        (None, None, " ", None),  # mol[12]
        ("AAA", None, " ", None),  # mol[13]
        (None, None, " ", "A"),  # mol[14]
    ]

    for atom, metadata_tuple in zip(top.atoms, atom_hier_info):
        residue_name = metadata_tuple[0]
        residue_number = metadata_tuple[1]
        insertion_code = metadata_tuple[2]
        chain_id = metadata_tuple[3]
        if residue_name is not None:
            atom.metadata["residue_name"] = residue_name
        if residue_number is not None:
            atom.metadata["residue_number"] = residue_number
        if insertion_code is not None:
            atom.metadata["insertion_code"] = insertion_code
        if chain_id is not None:
            atom.metadata["chain_id"] = chain_id

    return top
