import os

from openff.toolkit.utils import get_data_file_path
from _cif_to_substructure_dict import CifSubstructures

if not os.path.exists("aa-variants-v1.cif"):
    import requests

    r = requests.get("https://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif")
    print(r.ok)
    if r.ok:
        with open("aa-variants-v1.cif", "wb") as of:
            for chunk in r.iter_content(chunk_size=1024 * 8):
                if chunk:
                    of.write(chunk)
                    of.flush()
                    os.fsync(of.fileno())


cif_object = CifSubstructures()
cif_object.from_file(
    "aa-variants-v1.cif",
    replace_quadruple_bond_with_any=False,
    remove_charge_bond_order_resonant=False,
)

# Automatically patch known problems - better that this explodes when things are fixed
cif_object._patch_known_problems()
cif_object._add_common_substructures()
cif_object._add_common_linkages()

output_file = get_data_file_path("proteins/T4-protein.sdf").replace(
    "T4-protein.sdf", "aa_residues_substructures_explicit_bond_orders_with_caps.json"
)
cif_object.to_json_file(output_file)

# TODO: reformat and clean this
from rdkit import Chem
import json
from openff.toolkit.utils.utils import get_data_file_path
from collections import OrderedDict
from copy import deepcopy
import requests

def add_atom(rdmol, begin_atom_idx, new_atom_smarts):
    # adds atom with single bond, returns new mol
    emol = Chem.RWMol(deepcopy(rdmol))
    new_atom = Chem.AtomFromSmarts(new_atom_smarts)
    idx = emol.AddAtom(new_atom)
    assert idx == rdmol.GetNumAtoms()
    emol.AddBond(begin_atom_idx, idx)
    return emol.GetMol()

def _get_maximum_map_num(rdmol):
    return max([atom.GetAtomMapNum() for atom in rdmol.GetAtoms()])

def _fill_out_query(rdmol):
    current_map_num = 1 + _get_maximum_map_num(rdmol) # do not change existing atom map numbers! 
    for atom in rdmol.GetAtoms():
        if atom.GetAtomicNum() > 0:
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(current_map_num)
                current_map_num += 1
            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
            a_num = atom.GetAtomicNum()
            D_num = len([0 for _ in atom.GetBonds()])
            F_num = atom.GetFormalCharge()
            query_string = f"[#{a_num}D{D_num}{F_num:+}:{current_map_num}]"
            query = Chem.AtomFromSmarts(query_string)
            atom.SetQuery(query)
        elif atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(current_map_num)
            current_map_num += 1
            query_string = f"[*:{current_map_num}]"
            query = Chem.AtomFromSmarts(query_string)
            atom.SetQuery(query)
        else:
            raise Exception
    for bond in rdmol.GetBonds():
        if bond.GetIsAromatic():
            query = Chem.BondFromSmarts(":")
            bond.SetQuery(query)
        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            query = Chem.BondFromSmarts("-")
            bond.SetQuery(query)
        
        
    return rdmol

def get_atom_with_map_num(rdmol, map_num):
    for atom in rdmol.GetAtoms():
        if atom.GetAtomMapNum() == map_num:
            return atom
    return None

def clear_query(rdmol):
    for atom in rdmol.GetAtoms():
        atom.SetAtomMapNum(0)
        n = atom.GetAtomicNum()
        c = atom.GetFormalCharge()
        query_string = f"[#{n}{c:+}]"
        query = Chem.AtomFromSmarts(query_string)
        atom.SetQuery(query)

substructure_file_path = get_data_file_path(
    "proteins/aa_residues_substructures_explicit_bond_orders_with_caps.json"
)

with open(substructure_file_path, "r") as subfile:
    substructure_dictionary = json.load(
        subfile, object_pairs_hook=OrderedDict
    )

new_subs_dict = OrderedDict()
amino_acid_subs = [(Chem.MolFromSmarts("[N+:1](-[H])(-[H])1[C@@:2]([C:3](=[O:4]))([H:9])[C:5]([H:10])([H:11])[C:6]([H:12])([H:13])[C:7]1([H:14])[H:15]"), [], [], 4),
                   (Chem.MolFromSmarts("[N+0:1](-[H])1[C@@:2]([C:3](=[O:4]))([H:9])[C:5]([H:10])([H:11])[C:6]([H:12])([H:13])[C:7]1([H:14])[H:15]"), [], [], 3),
                   (Chem.MolFromSmarts("[N+1:1](-[H])1[C@@:2]([C:3](=[O:4]))([H:9])[C:5]([H:10])([H:11])[C:6]([H:12])([H:13])[C:7]1([H:14])[H:15]"), [0], [], 3),
                   (Chem.MolFromSmarts("[N+0:1]1[C@@:2]([C:3](=[O:4]))([H:9])[C:5]([H:10])([H:11])[C:6]([H:12])([H:13])[C:7]1([H:14])[H:15]"), [0], [], 2),
                   (Chem.MolFromSmarts("[N+](-[H])(-[H])(-[H])-[C@:2]([C:3](=[O:4]))([C:5]([S+0:6]([H]))([H:9])[H:10])[H:8]"), [], [], 5),
                   (Chem.MolFromSmarts("[N+](-[H])(-[H])(-[H])-[C@:2]([C:3](=[O:4]))([C:5]([S+0:6])([H:9])[H:10])[H:8]"), [8], [], 5),
                   (Chem.MolFromSmarts("[N](-[H])(-[H])-[C@:2]([C:3](=[O:4]))([C:5]([S+0:6]([H]))([H:9])[H:10])[H:8]"), [], [], 4),
                   (Chem.MolFromSmarts("[N](-[H])(-[H])-[C@:2]([C:3](=[O:4]))([C:5]([S+0:6])([H:9])[H:10])[H:8]"), [7], [], 4),
                   (Chem.MolFromSmarts("[N](-[H])-[C@:2]([C:3](=[O:4]))([C:5]([S+0:6]([H]))([H:9])[H:10])[H:8]"), [0], [], 3),
                   (Chem.MolFromSmarts("[N](-[H])-[C@:2]([C:3](=[O:4]))([C:5]([S+0:6])([H:9])[H:10])[H:8]"), [0, 6], [], 3),
                   (Chem.MolFromSmarts("[N:1][C@:2]([C:3](=[O:4]))([C:5]([S+0:6]([H]))([H:9])[H:10])[H:8]"), [0,0], [], 2),
                   (Chem.MolFromSmarts("[N:1][C@:2]([C:3](=[O:4]))([C:5]([S+0:6])([H:9])[H:10])[H:8]"), [0,0, 5], [], 2),
                   (Chem.MolFromSmarts("[N+](-[H])(-[H])(-[H])-[C@](-[H])(-[*])-[C](=[O])"), [], [], 7),
                   (Chem.MolFromSmarts("[N](-[H])(-[H])-[C@](-[H])(-[*])-[C](=[O])"), [], [], 6),
                   (Chem.MolFromSmarts("[N](-[H])-[C@](-[H])(-[*])-[C](=[O])"), [0], [], 5),
                   ]
special_cases = {'[C:1](=[O:2])[C:3]([H:4])([H:5])[H:6]': '[*]-[C:1](=[O:2])[C:3]([H:4])([H:5])[H:6]',
                 '[N:1]([C:2]([H:4])([H:5])[H:6])[H:3]': '[*]-[N:1]([C:2]([H:4])([H:5])[H:6])[H:3]',
                 '[N:1]([H:2])[H:3]': '[*]-[N:1]([H:2])[H:3]',
                 '[C:1](=[O:2])[N:3]([C:4])': '[*]-[C:1](=[O:2])[N:3](-[H])[C:4](-[*])(-[*])(-[*])',
                 '[S:1][S:2]': '[*]-[S:1]-[S:2]-[*]',
                 }
for res_name, subs in substructure_dictionary.items():
    res_dict = OrderedDict()
    res_num = 0
    for substruct, atom_names in subs.items():
        rdmol = Chem.MolFromSmarts(substruct)
        old_rdmol = deepcopy(rdmol)
        query = tuple()
        add_port = []
        add_Hs = []
        carboxyl_C = -1
        for query, add_port, add_Hs, carboxyl_C in amino_acid_subs:
            aa_match = rdmol.GetSubstructMatch(query)
            if aa_match:
                break
        if not aa_match and substruct not in special_cases.keys():
            raise Exception(f"{res_name}, {substruct} has no match")
        if substruct in special_cases.keys():
            rdmol = Chem.MolFromSmarts(special_cases[substruct])
        else:
            # add * where appropriate
            for port_idx in add_port:
                rdmol = add_atom(rdmol, aa_match[port_idx], "[*]")
            # add H where appropriate
            for H_idx in add_Hs:
                rdmol = add_atom(rdmol, aa_match[H_idx], "[H]")
            # add port to carboxyl atom if needed
            carb = rdmol.GetAtomWithIdx(aa_match[carboxyl_C])
            if carb.GetDegree() == 2:
                rdmol = add_atom(rdmol, aa_match[carboxyl_C], "[*]")
        formated_rdmol= _fill_out_query(rdmol)
        # In this document, the only new atoms I add are *, so edit the atom names
        num_new_atoms = formated_rdmol.GetNumAtoms() - len(atom_names)
        atom_names = atom_names + ["*"]*num_new_atoms
        sub_smarts = Chem.MolToSmarts(formated_rdmol)
        res_num += 1

        sub_smarts = sub_smarts.replace('&', '')
        res_dict[sub_smarts] = atom_names

    new_subs_dict[res_name] = res_dict
        
output_file = get_data_file_path("proteins/T4-protein.sdf").replace(
    "T4-protein.sdf", "aa_residues_substructures_explicit_bond_orders_with_caps_and_explicit_connectivity.json"
)

with open(output_file, 'w') as f:
    json.dump(new_subs_dict, f, indent=4)
