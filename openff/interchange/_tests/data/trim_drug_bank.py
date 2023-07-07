"""Script to generate a subset of MiniDrugBank."""
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import get_data_file_path
from rdkit import Chem

trimmed_mol_names = [
    "DrugBank_5354",
    "DrugBank_2800",
    "DrugBank_5418",
    "DrugBank_104",
    "DrugBank_5516",
    "DrugBank_5523",
    "DrugBank_2967",
    "DrugBank_246",
    "DrugBank_2991",
    "DrugBank_3028",
    "DrugBank_3046",
    "DrugBank_3087",
    "DrugBank_390",
    "DrugBank_5804",
    "DrugBank_3346",
    "DrugBank_3358",
    "DrugBank_3406",
    "DrugBank_5900",
    "DrugBank_5902",
    "DrugBank_3479",
    "DrugBank_3503",
    "DrugBank_3547",
    "DrugBank_3565",
    "DrugBank_6032",
    "DrugBank_914",
    "DrugBank_977",
    "DrugBank_6182",
    "DrugBank_3817",
    "DrugBank_3954",
    "DrugBank_6355",
    "DrugBank_4074",
    "DrugBank_1448",
    "DrugBank_1449",
    "DrugBank_4138",
    "DrugBank_4161",
    "DrugBank_1538",
    "DrugBank_6531",
    "DrugBank_1564",
    "DrugBank_6533",
    "DrugBank_4215",
    "DrugBank_4217",
    "DrugBank_1598",
    "DrugBank_1608",
    "DrugBank_4249",
    "DrugBank_1637",
    "DrugBank_1661",
    "DrugBank_6647",
    "DrugBank_4323",
    "DrugBank_1721",
    "DrugBank_1722",
    "DrugBank_1742",
    "DrugBank_6722",
    "DrugBank_4468",
    "DrugBank_4515",
    "DrugBank_6865",
    "DrugBank_4580",
    "DrugBank_1971",
    "DrugBank_4662",
    "DrugBank_7108",
    "DrugBank_7124",
    "DrugBank_2140",
    "DrugBank_2148",
    "DrugBank_2186",
    "DrugBank_2237",
    "DrugBank_4959",
    "DrugBank_2429",
    "DrugBank_5154",
    "DrugBank_2563",
    "DrugBank_2570",
    "DrugBank_2584",
    "DrugBank_2585",
    "DrugBank_2684",
]


mini_drug_bank = Molecule.from_file(
    get_data_file_path("molecules/MiniDrugBank.sdf"),
    allow_undefined_stereo=True,
)
mols = [mol for mol in mini_drug_bank if mol.name in trimmed_mol_names]

with Chem.SDWriter("MiniDrugBankTrimmed.sdf") as f:
    for openff_mol in mols:
        try:
            openff_mol.assign_partial_charges(partial_charge_method="am1bcc")
            rdmol = openff_mol.to_rdkit()
        except:  # noqa
            continue
        f.write(rdmol)
