import os

from _cif_to_substructure_dict import CifSubstructures

from openff.toolkit.utils.utils import get_data_file_path

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

# Make Old substructure file to keep some functions working
output_file = get_data_file_path("proteins/T4-protein.sdf").replace(
    "T4-protein.sdf", "aa_residues_substructures_explicit_bond_orders_with_caps.json"
)
cif_object.to_json_file(output_file)

# New substructure file with explicit connectivity and neighbor atom spec
cif_object._reformat_substructs()

output_file = get_data_file_path("proteins/T4-protein.sdf").replace(
    "T4-protein.sdf",
    "aa_residues_substructures_explicit_bond_orders_with_caps_explicit_connectivity.json",
)
cif_object.to_json_file(output_file)
