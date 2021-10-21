from openff.toolkit.utils._cif_to_substructure_dict import CifSubstructures
from openff.toolkit.utils import get_data_file_path


import os
if not os.path.exists('aa-variants-v1.cif'):
    import requests
    r = requests.get('https://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif')
    print(r.ok)
    if r.ok:
        with open('aa-variants-v1.cif', 'wb') as of:
            for chunk in r.iter_content(chunk_size=1024 * 8):
                if chunk:
                    of.write(chunk)
                    of.flush()
                    os.fsync(of.fileno())


cif_object = CifSubstructures()
cif_object.from_file('aa-variants-v1.cif',
                     replace_quadruple_bond_with_any=True,
                     remove_charge_bond_order_resonant=True)

# Automatically patch known problems - better that this explodes when things are fixed
cif_object._add_common_substructures()

output_file = get_data_file_path('proteins/T4-protein.sdf').replace('T4-protein.sdf', 'aa_residues_substructures_with_caps.json')
cif_object.to_json_file(output_file)
