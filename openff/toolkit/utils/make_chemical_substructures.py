from openff.toolkit.utils._cif_to_substructure_dict import CifSubstructures
from openff.toolkit.utils import get_data_file_path

cif_object = CifSubstructures()
cif_object.from_file('/home/ijpulidos/workdir/data/aa-variants-v1.cif', include_leaving=False, replace_quadruple_bond_with_any=False, remove_charge_bond_order_resonant=False)

# Automatically patch known problems - better that this explodes when things are fixed
cif_object._patch_known_problems()

output_file = get_data_file_path('proteins/T4-protein.sdf').replace('T4-protein.sdf', 'aa_variants_chemical_substructures.json')
cif_object.to_json_file(output_file)
