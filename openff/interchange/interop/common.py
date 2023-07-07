"""Utilities for interoperability with multiple packages."""

from openff.interchange import Interchange
from openff.interchange.models import VirtualSiteKey


def _build_typemap(interchange: "Interchange") -> dict[int, str]:
    typemap = dict()
    elements: dict[str, int] = dict()

    # TODO: Think about how this logic relates to atom name/type clashes
    for atom_index, atom in enumerate(interchange.topology.atoms):
        element_symbol = atom.symbol
        # TODO: Use this key to condense, see parmed.openmm._process_nobonded
        # parameters = _get_lj_parameters([*parameters.values()])
        # key = tuple([*parameters.values()])

        if element_symbol not in elements.keys():
            elements[element_symbol] = 1
        else:
            elements[element_symbol] += 1

        atom_type = f"{element_symbol}{elements[element_symbol]}"
        typemap[atom_index] = atom_type

    return typemap


def _build_virtual_site_map(interchange: "Interchange") -> dict[VirtualSiteKey, int]:
    """
    Construct a mapping between the VirtualSiteKey objects found in a SMIRNOFFVirtualSiteHandler and particle indices.
    """
    virtual_site_topology_index_map: dict[VirtualSiteKey, int] = dict()

    if "VirtualSites" not in interchange.collections:
        return virtual_site_topology_index_map

    n_atoms = interchange.topology.n_atoms

    for index, virtual_site_key in enumerate(
        interchange["VirtualSites"].key_map.keys(),
    ):
        assert isinstance(virtual_site_key, VirtualSiteKey)
        virtual_site_topology_index_map[virtual_site_key] = n_atoms + 1 + index

    return virtual_site_topology_index_map
