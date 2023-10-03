import io
import uuid

import nglview
from openff.toolkit import Topology


def visualize_topology(
    topology: Topology,
    ensure_correct_connectivity: bool = False,
) -> nglview.NGLWidget:
    """
    Visualize the trajectory with NGLView.

    NGLView is a 3D molecular visualization library for use in Jupyter
    notebooks. Note that for performance reasons, by default the
    visualized connectivity is inferred from positions and may not reflect
    the connectivity in the ``Topology``.

    https://github.com/openforcefield/openff-toolkit/pull/1623/files

    Parameters
    ==========
    ensure_correct_connectivity
        If ``True``, the visualization will be guaranteed to reflect the
        connectivity in the ``Topology``. Note that this will severely
        degrade performance, especially for topologies with many atoms.

    Examples
    ========
    Visualize a complex PDB file
    >>> from openff.toolkit import Topology
    >>> from openff.toolkit.utils.utils import get_data_file_path
    >>> pdb_filename = get_data_file_path("systems/test_systems/T4_lysozyme_water_ions.pdb")
    >>> topology = Topology.from_pdb(pdb_filename)
    >>> topology.visualize()
    """

    class TopologyNGLViewStructure(nglview.base_adaptor.Structure):
        """
        OpenFF Topology adaptor.

        https://github.com/openforcefield/openff-toolkit/pull/1623/files

        Parameters
        ----------
        topology
            The `Topology` object to display.
        ext
            The file extension to use to communicate with NGLView. The format must
            be supported for export by the Toolkit via the `Topology.to_file()
            <openff.toolkit.topology.Topology.to_file>` method, and import by
            NGLView. File formats supported by NGLView can be found at
            https://nglviewer.org/ngl/api/manual/file-formats.html

        Example
        -------
        >>> import nglview as nv
        >>>
        >>> top = Topology.from_pdb(pdb_filename)
        >>> nv.NGLWidget(TopologyNGLViewStructure(top))
        """

        def __init__(
            self,
            topology: Topology,
            ext: str = "PDB",
        ):
            self.topology = topology
            self.ext = ext.lower()
            self.params = {}
            self.id = str(uuid.uuid4())

        def get_structure_string(self):
            with io.StringIO() as f:
                self.topology.to_file(f, file_format=self.ext)
                structure_string = f.getvalue()
            return structure_string

    widget = nglview.NGLWidget(
        TopologyNGLViewStructure(
            topology,
            ext="sdf" if ensure_correct_connectivity else "pdb",
        ),
        representations=[
            dict(type="unitcell", params=dict()),
        ],
    )

    widget.add_representation("line", sele="water")
    widget.add_representation("spacefill", sele="ion")
    widget.add_representation("cartoon", sele="protein")
    widget.add_representation(
        "licorice",
        sele="not water and not ion and not protein",
        radius=0.25,
        multipleBond=bool(ensure_correct_connectivity),
    )

    return widget


def visualize_protein_ligand(
    filename: str,
    topology: Topology,
):
    """Assumes the protein is first, ligand second, does not visualize other molecules."""
    import mdtraj
    import nglview

    traj = mdtraj.load(
        filename,
        top=mdtraj.Topology.from_openmm(topology.to_openmm()),
    )

    off_atom_to_mdtraj = lambda atom: traj.topology.atom(topology.atom_index(atom))

    ligand_atom_indices = [
        off_atom_to_mdtraj(atom) for atom in topology.molecule(0).atoms
    ]

    traj.image_molecules(
        anchor_molecules=[ligand_atom_indices],
        inplace=True,
    )

    view = nglview.show_mdtraj(traj)

    view.add_line(selection="water", radius=0.1)
    return view
