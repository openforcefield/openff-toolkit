import uuid
from io import StringIO

from nglview.base_adaptor import Structure, Trajectory
from openff.units import unit

from openff.toolkit import Molecule, Topology

MOLECULE_DEFAULT_REPS = [
    dict(type="licorice", params=dict(radius=0.25, multipleBond=True))
]


class MoleculeNGLViewTrajectory(Structure, Trajectory):
    """
    OpenFF Molecule adaptor.

    The file format used to communicate with NGLView can make a big difference
    in the visualization. For instance, using PDB will require NGLView to infer
    bonds from positions, whereas using MOL2 or SDF may result in a loss of
    residue data.

    Parameters
    ----------
    molecule
        The `Molecule` object to display.
    ext
        The file extension to use to communicate with NGLView. The format must
        be supported for export by the Toolkit via the `Molecule.to_file()
        <openff.toolkit.topology.molecule.to_file>` method, and import by
        NGLView. File formats supported by NGLView can be found at
        https://nglviewer.org/ngl/api/manual/file-formats.html

    Example
    -------
    >>> import nglview as nv
    >>>
    >>> mol = Molecule.from_polymer_pdb("file.pdb")  # doctest: +SKIP
    >>> nv.NGLWidget(MoleculeNGLViewTrajectory(mol))  # doctest: +SKIP
    """

    def __init__(
        self,
        molecule: Molecule,
        ext: str = "MOL2",
    ):
        if not molecule.conformers:
            raise ValueError(
                "Cannot visualize a molecule without conformers with NGLView"
            )
        self.molecule = molecule
        self.ext = ext.lower()
        self.params: dict = dict()
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index: int = 0):
        return self.molecule.conformers[index].m_as(unit.angstrom)

    @property
    def n_frames(self):
        return len(self.molecule.conformers)

    def get_structure_string(self):
        with StringIO() as f:
            self.molecule.to_file(f, file_format=self.ext)
            structure_string = f.getvalue()
        return structure_string


class TopologyNGLViewStructure(Structure):
    """
    OpenFF Topology adaptor.

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
    >>> top = Topology.from_pdb("file.pdb")  # doctest: +SKIP
    >>> nv.NGLWidget(TopologyNGLViewStructure(top))  # doctest: +SKIP
    """

    def __init__(
        self,
        topology: Topology,
        ext: str = "PDB",
    ):
        self.topology = topology
        self.ext = ext.lower()
        self.params: dict = dict()
        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        with StringIO() as f:
            self.topology.to_file(f, file_format=self.ext)
            structure_string = f.getvalue()
        return structure_string
