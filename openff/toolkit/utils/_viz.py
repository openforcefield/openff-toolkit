import math
import uuid
from io import StringIO
from typing import TYPE_CHECKING, TextIO

from nglview.base_adaptor import Structure, Trajectory
from openff.units import unit

if TYPE_CHECKING:
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
        The ``Molecule`` object to display.
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
        molecule: "Molecule",
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
        return self.molecule.conformers[index].m_as("angstrom")

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

    Communicates with NGLView via PDB, using RDKit to write redundant CONECT
    records indicating multiple bonds. If RDKit is unavailable, falls back
    to ``Topology.to_file``.

    Parameters
    ----------
    topology
        The ``Topology`` object to display.

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
    ):
        self.topology = topology
        self.ext = "pdb"
        self.params: dict = dict()
        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        with StringIO() as f:
            if RDKIT_AVAILABLE:
                from rdkit.Chem.rdmolfiles import PDBWriter  # type: ignore

                write_box_vectors(f, self.topology)

                writer: PDBWriter = PDBWriter(f)
                for mol in self.topology.molecules:
                    writer.write(mol.to_rdkit())

                writer.close()

                structure_string = f.getvalue()
            else:
                self.topology.to_file(f, file_format="pdb")
                structure_string = f.getvalue()

        return structure_string


def write_box_vectors(file_obj: TextIO, topology: Topology):
    if topology.box_vectors is not None:
        a, b, c = topology.box_vectors.m_as(unit.nanometer)
        a_length = a.norm()
        b_length = b.norm()
        c_length = c.norm()

        alpha = math.acos(b.dot(c) / (b_length * c_length))
        beta = math.acos(c.dot(a) / (c_length * a_length))
        gamma = math.acos(a.dot(b) / (a_length * b_length))

        RAD_TO_DEG = 180 / math.pi
        print(
            "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 "
            % (
                a_length * 10,
                b_length * 10,
                c_length * 10,
                alpha * RAD_TO_DEG,
                beta * RAD_TO_DEG,
                gamma * RAD_TO_DEG,
            ),
            file=file_obj,
        )
