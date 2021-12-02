import uuid
from io import StringIO

from openmm import unit

try:
    from nglview import Trajectory as _NGLViewTrajectory
except ImportError:
    _NGLViewTrajectory = object


class _OFFTrajectoryNGLView(_NGLViewTrajectory):
    """
    Handling conformers of an OpenFF Molecule as frames in a trajectory. Only
    to be used for NGLview visualization.

    Parameters
    ----------
    molecule : openff.toolkit.topology.Molecule
        The molecule (with conformers) to visualize
    """

    def __init__(self, molecule):
        self.molecule = molecule
        self.ext = "pdb"
        self.params = {}
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index):
        return self.molecule.conformers[index].value_in_unit(unit.angstrom)

    @property
    def n_frames(self):
        return len(self.molecule.conformers)

    def get_structure_string(self):
        memfile = StringIO()
        self.molecule.to_file(memfile, "pdb")
        memfile.seek(0)
        block = memfile.getvalue()
        # FIXME: Prevent multi-model PDB export with a keyword in molecule.to_file()?
        models = block.split("END\n", 1)
        return models[0]
