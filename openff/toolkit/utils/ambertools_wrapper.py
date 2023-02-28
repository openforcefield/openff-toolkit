"""
Wrapper class providing a minimal consistent interface to `AmberTools <http://ambermd.org/AmberTools.php>`_.
"""

__all__ = ("AmberToolsToolkitWrapper",)

import subprocess
import tempfile
from collections import defaultdict
from shutil import which
from typing import TYPE_CHECKING, Dict, List, Optional, Union

import numpy as np
from openff.units import Quantity, unit

from openff.toolkit.utils import base_wrapper, rdkit_wrapper
from openff.toolkit.utils.exceptions import (
    AntechamberNotFoundError,
    ChargeCalculationError,
    ChargeMethodUnavailableError,
    ToolkitUnavailableException,
)
from openff.toolkit.utils.utils import temporary_cd

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import Molecule


class AmberToolsToolkitWrapper(base_wrapper.ToolkitWrapper):
    """
    AmberTools toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "AmberTools"
    _toolkit_installation_instructions = (
        "The AmberTools toolkit (free and open source) can be found at "
        "https://anaconda.org/conda-forge/ambertools"
    )

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = []
        self._toolkit_file_write_formats = []

        if not self.is_available():
            raise ToolkitUnavailableException(
                f"The required toolkit {self._toolkit_name} is not "
                f"available. {self._toolkit_installation_instructions}"
            )

        # TODO: More reliable way to extract AmberTools version
        out = subprocess.check_output(["antechamber", "-L"])
        ambertools_version = out.decode("utf-8").split("\n")[1].split()[3].strip(":")
        self._toolkit_version = ambertools_version

        # TODO: Find AMBERHOME or executable home, checking miniconda if needed
        # Store an instance of an RDKitToolkitWrapper for file I/O
        self._rdkit_toolkit_wrapper = rdkit_wrapper.RDKitToolkitWrapper()

    @staticmethod
    def is_available() -> bool:
        """
        Check whether the AmberTools toolkit is installed

        Returns
        -------
        is_installed : bool
            True if AmberTools is installed, False otherwise.

        """
        # TODO: Check all tools needed
        ANTECHAMBER_PATH = which("antechamber")
        if ANTECHAMBER_PATH is None:
            return False
        # AmberToolsToolkitWrapper needs RDKit to do basically anything, since its interface requires SDF I/O
        if not (rdkit_wrapper.RDKitToolkitWrapper.is_available()):
            return False
        return True

    def assign_partial_charges(
        self,
        molecule: "Molecule",
        partial_charge_method: Optional[str] = None,
        use_conformers: Optional[List[Quantity]] = None,
        strict_n_conformers: bool = False,
        normalize_partial_charges: bool = True,
        _cls=None,
    ):
        """
        Compute partial charges with AmberTools using antechamber/sqm, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The charge model to use. One of ['gasteiger', 'am1bcc', 'am1-mulliken'].
            If None, 'am1-mulliken' will be used.
        use_conformers : iterable of unit-wrapped numpy arrays, each
            with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            List of unit-wrapped numpy arrays to use for partial charge calculation.
            If None, an appropriate number of conformers will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for
            the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        normalize_partial_charges : bool, default=True
            Whether to offset partial charges so that they sum to the total formal charge of the molecule.
            This is used to prevent accumulation of rounding errors when the partial charge generation method has
            low precision.
        _cls : class
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        ChargeCalculationError if the charge method is supported by this toolkit, but fails
        """

        import os
        import subprocess

        from openff.toolkit.topology import Molecule

        if partial_charge_method is None:
            partial_charge_method = "am1-mulliken"
        else:
            # Standardize method name for string comparisons
            partial_charge_method = partial_charge_method.lower()

        SUPPORTED_CHARGE_METHODS: Dict[str, Dict[str, Union[int, str]]] = {
            "am1bcc": {
                "antechamber_keyword": "bcc",
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "am1-mulliken": {
                "antechamber_keyword": "mul",
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "gasteiger": {
                "antechamber_keyword": "gas",
                "min_confs": 0,
                "max_confs": 0,
                "rec_confs": 0,
            },
        }

        if partial_charge_method not in SUPPORTED_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f"partial_charge_method '{partial_charge_method}' is not available from AmberToolsToolkitWrapper. "
                f"Available charge methods are {list(SUPPORTED_CHARGE_METHODS.keys())} "
            )

        charge_method = SUPPORTED_CHARGE_METHODS[partial_charge_method]

        if _cls is None:
            _cls = Molecule

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = _cls(molecule)

        if use_conformers is None:
            if charge_method["rec_confs"] == 0:
                mol_copy._conformers = None
            else:
                mol_copy.generate_conformers(
                    n_conformers=charge_method["rec_confs"],
                    rms_cutoff=0.25 * unit.angstrom,
                    toolkit_registry=rdkit_wrapper.RDKitToolkitWrapper(),
                )
            # TODO: What's a "best practice" RMS cutoff to use here?
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy._add_conformer(conformer)
            self._check_n_conformers(
                mol_copy,
                partial_charge_method=partial_charge_method,
                min_confs=charge_method["min_confs"],  # type: ignore[arg-type]
                max_confs=charge_method["max_confs"],  # type: ignore[arg-type]
                strict_n_conformers=strict_n_conformers,
            )

        ANTECHAMBER_PATH = which("antechamber")
        if ANTECHAMBER_PATH is None:
            raise AntechamberNotFoundError(
                "Antechamber not found, cannot run charge_mol()"
            )

        # Compute charges
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = mol_copy.total_charge.m_as(unit.elementary_charge)
                # Write out molecule in SDF format
                # TODO: How should we handle multiple conformers?
                self._rdkit_toolkit_wrapper.to_file(
                    mol_copy, "molecule.sdf", file_format="sdf"
                )
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                short_charge_method = charge_method["antechamber_keyword"]
                subprocess.check_output(
                    [
                        "antechamber",
                        "-i",
                        "molecule.sdf",
                        "-fi",
                        "sdf",
                        "-o",
                        "charged.mol2",
                        "-fo",
                        "mol2",
                        "-pf",
                        "yes",
                        "-dr",
                        "n",
                        "-c",
                        str(short_charge_method),
                        "-nc",
                        str(net_charge),
                    ]
                )
                # Write out just charges
                subprocess.check_output(
                    [
                        "antechamber",
                        "-dr",
                        "n",
                        "-i",
                        "charged.mol2",
                        "-fi",
                        "mol2",
                        "-o",
                        "charges2.mol2",
                        "-fo",
                        "mol2",
                        "-c",
                        "wc",
                        "-cf",
                        "charges.txt",
                        "-pf",
                        "yes",
                    ]
                )
                # Check to ensure charges were actually produced
                if not os.path.exists("charges.txt"):
                    # TODO: copy files into local directory to aid debugging?
                    raise ChargeCalculationError(
                        "Antechamber/sqm partial charge calculation failed on "
                        "molecule {} (SMILES {})".format(
                            molecule.name, molecule.to_smiles()
                        )
                    )
                # Read the charges
                with open("charges.txt", "r") as infile:
                    contents = infile.read()
                text_charges = contents.split()
                charges = np.zeros([molecule.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)
                # TODO: Ensure that the atoms in charged.mol2 are in the same order as in molecule.sdf
        charges = unit.Quantity(charges, unit.elementary_charge)
        molecule.partial_charges = charges

        if normalize_partial_charges:
            molecule._normalize_partial_charges()

    def _modify_sqm_in_to_request_bond_orders(self, file_path):
        """
        Modify a sqm.in file produced by antechamber to include the "printbondorders=1" directive
        in the header. This method will overwrite the original file.

        Parameters
        ----------
        file_path : str
            The path to sqm.in
        """

        data = open(file_path).read()

        # Original sqm.in file headerlooks like:

        # Run semi-empirical minimization
        #  &qmmm
        #    qm_theory='AM1', grms_tol=0.0005,
        #    scfconv=1.d-10, ndiis_attempts=700,   qmcharge=0,
        #  /
        # ... (atom coordinates in something like XYZ format) ...

        # To get WBOs, we need to add "printbondorders=1" to the list of keywords

        # First, split the sqm.in text at the "/" mark at the end of the header
        datasp = data.split("/")
        # Insert the "printbondorders" directive in a new line and re-add the "/"
        datasp.insert(1, "printbondorders=1, \n /")
        # Reassemble the file text
        new_data = "".join(datasp)
        # Write the new file contents, overwriting the original file.
        with open(file_path, "w") as of:
            of.write(new_data)

    def _get_fractional_bond_orders_from_sqm_out(
        self, file_path, validate_elements=None
    ):
        """
        Process a SQM output file containing bond orders, and return a dict of the form
        dict[atom_1_index, atom_2_index] = fractional_bond_order

        Parameters
        ----------
        file_path : str
            File path for sqm output file
        validate_elements : iterable of str
            The element symbols expected in molecule index order. A ValueError will be raised
            if the elements are not found in this order.

        Returns
        -------
        bond_orders : dict[(int, int)]: float
            A dictionary where the keys are tuples of two atom indices and the values are
            floating-point bond orders. The keys are sorted in ascending order, such that
            the lower atom index is key[0] and the higher is key[1].
        """

        # Example sqm.out section with WBOs:
        #  Bond Orders
        #
        #   QMMM:    NUM1 ELEM1 NUM2 ELEM2      BOND_ORDER
        #   QMMM:       2   C      1   C        1.41107532
        #   QMMM:       3   C      1   C        1.41047804
        # ...
        #   QMMM:      15   H     13   H        0.00000954
        #   QMMM:      15   H     14   H        0.00000813
        #
        #            --------- Calculation Completed ----------

        data = open(file_path).read()

        begin_sep = """ Bond Orders
 
  QMMM:    NUM1 ELEM1 NUM2 ELEM2      BOND_ORDER
"""
        end_sep = """

           --------- Calculation Completed ----------
"""
        # Extract the chunk of text between begin_sep and end_sep, and split it by newline
        fbo_lines = data.split(begin_sep)[1].split(end_sep)[0].split("\n")

        # Iterate over the lines and populate the dict to return
        bond_orders = dict()
        for line in fbo_lines:
            linesp = line.split()
            atom_index_1 = int(linesp[1])
            atom_element_1 = linesp[2]
            atom_index_2 = int(linesp[3])
            atom_element_2 = linesp[4]
            bond_order = float(linesp[5])

            # If validate_elements was provided, ensure that the ordering of element symbols is what we expected
            if validate_elements is not None:
                if (atom_element_1 != validate_elements[atom_index_1 - 1]) or (
                    atom_element_2 != validate_elements[atom_index_2 - 1]
                ):
                    # raise ValueError('\n'.join(fbo_lines))
                    raise ValueError(
                        f"Elements or indexing in sqm output differ from expectation. "
                        f"Expected {validate_elements[atom_index_1]} with index {atom_index_1} and "
                        f"{validate_elements[atom_index_2]} with index {atom_index_2}, "
                        f"but SQM output has {atom_element_1} and {atom_element_2} for the same atoms."
                    )

            # To make lookup easier, we identify bonds as integer tuples with the lowest atom index
            # first and the highest second.
            index_tuple = tuple(sorted([atom_index_1, atom_index_2]))
            bond_orders[index_tuple] = bond_order
        return bond_orders

    def assign_fractional_bond_orders(
        self,
        molecule: "Molecule",
        bond_order_model: Optional[str] = None,
        use_conformers: Optional[List[str]] = None,
        _cls=None,
    ):
        """
        Update and store list of bond orders this molecule. Bond orders are stored on each
        bond, in the `bond.fractional_bond_order` attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.molecule Molecule
            The molecule to assign wiberg bond orders to
        bond_order_model : str, optional, default=None
            The charge model to use. Only allowed value is 'am1-wiberg'. If None, 'am1-wiberg' will be used.
        use_conformers : iterable of unit-wraapped np.array with shape (n_atoms, 3)
            and dimension of distance, optional, default=None
            The conformers to use for fractional bond order calculation. If None, an appropriate
            number of conformers will be generated by an available ToolkitWrapper.
        _cls : class
            Molecule constructor
        """
        from openff.toolkit.topology import Molecule

        ANTECHAMBER_PATH = which("antechamber")
        if ANTECHAMBER_PATH is None:
            raise AntechamberNotFoundError(
                "Antechamber not found, cannot run "
                "AmberToolsToolkitWrapper.assign_fractional_bond_orders()"
            )

        if _cls is None:
            _cls = Molecule

        # Make a copy since we'll be messing with this molecule's conformers
        temp_mol = _cls(molecule)

        if use_conformers is None:
            temp_mol.generate_conformers(
                n_conformers=1,
                toolkit_registry=self._rdkit_toolkit_wrapper,
            )
        else:
            temp_mol._conformers = None
            for conformer in use_conformers:
                temp_mol._add_conformer(conformer)

        if len(temp_mol.conformers) == 0:
            raise ValueError(
                "No conformers present in molecule submitted for fractional bond order calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.assign_fractional_bond_orders"
            )

        # Compute bond orders
        bond_order_model_to_antechamber_keyword = {"am1-wiberg": "mul"}
        supported_bond_order_models = list(
            bond_order_model_to_antechamber_keyword.keys()
        )
        if bond_order_model is None:
            bond_order_model = "am1-wiberg"

        bond_order_model = bond_order_model.lower()

        if bond_order_model not in supported_bond_order_models:
            raise ValueError(
                f"Bond order model '{bond_order_model}' is not supported by AmberToolsToolkitWrapper. "
                f"Supported models are {supported_bond_order_models}"
            )
        ac_charge_keyword = bond_order_model_to_antechamber_keyword[bond_order_model]

        bond_orders = defaultdict(list)

        for conformer in [*temp_mol.conformers]:
            with tempfile.TemporaryDirectory() as tmpdir:
                with temporary_cd(tmpdir):
                    net_charge = temp_mol.total_charge
                    # Write out molecule in SDF format
                    temp_mol._conformers = [conformer]
                    self._rdkit_toolkit_wrapper.to_file(
                        temp_mol, "molecule.sdf", file_format="sdf"
                    )
                    # Prepare sqm.in file as if we were going to run charge calc
                    # TODO: Add error handling if antechamber chokes
                    subprocess.check_output(
                        [
                            "antechamber",
                            "-i",
                            "molecule.sdf",
                            "-fi",
                            "sdf",
                            "-o",
                            "sqm.in",
                            "-fo",
                            "sqmcrt",
                            "-pf",
                            "yes",
                            "-c",
                            ac_charge_keyword,
                            "-nc",
                            str(net_charge),
                        ]
                    )
                    # Modify sqm.in to request bond order calculation
                    self._modify_sqm_in_to_request_bond_orders("sqm.in")
                    # Run sqm to get bond orders
                    subprocess.check_output(
                        ["sqm", "-i", "sqm.in", "-o", "sqm.out", "-O"]
                    )
                    # Ensure that antechamber/sqm did not change the indexing by checking against
                    # an ordered list of element symbols for this molecule
                    expected_elements = [atom.symbol for atom in molecule.atoms]
                    conformer_bond_orders = (
                        self._get_fractional_bond_orders_from_sqm_out(
                            "sqm.out", validate_elements=expected_elements
                        )
                    )

                    for bond_indices, value in conformer_bond_orders.items():
                        bond_orders[bond_indices].append(value)

        # Note that sqm calculate WBOs for ALL PAIRS of atoms, not just those that have
        # bonds defined in the original molecule. So here we iterate over the bonds in
        # the original molecule and only nab the WBOs for those.
        for bond in molecule.bonds:
            # The atom index tuples that act as bond indices are ordered from lowest to highest by
            # _get_fractional_bond_orders_from_sqm_out, so here we make sure that we look them up in
            # sorted order as well
            sorted_atom_indices = sorted(
                tuple([bond.atom1_index + 1, bond.atom2_index + 1])
            )
            bond.fractional_bond_order = np.mean(
                bond_orders[tuple(sorted_atom_indices)]
            )
