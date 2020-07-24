import pytest

from openforcefield.topology import Molecule
from openforcefield.cli.generate_conformers import generate_conformers
from openforcefield.tests.utils import get_data_file_path


class TestCLIFunctions:

    @pytest.mark.parametrize('toolkit', ['rdkit', 'openeye'])
    def test_generate_conformers_function(self, toolkit):
        """
        loading one molecule from a .sdf file WITH charges
        (OE only) loading one molecule from a .mol2 file WITH charges
        (OE only) loading one molecule from a .mol2 file WITHOUT charges (is this even possible)
        loading a molecule with a defined name (from any format) and ensuring the output file has that prefix
        loading a molecule with a defined name (from any format), and providing a -f option and ensuring the output file has the -f prefix
        loading a molecule with a defined name (from any format) and ensuring the output file has that prefix
        loading multiple molecules from a .sdf file
        loading one molecule from SMILES in a .smi file
        loading multiple molecules from SMILES in a .smi file
        loading a molecule with ambiguous stereo from SMILES and enumerating stereoisomers
        loading a molecule with defined stereo from SMILES and preserving that stereochemistry

        """
        # loading one molecule from a .sdf file WITHOUT charges
        ethanol = get_data_file_path('molecules/ethanol.sdf')
        assert Molecule.from_file(ethanol).partial_charges is None

        generate_conformers(molecule=ethanol, forcefield='openff-1.0.0.offxml', toolkit=toolkit)

        # loading one molecule from a .sdf file WITH charges
        ethanol_partial_charges = get_data_file_path('molecules/ethanol_partial_charges.sdf')
        assert Molecule.from_file(ethanol_partial_charges).partial_charges is not None

        generate_conformers(molecule=ethanol_partial_charges, forcefield='openff-1.0.0.offxml', toolkit=toolkit)

        if toolkit == 'openeye':
            # (OE only) loading one molecule from a .mol2 file WITH charges
            toluene_partial_charges = get_data_file_path('molecules/toluene_charged.mol2')
            generate_conformers(molecule=toluene_partial_charges, forcefield='openff-1.0.0.offxml', toolkit=toolkit)

            # (OE only) loading one molecule from a .mol2 file WITHOUT charges (is this even possible)

        # loading a molecule with a defined name (from any format) and ensuring the output file has that prefix
        # loading a molecule with a defined name (from any format), and providing a -f option and ensuring the output file has the -f prefix
        # loading a molecule with a defined name (from any format) and ensuring the output file has that prefix
        # loading multiple molecules from a .sdf file
        # loading one molecule from SMILES in a .smi file
        ebastine = get_data_file_path('molecules/ebastine.smi')
        generate_conformers(molecule=ebastine, forcefield='openff-1.0.0.offxml', toolkit=toolkit)

        # TODO:
        # loading multiple molecules from SMILES in a .smi file
        dyes = get_data_file_path('molecules/dyes.smi')
        generate_conformers(molecule=dyes, forcefield='openff-1.0.0.offxml', toolkit=toolkit)

        # loading a molecule with ambiguous stereo from SMILES and enumerating stereoisomers
        # loading a molecule with defined stereo from SMILES and preserving that stereochemistry