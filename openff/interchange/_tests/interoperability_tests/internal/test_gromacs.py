import sys
from math import exp

import numpy
import openmm
import parmed
import pytest
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField, VirtualSiteHandler
from openff.toolkit.utils import get_data_file_path
from openff.units import unit
from openff.units.openmm import ensure_quantity
from openff.utilities.testing import skip_if_missing
from openmm import app
from openmm import unit as openmm_unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest, needs_gmx
from openff.interchange.components.nonbonded import BuckinghamvdWCollection
from openff.interchange.components.potentials import Potential
from openff.interchange.drivers import get_gromacs_energies, get_openmm_energies
from openff.interchange.exceptions import (
    GMXMdrunError,
    UnsupportedExportError,
    VirtualSiteTypeNotImplementedError,
)
from openff.interchange.interop.gromacs._import._import import (
    _read_box,
    _read_coordinates,
)
from openff.interchange.models import PotentialKey, TopologyKey

if sys.version_info >= (3, 10):
    from importlib import resources
else:
    import importlib_resources as resources


@needs_gmx
class TestGROMACSGROFile(_BaseTest):
    _INTERMOL_PATH = resources.files(
        "intermol.tests.gromacs.unit_tests",
    )

    @skip_if_missing("intermol")
    def test_load_gro(self):
        file = self._INTERMOL_PATH / "angle10_vacuum/angle10_vacuum.gro"

        positions = _read_coordinates(file)
        box = _read_box(file)

        openmm_gro = openmm.app.GromacsGroFile(file)
        openmm_positions = ensure_quantity(
            openmm_gro.getPositions(),
            "openff",
        )
        openmm_box = ensure_quantity(
            openmm_gro.getPeriodicBoxVectors(),
            "openff",
        )

        assert numpy.allclose(positions, openmm_positions)
        assert numpy.allclose(box, openmm_box)

    @skip_if_missing("intermol")
    def test_load_gro_nonstandard_precision(self):
        file = self._INTERMOL_PATH / "lj3_bulk/lj3_bulk.gro"

        coords = _read_coordinates(file).m_as(unit.nanometer)

        # OpenMM seems to assume a precision of 3. Use InterMol instead here.
        from intermol.gromacs.grofile_parser import GromacsGroParser

        intermol_gro = GromacsGroParser(file)
        intermol_gro.read()

        # InterMol stores positions as NumPy arrays _of openmm.unit.Quantity_
        # objects, so need to carefully convert everything through unitless
        # floats; converting through Pint-like quantities might produce a
        # funky double-wrapped thing
        def converter(x):
            return x.value_in_unit(openmm_unit.nanometer)

        other_coords = numpy.frompyfunc(converter, 1, 1)(intermol_gro.positions).astype(
            float,
        )

        assert numpy.allclose(coords, other_coords)

        # This file happens to have 12 digits of preicion; what really matters is that
        # the convential precision of 3 was not used.
        n_decimals = len(str(coords[0, 0]).split(".")[1])
        assert n_decimals == 12

    def test_vaccum_warning(self, sage):
        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(force_field=sage, topology=[molecule])

        assert out.box is None

        with pytest.warns(UserWarning, match="gitlab"):
            out.to_gro("tmp.gro")

    @pytest.mark.slow()
    def test_residue_info(self, sage):
        """Test that residue information is passed through to .gro files."""
        import mdtraj

        protein = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_HIE.pdb"),
        )

        ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

        out = Interchange.from_smirnoff(
            force_field=ff14sb,
            topology=[protein],
        )

        out.to_gro("tmp.gro")

        mdtraj_topology = mdtraj.load("tmp.gro").topology

        for found_residue, original_residue in zip(
            mdtraj_topology.residues,
            out.topology.hierarchy_iterator("residues"),
        ):
            assert found_residue.name == original_residue.residue_name
            assert str(found_residue.resSeq) == original_residue.residue_number

    def test_atom_names_pdb(self):
        peptide = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_ALA_ALA.pdb"),
        )
        ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

        Interchange.from_smirnoff(ff14sb, peptide.to_topology()).to_gro(
            "atom_names.gro",
        )

        pdb_object = app.PDBFile(get_data_file_path("proteins/MainChain_ALA_ALA.pdb"))
        pdb_atom_names = [atom.name for atom in pdb_object.topology.atoms()]

        openmm_atom_names = app.GromacsGroFile("atom_names.gro").atomNames

        assert openmm_atom_names == pdb_atom_names


@needs_gmx
class TestGROMACS(_BaseTest):
    @pytest.mark.slow()
    @pytest.mark.skip("from_top is not yet refactored for new Topology API")
    @pytest.mark.parametrize("reader", ["intermol", "internal"])
    @pytest.mark.parametrize(
        "smiles",
        [
            "C",
            "O=C=O",  # Adds unconstrained bonds without torsion(s)
            "CC",  # Adds a proper torsion term(s)
            # "C=O",  # Simplest molecule with any improper torsion
            "OC=O",  # Simplest molecule with a multi-term torsion
            # "CCOC",  # This hits t86, which has a non-1.0 idivf
            # "C1COC(=O)O1",  # This adds an improper, i2
        ],
    )
    def test_simple_roundtrip(self, sage, smiles, reader):
        molecule = Molecule.from_smiles(smiles)
        molecule.name = molecule.to_hill_formula()
        molecule.generate_conformers(n_conformers=1)
        topology = molecule.to_topology()

        out = Interchange.from_smirnoff(force_field=sage, topology=topology)
        out.box = [4, 4, 4]
        out.positions = molecule.conformers[0]

        out.to_top("out.top")
        out.to_gro("out.gro")

        converted = Interchange.from_gromacs("out.top", "out.gro", reader=reader)

        assert numpy.allclose(out.positions, converted.positions)
        assert numpy.allclose(out.box, converted.box)

        get_gromacs_energies(out).compare(
            get_gromacs_energies(converted),
            tolerances={
                "Bond": 0.002 * molecule.n_bonds * unit.kilojoule / unit.mol,
                "Electrostatics": 0.05 * unit.kilojoule / unit.mol,
            },
        )

    @skip_if_missing("parmed")
    def test_num_impropers(self, sage):
        top = Molecule.from_smiles("CC1=CC=CC=C1").to_topology()
        out = Interchange.from_smirnoff(sage, top)
        out.box = unit.Quantity(4 * numpy.eye(3), units=unit.nanometer)
        out.to_top("tmp.top")

        # Sanity check; toluene should have some improper(s)
        assert len(out["ImproperTorsions"].key_map) > 0

        struct = parmed.load_file("tmp.top")
        n_impropers_parmed = len([d for d in struct.dihedrals if d.improper])
        assert n_impropers_parmed == len(out["ImproperTorsions"].key_map)

    @pytest.mark.slow()
    @skip_if_missing("intermol")
    @pytest.mark.skip(reason="Re-implement when SMIRNOFF supports more mixing rules")
    def test_set_mixing_rule(self, ethanol_top, sage):
        from intermol.gromacs.gromacs_parser import GromacsParser

        interchange = Interchange.from_smirnoff(force_field=sage, topology=ethanol_top)
        interchange.positions = numpy.zeros((ethanol_top.n_atoms, 3))
        interchange.to_gro("tmp.gro")

        interchange.box = [4, 4, 4]
        interchange.to_top("lorentz.top")
        lorentz = GromacsParser("lorentz.top", "tmp.gro").read()
        assert lorentz.combination_rule == "Lorentz-Berthelot"

        interchange["vdW"].mixing_rule = "geometric"

        interchange.to_top("geometric.top")
        geometric = GromacsParser("geometric.top", "tmp.gro").read()
        assert geometric.combination_rule == "Multiply-Sigeps"

    @pytest.mark.skip(reason="Re-implement when SMIRNOFF supports more mixing rules")
    def test_unsupported_mixing_rule(self, ethanol_top, sage):
        # TODO: Update this test when the model supports more mixing rules than GROMACS does
        interchange = Interchange.from_smirnoff(force_field=sage, topology=ethanol_top)
        interchange["vdW"].mixing_rule = "kong"

        with pytest.raises(UnsupportedExportError, match="rule `geometric` not compat"):
            interchange.to_top("out.top")

    @pytest.mark.slow()
    def test_residue_info(self, sage):
        """Test that residue information is passed through to .top files."""
        import parmed
        from openff.units.openmm import from_openmm

        pdb_path = get_data_file_path("proteins/MainChain_HIE.pdb")

        protein = Molecule.from_polymer_pdb(pdb_path)

        box_vectors = from_openmm(
            openmm.app.PDBFile(pdb_path).topology.getPeriodicBoxVectors(),
        )

        ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

        out = Interchange.from_smirnoff(
            force_field=ff14sb,
            topology=[protein],
            box=box_vectors,
        )

        out.to_top("tmp.top")

        parmed_structure = parmed.load_file("tmp.top")

        for found_residue, original_residue in zip(
            parmed_structure.residues,
            out.topology.hierarchy_iterator("residues"),
        ):
            assert found_residue.name == original_residue.residue_name
            assert str(found_residue.number + 1) == original_residue.residue_number

    @pytest.mark.slow()
    @pytest.mark.skip(
        reason="Update when energy reports intentionally support non-vdW handlers",
    )
    def test_argon_buck(self):
        """Test that Buckingham potentials are supported and can be exported"""
        from openff.interchange.smirnoff._nonbonded import (
            SMIRNOFFElectrostaticsCollection,
        )

        mol = Molecule.from_smiles("[#18]")
        mol.name = "Argon"
        top = Topology.from_molecules([mol, mol])

        # http://www.sklogwiki.org/SklogWiki/index.php/Argon#Buckingham_potential
        erg_mol = unit.erg / unit.mol * float(unit.avogadro_number)
        A = 1.69e-8 * erg_mol
        B = 1 / (0.273 * unit.angstrom)
        C = 102e-12 * erg_mol * unit.angstrom**6

        r = 0.3 * unit.nanometer

        buck = BuckinghamvdWCollection()
        coul = SMIRNOFFElectrostaticsCollection(method="pme")

        pot_key = PotentialKey(id="[#18]")
        pot = Potential(parameters={"A": A, "B": B, "C": C})

        for atom in top.atoms:
            top_key = TopologyKey(atom_indices=(top.atom_index(atom),))
            buck.key_map.update({top_key: pot_key})

            coul.key_map.update({top_key: pot_key})
            coul.potentials.update(
                {pot_key: Potential(parameters={"charge": 0 * unit.elementary_charge})},
            )

        for molecule in top.molecules:
            molecule.partial_charges = unit.Quantity(
                molecule.n_atoms * [0],
                unit.elementary_charge,
            )

        buck.potentials[pot_key] = pot

        out = Interchange()
        out.collections["Buckingham-6"] = buck
        out.collections["Electrostatics"] = coul
        out.topology = top
        out.box = [10, 10, 10] * unit.nanometer
        out.positions = [[0, 0, 0], [0.3, 0, 0]] * unit.nanometer
        out.to_gro("out.gro", writer="internal")
        out.to_top("out.top", writer="internal")

        omm_energies = get_openmm_energies(out, combine_nonbonded_forces=True)
        by_hand = A * exp(-B * r) - C * r**-6

        resid = omm_energies.energies["vdW"] - by_hand
        assert resid < 1e-5 * unit.kilojoule / unit.mol

        # TODO: Add back comparison to GROMACS energies once GROMACS 2020+
        # supports Buckingham potentials
        with pytest.raises(GMXMdrunError):
            get_gromacs_energies(out, mdp="cutoff_buck")

    @pytest.mark.skip("Broken, unclear if cases like these are worth supporting")
    def test_nonconsecutive_isomorphic_molecules(self, sage_unconstrained):
        molecules = [Molecule.from_smiles(smiles) for smiles in ["CC", "CCO", "CC"]]

        for index, molecule in enumerate(molecules):
            molecule.generate_conformers(n_conformers=1)
            molecule.conformers[0] += unit.Quantity(3 * [5 * index], unit.angstrom)

        topology = Topology.from_molecules(molecules)
        topology.box_vectors = unit.Quantity([4, 4, 4], unit.nanometer)

        out = Interchange.from_smirnoff(sage_unconstrained, topology)

        get_gromacs_energies(out).compare(
            get_openmm_energies(out, combine_nonbonded_forces=True),
            {"Nonbonded": 0.5 * unit.kilojoule_per_mole},
        )


class TestGROMACSMetadata(_BaseTest):
    def test_atom_names_pdb(self):
        peptide = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_ALA_ALA.pdb"),
        )
        ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

        Interchange.from_smirnoff(ff14sb, peptide.to_topology()).to_gro(
            "atom_names.gro",
        )
        Interchange.from_smirnoff(ff14sb, peptide.to_topology()).to_top(
            "atom_names.top",
        )

        pdb_object = app.PDBFile(get_data_file_path("proteins/MainChain_ALA_ALA.pdb"))
        openmm_object = app.GromacsTopFile("atom_names.top")

        pdb_atom_names = [atom.name for atom in pdb_object.topology.atoms()]

        openmm_atom_names = [atom.name for atom in openmm_object.topology.atoms()]

        assert openmm_atom_names == pdb_atom_names


@needs_gmx
@pytest.mark.skip("Needs rewrite")
class TestGROMACSVirtualSites(_BaseTest):
    @pytest.fixture()
    def sigma_hole_type(self, sage):
        """A handler with a bond charge virtual site on a C-Cl bond."""
        return VirtualSiteHandler.VirtualSiteBondChargeType(
            name="EP",
            smirks="[#6:1]-[#17:2]",
            distance=1.4 * unit.angstrom,
            type="BondCharge",
            match="once",
            charge_increment1=0.1 * unit.elementary_charge,
            charge_increment2=0.2 * unit.elementary_charge,
        )

    @pytest.fixture()
    def sage_with_monovalent_lone_pair(self, sage):
        """Fixture that loads an SMIRNOFF XML for argon"""
        virtual_site_handler = VirtualSiteHandler(version=0.3)

        carbonyl_type = VirtualSiteHandler.VirtualSiteType(
            name="EP",
            smirks="[O:1]=[C:2]-[*:3]",
            distance=0.3 * unit.angstrom,
            type="MonovalentLonePair",
            match="all_permutations",
            outOfPlaneAngle=0.0 * unit.degree,
            inPlaneAngle=120.0 * unit.degree,
            charge_increment1=0.05 * unit.elementary_charge,
            charge_increment2=0.1 * unit.elementary_charge,
            charge_increment3=0.15 * unit.elementary_charge,
        )

        virtual_site_handler.add_parameter(parameter=carbonyl_type)
        sage.register_parameter_handler(virtual_site_handler)

        return sage

    @pytest.mark.xfail()
    @skip_if_missing("parmed")
    def test_sigma_hole_example(self, sage_with_sigma_hole):
        """Test that a single-molecule sigma hole example runs"""
        mol = Molecule.from_smiles("CCl")
        mol.name = "Chloromethane"
        mol.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(
            force_field=sage_with_sigma_hole,
            topology=mol.to_topology(),
        )
        out.box = [4, 4, 4]
        out.positions = mol.conformers[0]

        # TODO: Sanity-check reported energies
        get_gromacs_energies(out)

        out.to_top("sigma.top")
        gmx_top = parmed.load_file("sigma.top")

        assert abs(numpy.sum([p.charge for p in gmx_top.atoms])) < 1e-3

    def test_carbonyl_example(self, sage_with_monovalent_lone_pair):
        """Test that a single-molecule DivalentLonePair example runs"""
        mol = Molecule.from_smiles("C=O")
        mol.name = "Carbon_monoxide"
        mol.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(
            force_field=sage_with_monovalent_lone_pair,
            topology=mol.to_topology(),
        )
        out.box = [4, 4, 4]
        out.positions = mol.conformers[0]

        with pytest.raises(
            VirtualSiteTypeNotImplementedError,
            match="MonovalentLonePair not implemented.",
        ):
            # TODO: Sanity-check reported energies
            get_gromacs_energies(out)
