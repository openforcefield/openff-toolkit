import pytest

from openff.toolkit import ForceField


class StaleForceFieldTests:
    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_forcefield_from_file_list(self):
        # These offxml files are located in package data path, which is automatically installed and searched
        file_paths = [smirnoff99Frosst_offxml_file_path, tip3p_offxml_file_path]
        # Create a forcefield from multiple offxml files
        ForceField(file_paths)

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_forcefield_from_file_path_iterator(self):
        # These offxml files are located in package data path, which is automatically installed and searched
        file_paths = [smirnoff99Frosst_offxml_file_path, tip3p_offxml_file_path]
        # A generator should work as well
        ForceField(iter(file_paths))

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_gbsa(self):
        """Test reading of ffxml files with GBSA support."""
        ForceField("test_forcefields/Frosst_AlkEthOH_GBSA.offxml")

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_forcefield_from_url(self):
        urls = [
            "https://raw.githubusercontent.com/openforcefield/openff-toolkit/master/openff/toolkit/data/test_forcefields/test_forcefield.offxml",
            "https://raw.githubusercontent.com/openforcefield/openff-toolkit/master/openff/toolkit/data/test_forcefields/tip3p.offxml",
        ]
        # Test creation with smirnoff99frosst URL
        ForceField(urls[0])

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_forcefield_from_url_list(self):
        urls = [
            "https://raw.githubusercontent.com/openforcefield/openff-toolkit/master/openff/toolkit/data/test_forcefields/test_forcefield.offxml",
            "https://raw.githubusercontent.com/openforcefield/openff-toolkit/master/openff/toolkit/data/test_forcefields/tip3p.offxml",
        ]
        # Test creation with multiple URLs
        ForceField(urls)

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_forcefield_from_url_iterator(self):
        urls = [
            "https://raw.githubusercontent.com/openforcefield/openff-toolkit/master/openff/toolkit/data/test_forcefields/test_forcefield.offxml",
            "https://raw.githubusercontent.com/openforcefield/openff-toolkit/master/openff/toolkit/data/test_forcefields/tip3p.offxml",
        ]
        # A generator should work as well
        ForceField(iter(urls))

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_charge_increment(self):
        """Test parameter assignment using smirnoff99Frosst on laromustine with ChargeIncrementModel."""
        molecules_file_path = get_data_file_path("molecules/laromustine_tripos.mol2")
        molecule = Molecule.from_file(molecules_file_path)
        forcefield = ForceField(
            ["test_forcefields/test_forcefield.offxml", "chargeincrement-test"]
        )
        check_system_creation_from_molecule(forcefield, molecule)
        # TODO: We can't implement a test for chargeincrement yet because we
        #       haven't settled on a SMIRNOFF spec for chargeincrementmodel

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_create_system_molecules_parmatfrosst_gbsa(self):
        """Test creation of a System object from small molecules to test parm@frosst force field with GBSA support."""
        molecules_file_path = get_data_file_path(
            "molecules/AlkEthOH_test_filt1_tripos.mol2"
        )
        check_parameter_assignment(
            offxml_file_path="test_forcefields/Frosst_AlkEthOH_GBSA.offxml",
            molecules_file_path=molecules_file_path,
        )
        # TODO: Figure out if we just want to check that energy is finite (this is what the original test did,
        #       or compare numerically to a reference system.

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    def test_deep_copy(self):
        force_field = ForceField(smirnoff99Frosst_offxml_file_path)
        # Deep copy
        force_field2 = copy.deepcopy(force_field)
        assert_forcefields_equal(
            force_field,
            force_field2,
            "ForceField deep copy does not match original ForceField",
        )

    @pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
    # TODO: This should check the output of forcefield.to_dict
    def test_serialize(self):
        force_field = ForceField(smirnoff99Frosst_offxml_file_path)
        # Serialize/deserialize
        serialized_forcefield = force_field.__getstate__()
        force_field2 = ForceField.__setstate__(serialized_forcefield)
        assert_forcefields_equal(
            force_field,
            force_field2,
            "Deserialized serialized ForceField does not match original ForceField",
        )


@pytest.mark.skip(reason="Needs to be updated for 0.2.0 syntax")
def test_electrostatics_options(self):
    """Test parameter assignment using smirnoff99Frosst on laromustine with various long-range electrostatics options."""
    from functools import partial

    molecules_file_path = get_data_file_path("molecules/laromustine_tripos.mol2")
    molecule = Molecule.from_file(molecules_file_path)
    forcefield = ForceField(
        [smirnoff99Frosst_offxml_file_path, charge_increment_offxml_file_path]
    )
    for method in ["PME", "reaction-field", "Coulomb"]:
        # Change electrostatics method
        forcefield.forces["Electrostatics"].method = method
        f = partial(check_system_creation_from_molecule, forcefield, molecule)
        f.description = "Testing {} parameter assignment using molecule {}".format(
            offxml_file_path, molecule.name
        )
        # yield f
    # TODO: Implement a similar test, where we compare OpenMM energy evals from an
    #       AMBER-parameterized system to OFF-parameterized systems


class StaleIOTests:
    @pytest.mark.skip(reason="Needs to be updated for 1.0.0 syntax")
    def test_to_xml(self):
        forcefield = ForceField(smirnoff99Frosst_offxml_filename)
        # Retrieve XML as a string
        xml = forcefield.to_xml()
        # Restore ForceField from XML
        forcefield2 = ForceField(xml)
        assert_forcefields_equal(
            cls.forcefield,
            forcefield2,
            "ForceField serialized to XML does not match original ForceField",
        )

    # TODO: Remove ForceField from this whole file. All tests should be for converting between hierarchical SMIRNOFF
    #       dicts and XML
    @pytest.mark.skip(reason="Needs to be updated for 1.0.0 syntax")
    def test_save(self):
        """Test writing and reading of SMIRNOFF in XML format."""
        forcefield = ForceField(smirnoff99Frosst_offxml_filename)
        # Write XML to a file
        with TemporaryDirectory() as tmpdir:
            offxml_tmpfile = os.path.join(tmpdir, "forcefield.offxml")
            forcefield.save(offxml_tmpfile)
            forcefield2 = ForceField(offxml_tmpfile)
            assert_forcefields_equal(
                cls.forcefield,
                forcefield2,
                "ForceField written to .offxml does not match original ForceField",
            )
