from functools import partial
from openforcefield.typing.chemistry import *
from openforcefield.utils import get_data_file_path
from unittest import TestCase
import pytest
from openforcefield.utils.toolkits import OPENEYE_AVAILABLE

# TODO: Bring these back online once OpenEye dependence is resolved
@pytest.mark.skip
class TestChemicalEnvironments(TestCase):
    def test_createEnvironments(self):
        """
        Test all types of ChemicalEnvironment objects with defined atoms and bonds
        Each will be tetrahedral carbons connected by ring single bonds
        """
        carbon = [['#6'], ['X4']]
        singleBond = [['-'], ['@']]
        atom = AtomChemicalEnvironment('[#6X4:1]','CT')
        bond = BondChemicalEnvironment('[#6X4:1]-[#6X4:2]', 'CT-CT')
        angle = AngleChemicalEnvironment('[#6X4:1]-[#6X4:2]-[#6X4:3]', 'CT-CT-CT')
        torsion = TorsionChemicalEnvironment('[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]', 'CT-CT-CT-CT')
        improper = ImproperChemicalEnvironment('[#6X4:1]-[#6X4:2](-[#6X4:3])-[#6X4:4]', 'CT-CT(-CT)-CT')

    # TODO: Can we remove explicit OE dependence from this?
    @pytest.mark.skip
    def test_complicatedTorsion(self):
        """
        Test ChemicalEnvironment objects with complicated torsion
        test methods that add atoms, remove atoms
        add ORtypes and ANDtypes to existing atoms

        This is the SMIRK for the final torsion
        "[*:1] - [#6:2](=[#8,#7;H0]) - [#6:3](-[#7X3,#8X2;+0]-[#1])(-[#1]) - [*:4]"
        """
        #from openeye.oechem import *
        torsion_smirks = "[*:1]-[#6:2]-[#6:3]-[*:4]"
        torsion = TorsionChemicalEnvironment(torsion_smirks)
        # save atoms (use selectAtom)
        atom1 = torsion.selectAtom(1)
        atom2 = torsion.selectAtom(2)
        atom3 = torsion.selectAtom(3)

        # Add atoms with names so I can try to remove them
        atom2alpha = torsion.addAtom(atom2, [('=',[])], None, [('#8',[]),('#7',[])], ['H0'])
        atom3alpha1 = torsion.addAtom(atom3)
        atom3beta1 = torsion.addAtom(atom3alpha1, [('-',[])], None, [('#1',[])])
        atom3alpha2 = torsion.addAtom(atom3, [('-',[])], None, [('#1',[])])

        # Get bond for atom3 and alpha and add ANDtype
        bond = torsion.getBond(atom3, atom3alpha1)
        self.assertIsNotNone(bond, "could not find bond between atom3 and it's alpha atom")
        bond.addORtype('-', [])

        # Add ORtypes and ANDtypes to atom3 alpha atom
        atom3alpha1.addORtype('#7', ['X3'])
        atom3alpha1.addORtype('#8', ['X2'])
        atom3alpha1.addANDtype('+0')

        # Call getAtoms and getBonds just to make sure they work
        torsion.getAtoms()
        torsion.getBonds()

        # get smarts and smirks for the large torsion
        smarts = torsion.asAtomtypeSMARTS()
        smirks = torsion.asSMIRKS()
        qmol = OEQMol()
        if not OEParseSmarts(qmol, smirks):
            raise Exception("could not parse created SMIRKS %s" % smirks)

        # Try removing atoms
        # if it was labeled:
        removed = torsion.removeAtom(atom1)
        self.assertFalse(removed,"removed labeled atom (%s) in %s" % (atom1.asSMIRKS(), torsion.asSMIRKS()))
        removed = torsion.removeAtom(atom3alpha1)
        self.assertFalse(removed, "removed bridging atom (%s) in %s" % (atom3alpha1.asSMIRKS(), torsion.asSMIRKS()))
        removed = torsion.removeAtom(atom3beta1)
        self.assertTrue(removed, "failed to remove atom (%s) in %s which is removeable" % (atom3beta1.asSMIRKS(), torsion.asSMIRKS()))

    def test_parseSMIRKS(self):
        """
        Test creating environments with SMIRKS
        """
        replacements = [ ('ewg1', '[#7!-1,#8!-1,#16!-1,#9,#17,#35,#53]'),
                ('ewg2', '[#7!-1,#8,#16]') ]

        smirksList = [ ["[#6](-[#1])-[#8]", None, ChemicalEnvironment],
                ["[#6&X4&H0:1](-[#1])-[#6&X4]", 'Atom', AtomChemicalEnvironment],
                [ "[#6&X4&H0:1](-[#1])-[#6&X4:2]", 'Bond', BondChemicalEnvironment],
                [ "[*:1]-[*:2](-[#6&X4])-[*:3]", 'Angle', AngleChemicalEnvironment],
                [ "[#6&X4&H0:1](-[#1])-[#6&X4:2]-[#6&X4&H0:3](-[#1])-[#6&X4:4]", 'ProperTorsion', TorsionChemicalEnvironment],
                [ "[#1:1]-[#6&X4:2](-[#8:3])-[#1:4]", 'ImproperTorsion', ImproperChemicalEnvironment],
                [ "[#1:1]-[#6&X4:2](-[#8:3])-[*:4](-[#6&H1])-[#8:5]", None, ChemicalEnvironment],
                [ "[#6$(*~[#6]=[#8])$(*-,=[$ewg2,#7])]", None, ChemicalEnvironment],
                [ "CCC", None, ChemicalEnvironment],
                [ "[#6:1]1(-;!@[#1,#6])=;@[#6]-;@[#6]1", 'Atom', ChemicalEnvironment],
                [ "C(O-[#7,#8])CC=[*]", None, ChemicalEnvironment],
                [ "[#6$([#6X4](~[$ewg1])(~[#8]~[#1])):1]-[#6X2H2;+0:2]-,=,:;!@;!#[$ewg2:3]-[#4:4]", 'ProperTorsion', TorsionChemicalEnvironment],
                [ "[#6$([#6X4](~[$ewg1])(~[#8]~[#1])):1]1=CCCC1", 'Atom', AtomChemicalEnvironment],
                [ "[*:1]-[#7X3:2](-[#6a$(*1ccc(-[#8-1X1])cc1):3])-[*:4]", 'ImproperTorsion', ImproperChemicalEnvironment],
                [ "[#6X4:1]1~[*:2]~[*$(*~[#1]):3]1", 'Angle', AngleChemicalEnvironment],
                [ "[$([#7]1~[#6]-CC1)]", None, ChemicalEnvironment],
                [ "[$(c1ccccc1)]", None, ChemicalEnvironment],
                ]

        for toolkit in ['openeye', 'rdkit']:
            for [smirks, checkType, chemEnv] in smirksList:
                env = chemEnv(smirks = smirks, replacements = replacements, toolkit=toolkit)
                Type = env.getType()
                self.assertEqual(Type,checkType,
                        "SMIRKS (%s) clasified as %s instead of %s using %s toolkit" % (smirks, Type, checkType, toolkit))

    def test_environment_functions(self):
        """
        Test all minor functions in environments
        """
        # example 3 indexed, 1 unindexed, 1 alpha, 0 beta
        angle_smirks = "[#6X3;R1:1]=,:;@[#6X3;R1;a:2](-,:;@[#7])-;!@[#8X2H1;!R:3]"
        angle = AngleChemicalEnvironment(angle_smirks)

        # check selectAtom and selectBond
        testSelect = [(4,True), ('Beta', True), (2, False), ('Indexed', False),
                ('Unindexed', False), ('Alpha', False), (None, False)]
        for (descriptor, isNone) in testSelect:
            atom = angle.selectAtom(descriptor)
            bond = angle.selectBond(descriptor)
            if isNone:
                self.assertIsNone(atom, "Found atom with descriptor %s in angle %s" % (descriptor, angle_smirks))
                self.assertIsNone(bond, "Found bond with descriptor %s in angle %s" % (descriptor, angle_smirks))

            else:
                self.assertIsNotNone(atom, "Could not find a %s atom in angle %s" % (descriptor, angle_smirks))
                self.assertIsNotNone(bond, "Could not find a %s bond in angle %s" % (descriptor, angle_smirks))

        # Check function getComponentList which uses all functions of the form
        # get_Atoms and get_Bonds
        compOptions = {'atom': [ ('all', 4), ('Indexed', 3), ('Unindexed', 1), ('Alpha', 1), ('Beta', 0)], 'bond': [ ('all', 3), ('Indexed', 2), ('Unindexed', 1), ('Alpha', 1), ('Beta', 0)]}

        for component, desList in compOptions.items():
            for (descriptor, amount) in desList:
                compList = angle.getComponentList(component, descriptor)
                ncomps = len(compList)
                assert ncomps == amount, "Found %i %s %ss instead of %i in angle %s" % (ncomps, descriptor, component, amount, angle_smirks)

        # Check is__ descriptors
        atom2 = angle.selectAtom(2)
        bond1 = angle.selectBond(1)
        alpha_atom = angle.selectAtom('Alpha')
        beta_atom = angle.addAtom(alpha_atom)
        alpha_bond = angle.getBond(atom2, alpha_atom)
        beta_bond = angle.getBond(alpha_atom, beta_atom)

        # list of lists:
        # [ [[components], [(method, expected)]], [...]]
        checkIsMethods = [
                [[atom2,bond1], [(angle.isAlpha, False), (angle.isBeta, False),
                    (angle.isIndexed, True), (angle.isUnindexed, False)]],
                [[alpha_atom, alpha_bond], [(angle.isAlpha, True), (angle.isIndexed, False),
                    (angle.isUnindexed, True)]],
                [[beta_atom, beta_bond], [(angle.isBeta, True)]]]

        for compSet, methodList in checkIsMethods:
            for comp in compSet:
                for (method, expected) in methodList:
                    # same message
                    classify = method(comp)
                    msg = "%s wrongly returned %s for %s in angle %s" % (
                            method, classify, comp.asSMIRKS(), angle.asSMIRKS())
                    if expected:
                        self.assertTrue(classify, msg)
                    else:
                        self.assertFalse(classify, msg)

        # Check getBond when atoms aren't bonded
        atom1 = angle.selectAtom(1)
        beta_to_atom1 = angle.getBond(beta_atom, atom1)
        self.assertIsNone( beta_to_atom1, "Incorrect bond was returned connecting the beta atom (%s) and atom 1 (%s) in angle %s" % (beta_atom.asSMIRKS(), atom1.asSMIRKS(), angle_smirks))

        # Check valence: should be 3 for atom2
        val = angle.getValence(atom2)
        self.assertEqual(val, 3, "Atom 2 in angle %s should have a valence of 3, but a valence of %i was found" % (angle_smirks, val))

        # Check bond order
        # For bond1 =,:;@ it should be 1.5 because order returns lowest possible
        order = bond1.getOrder()
        self.assertEqual(order, 1.5, "Bond 1 (%s) should have a minimum order of 1.5, but %.1f was assigned" % (bond1.asSMIRKS(), order))

        # For atom
        order = angle.getBondOrder(atom2)
        self.assertEqual(order, 3.5, "Atom 2 in angle (%s) should have a minimum bond order of 3.5, but %.1f was assigned" % (angle_smirks, order))

    def test_creating_wrong_environments(self):
        """
        Test exceptions for making environments with the wrong smirks
        """

        wrongDict = {'[*]': [AtomChemicalEnvironment, BondChemicalEnvironment, AngleChemicalEnvironment, TorsionChemicalEnvironment, ImproperChemicalEnvironment],
                "[*:1]": [BondChemicalEnvironment, AngleChemicalEnvironment, TorsionChemicalEnvironment, ImproperChemicalEnvironment],
                "[*:1]~[*:2]": [AtomChemicalEnvironment, AngleChemicalEnvironment, TorsionChemicalEnvironment, ImproperChemicalEnvironment],
                "[*:3]~[*:2]~[*:1]":[AtomChemicalEnvironment, BondChemicalEnvironment, TorsionChemicalEnvironment, ImproperChemicalEnvironment],
                "[*:1]~[*:2]~[*:3]~[*:4]": [AtomChemicalEnvironment, BondChemicalEnvironment, AngleChemicalEnvironment, ImproperChemicalEnvironment],
                "[*:1]~[*:2](~[*:3])~[*:4]": [AtomChemicalEnvironment, BondChemicalEnvironment, AngleChemicalEnvironment, TorsionChemicalEnvironment],
                "[*:1]~[*:2]~[*:3]~[*:4]~[*:5]": [AtomChemicalEnvironment, BondChemicalEnvironment, AngleChemicalEnvironment, TorsionChemicalEnvironment, ImproperChemicalEnvironment] }

        for smirks, Environments in wrongDict.items():
            for Environment in Environments:
                msg = "SMIRKS (%s) is in appropriate input for %s, but no error was raised" % (smirks, str(Environment))
                with self.assertRaises(SMIRKSMismatchError, msg = msg):
                    env = Environment(smirks)

    def test_wrong_smirks_error(self):
        """
        Check that an imparseable SMIRKS raises errors
        """
        smirks = "[*;X:1]"
        msg = "SMIRKS (%s) should not be parseable, but an environment was successfully created"
        with self.assertRaises(SMIRKSParsingError, msg = msg):
            env = ChemicalEnvironment(smirks)

    def test_embedded_atoms_smirks(self):
        """
        Check embedded atom parsing works
        """
        smirks = "[#1$(*-[#6](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]):1]~[$([#1]~[#6])]"
        env = ChemicalEnvironment(smirks)
