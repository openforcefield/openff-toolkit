from functools import partial
import smarty
import smarty.environment
from smarty.environment import *
from smarty.utils import get_data_filename
from unittest import TestCase
import openeye.oechem
from openeye.oechem import *

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

    def test_complicatedTorsion(self):
        """
        Test ChemicalEnvironment objects with complicated torsion
        test methods that add atoms, remove atoms
        add ORtypes and ANDtypes to existing atoms

        This is the SMIRK for the final torsion
        "[*:1] - [#6:2](=[#8,#7;H0]) - [#6:3](-[#7X3,#8X2;+0]-[#1])(-[#1]) - [*:4]"
        """
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
        if bond == None:
            # If None, bond wasn't found correctly
            raise Exception("could not find bond between atom3 and it's alpha atom")
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
        if removed:
            raise Exception("removed labeled atom")
        removed = torsion.removeAtom(atom3alpha1)
        if removed:
            raise Exception("removed bridging atom")
        removed = torsion.removeAtom(atom3beta1)
        if not removed:
            raise Exception("did not remove an atom that should be allowed")

    def test_parseSMIRKS(self):
        """
        Test creating environments with SMIRKS
        """
        smirksList = [ ["[#6](-[#1])-[#8]", None],
                ["[#6&X4&H0:1](-[#1])-[#6&X4]", 'VdW'],
                [ "[#6&X4&H0:1](-[#1])-[#6&X4:2]", 'Bond'],
                [ "[*:1]-[*:2](-[#6&X4])-[*:3]", 'Angle'],
                [ "[#6&X4&H0:1](-[#1])-[#6&X4:2]-[#6&X4&H0:3](-[#1])-[#6&X4:4]", 'Torsion'],
                [ "[#1:1]-[#6&X4:2](-[#8:3])-[#1:4]", 'Improper'],
                [ "[#1:1]-[#6&X4:2](-[#8:3])-[*:4](-[#6&H1])-[#8:5]", None] ]

        for [smirks, checkType] in smirksList:
            env = ChemicalEnvironment(smirks)
            Type = env.getType()
            if Type != checkType:
                raise Exception("SMIRKS (%s) clasified as %s instead of %s" % (smirks, Type, checkType))
