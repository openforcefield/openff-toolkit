from functools import partial
from smarty import environment
import smarty
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
        atom = environment.AtomChemicalEnvironment(carbon)
        bond = environment.BondChemicalEnvironment(carbon, singleBond, carbon)
        angle = environment.AngleChemicalEnvironment(
                carbon, singleBond, carbon, singleBond, carbon)
        torsion = environment.TorsionChemicalEnvironment(
                carbon, singleBond, carbon, singleBond, carbon, singleBond, carbon)
        improper = environment.ImproperChemicalEnvironment(
                carbon, singleBond, carbon, singleBond, carbon, singleBond, carbon)

    def test_complicatedTorsion(self):
        """
        Test ChemicalEnvironment objects with complicated torsion
        test methods that add atoms, remove atoms
        add ORtypes and ANDtypes to existing atoms

        This is the SMIRK for the final torsion
        "[*:1] - [#6:2](=[#8,#7;H0]) - [#6:3](-[#7X3,#8X2;+0]-[#1])(-[#1]) - [*:4]"
        """
        carbon = [['#6'], None]
        single = [['-'], None]
        torsion = environment.TorsionChemicalEnvironment(Atom2Info = carbon,
                Bond2Info = single, Atom3Info = carbon, Bond3Info = single)

        # save atoms (use selectAtom)
        atom1 = torsion.selectAtom(1)
        atom2 = torsion.selectAtom(2)
        atom3 = torsion.selectAtom(3)

        # Add atoms with names so I can try to remove them
        atom2alpha = torsion.addAtom(atom2, ['='], None, ['#8','#7'], ['H0'])
        atom3alpha1 = torsion.addAtom(atom3)
        atom3beta1 = torsion.addAtom(atom3alpha1, ['-'], None, ['#1'])
        atom3alpha2 = torsion.addAtom(atom3, ['-'], None, ['#1'])

        # Get bond for atom3 and alpha and add ANDtype
        bond = torsion.getBond(atom3, atom3alpha1)
        if bond == None:
            # If None, bond wasn't found correctly
            raise Exception("could not find bond between atom3 and it's alpha atom")
        bond.addORtype('-')

        # Add ORtypes and ANDtypes to atom3 alpha atom
        atom3alpha1.addORtype('#7X3')
        atom3alpha1.addORtype('#8X2')
        atom3alpha1.addANDtype('+0')

        # Call getAtoms and getBonds just to make sure they work
        torsion.getAtoms()
        torsion.getBonds()

        # get smarts and smirks for the large torsion
        smarts = torsion.asAtomtypeSMARTS()
        smirks = torsion.asSMIRKS()
        # TODO: add test that these are relevant

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
            env = environment.ChemicalEnvironment(smirks)
            Type = env.getType()
            if Type != checkType:
                raise Exception("SMIRKS (%s) clasified as %s instead of %s" % (smirks, Type, checkType))
