#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
forcefield.py

OpenMM ForceField replacement using SMIRKS-based matching.

AUTHORS

John D. Chodera <john.chodera@choderalab.org>

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

import xml.etree.ElementTree as etree

from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye import oechem

import time

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

def _convertParameterToNumber(param):
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)

#=============================================================================================
# FORCEFIELD
#=============================================================================================

class ForceField(object):
    """A ForceField constructs OpenMM System objects based on a Topology.
    """

    def __init__(self, *files):
        """Load one or more XML parameter definition files and create a SMARTY ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the SMARTY force field.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.

        """
        self._forces = []
        self.loadFile(files)

    def loadFile(self, files):
        """Load a SMARTY XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing SMARTY force field definitions.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.
        """

        # Ensure that we are working with a tuple of files.
        if not isinstance(files, tuple):
            files = (files,)

        # Load in all XML trees.
        trees = list()
        for file in files:
            try:
                # this handles either filenames or open file-like objects
                tree = etree.parse(file)
            except IOError:
                tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', file))
            except Exception as e:
                # Fail with an error message about which file could not be read.
                # TODO: Also handle case where fallback to 'data' directory encounters problems,
                # but this is much less worrisome because we control those files.
                msg  = str(e) + '\n'
                if hasattr(file, 'name'):
                    filename = file.name
                else:
                    filename = str(file)
                msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
                raise Exception(msg)

            trees.append(tree)

        # Load the atom masses.
        for tree in trees:
            if tree.getroot().find('AtomTypes') is not None:
                for type in tree.getroot().find('AtomTypes').findall('Type'):
                    self.registerAtomType(type.attrib)

        # Load force definitions
        for tree in trees:
            for child in tree.getroot():
                if child.tag in parsers:
                    parsers[child.tag](child, self)

    @classmethod
    def generateOEMolFromTopologyResidue(cls, residue, geometry=False, tripos_atom_names=False):
        """
        Generate an OpenEye OEMol molecule from an OpenMM Topology Residue.

        Parameters
        ----------
        residue : simtk.openmm.app.topology.Residue
            The topology Residue from which an OEMol is to be created.
            An Exception will be thrown if this residue has external bonds.
        geometry : bool, optional, default=False
            If True, will generate a single configuration with OEOmega.
            Note that stereochemistry will be *random*.
        tripos_atom_names : bool, optional, default=False
            If True, will generate and assign Tripos atom names.

        Returns
        -------
        molecule : openeye.oechem.OEMol
            The OEMol molecule corresponding to the topology.
            Atom order will be preserved and bond orders assigned.

        The Antechamber `bondtype` program will be used to assign bond orders, and these
        will be converted back into OEMol bond type assignments.

        Note that there is no way to preserve stereochemistry since `Residue` does
        not note stereochemistry in any way.

        """
        # Raise an Exception if this residue has external bonds.
        if len(list(residue.external_bonds())) > 0:
            raise Exception("Cannot generate an OEMol from residue '%s' because it has external bonds." % residue.name)

        # Create OEMol where all atoms have bond order 1.
        molecule = oechem.OEMol()
        molecule.SetTitle(residue.name) # name molecule after first residue
        for atom in residue.atoms():
            oeatom = molecule.NewAtom(atom.element.atomic_number)
            oeatom.SetName(atom.name)
            oeatom.AddData("topology_index", atom.index)
        oeatoms = { oeatom.GetName() : oeatom for oeatom in molecule.GetAtoms() }
        for (atom1, atom2) in residue.bonds():
            order = 1
            molecule.NewBond(oeatoms[atom1.name], oeatoms[atom2.name], order)

        # Write out a mol2 file without altering molecule.
        import tempfile
        tmpdir = tempfile.mkdtemp()
        mol2_input_filename = os.path.join(tmpdir,'molecule-before-bond-perception.mol2')
        ac_output_filename = os.path.join(tmpdir,'molecule-after-bond-perception.ac')
        ofs = oechem.oemolostream(mol2_input_filename)
        m2h = True
        substruct = False
        oechem.OEWriteMol2File(ofs, molecule, m2h, substruct)
        ofs.close()
        # Run Antechamber bondtype
        import subprocess
        #command = 'bondtype -i %s -o %s -f mol2 -j full' % (mol2_input_filename, ac_output_filename)
        command = 'antechamber -i %s -fi mol2 -o %s -fo ac -j 2' % (mol2_input_filename, ac_output_filename)
        [status, output] = getstatusoutput(command)

        # Define mapping from GAFF bond orders to OpenEye bond orders.
        order_map = { 1 : 1, 2 : 2, 3: 3, 7 : 1, 8 : 2, 9 : 5, 10 : 5 }
        # Read bonds.
        infile = open(ac_output_filename)
        lines = infile.readlines()
        infile.close()
        antechamber_bond_types = list()
        for line in lines:
            elements = line.split()
            if elements[0] == 'BOND':
                antechamber_bond_types.append(int(elements[4]))
        oechem.OEClearAromaticFlags(molecule)
        for (bond, antechamber_bond_type) in zip(molecule.GetBonds(), antechamber_bond_types):
            #bond.SetOrder(order_map[antechamber_bond_type])
            bond.SetIntType(order_map[antechamber_bond_type])
        oechem.OEFindRingAtomsAndBonds(molecule)
        oechem.OEKekulize(molecule)
        oechem.OEAssignFormalCharges(molecule)
        oechem.OEAssignAromaticFlags(molecule, oechem.OEAroModelOpenEye)

        # Clean up.
        os.unlink(mol2_input_filename)
        os.unlink(ac_output_filename)
        os.rmdir(tmpdir)

        # Generate Tripos atom names if requested.
        if tripos_atom_names:
            oechem.OETriposAtomNames(molecule)

        # Assign geometry
        if geometry:
            from openeye import oeomega
            omega = oeomega.OEOmega()
            omega.SetMaxConfs(1)
            omega.SetIncludeInput(False)
            omega.SetStrictStereo(False)
            omega(molecule)

        return molecule

    def createSystem(self, topology, molecules=nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, removeCMMotion=True, hydrogenMass=None, residueTemplates=dict(), **args):
        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        residueTemplates : dict=dict()
           Key: Topology Residue object
           Value: string, name of _TemplateData residue template object to use for
                  (Key) residue
           This allows user to specify which template to apply to particular Residues
           in the event that multiple matching templates are available (e.g Fe2+ and Fe3+
           templates in the ForceField for a monoatomic iron ion in the topology).
        args
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        # TODO: Create an OEMol object from each molecule in the Topology.
        molecules = ForceField._extractMolecules(topology)

        # Create the System and add atoms
        sys = mm.System()
        for atom in topology.atoms():
            # Add the particle to the OpenMM system.
            sys.addParticle(atom.element.mass) # TODO: Allow option to use a different mass than the integral elemental mass?

        # Adjust hydrogen masses if requested.
        if hydrogenMass is not None:
            if not unit.is_quantity(hydrogenMass):
                hydrogenMass *= unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.
        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # TODO: Convert requested bonds and angles to use constraints

        # Add forces to the System
        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing
        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(sys, data, args)

        return sys

#=============================================================================================
# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.
#=============================================================================================

## @private
class HarmonicBondGenerator(object):
    """A HarmonicBondGenerator constructs a HarmonicBondForce."""

    class Bond(object):
        def __init__(self, node):
            self.smirks = node.attrib['smirks']
            self.length = _convertParameterToNumber(node.attrib['length'])
            self.k = _convertParameterToNumber(node.attrib['k'])

            # Check SMIRKS for validity
            qmol = oechem.OEQMol()
            if not oechem.OEParseSmarts( qmol, bond.smirks ):
                raise Exception('HarmonicBondForce error in ')



    def __init__(self, forcefield):
        self.ff = forcefield
        self._bonds = list()

    def registerBond(self, parameters):
        self._bonds.append(HarmonicBondGenerator.Bond(bond))

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, HarmonicBondGenerator)]
        if len(existing) == 0:
            generator = HarmonicBondGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        for bond in element.findall('Bond'):
            generator.registerBond(bond)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicBondForce]
        if len(existing) == 0:
            force = mm.HarmonicBondForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # Match all bonds of each type.
        for bond in self._bonds:
            qmol = oechem.OEQMol()
            if not oechem.OEParseSmarts( qmol, bond.smirks ):
                raise Exception('HarmonicBondForce error in ')
            ss = oechem.OESubSearch( qmol)
                for match in ss.Match(mol, unique):
        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    bond.length = self.length[i]
                    if bond.isConstrained:
                        data.addConstraint(sys, bond.atom1, bond.atom2, self.length[i])
                    elif self.k[i] != 0:
                        force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
                    break

parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement
