#!/usr/bin/env python

from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField

# Create a propane molecule and its topology.
molecule = Molecule.from_smiles('CCC')
topology = Topology.from_molecules([molecule])

# Load the Frosst AlkEthOH force field.
ff = ForceField('Frosst_AlkEthOH.offxml')

labels = ff.label_molecules(topology, verbose=True)
for mol_entry in range(len(labels)):
    for force in labels[mol_entry].keys():
        print("\n%s:" % force)
        for (atom_indices, parameter) in labels[mol_entry][force].items():
            atomstr=''
            for idx in atom_indices:
                atomstr += '%6s' % idx
            print("%s : %s \t smirks %s" % (atomstr, parameter.id, parameter.smirks) )
