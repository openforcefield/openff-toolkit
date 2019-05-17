#!/usr/bin/env python

import time

# Import OpenMM tools
from simtk import openmm, unit
from simtk.openmm.app import *

# Use MDTraj to write simulation trajectories
from mdtraj.reporters import NetCDFReporter

# Import the SMIRNOFF forcefield engine and some useful tools
from openforcefield.typing.engines.smirnoff import ForceField
# LPW: openforcefield's PME is different from openmm's PME
from openforcefield.typing.engines.smirnoff.forcefield import PME
from openforcefield.utils import get_data_file_path

# Import the OpenEye toolkit
from openeye import oechem

pdb_filename = 'ETH-box.pdb'
mol_filename = 'ETH.mol2'

# Define a few simulation parameters
time_step = 2*unit.femtoseconds # simulation timestep
temperature = 300*unit.kelvin # simulation temperature
friction = 1/unit.picosecond # collision rate
num_steps = 1000000 # number of steps to run
trj_freq = 1000 # number of steps per written trajectory frame
data_freq = 1000 # number of steps per written simulation statistics

# Load molecule and create pdb object
pdb = PDBFile(pdb_filename)

# Load a SMIRNOFF forcefield
forcefield = ForceField(get_data_file_path('test_forcefields/Frosst_AlkEthOH_parmAtFrosst.offxml'))

# Load molecule using OpenEye tools
mol = oechem.OEGraphMol()
ifs = oechem.oemolistream(mol_filename)
# LPW: I don't understand the meaning of these lines.
# flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
# ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol)
oechem.OETriposAtomNames(mol)

# Create the OpenMM system
system = forcefield.createSystem(pdb.topology, [mol], nonbondedMethod=PME, nonbondedCutoff=1.0*unit.nanometers, rigidWater=True)

# Set up an OpenMM simulation
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
platform = openmm.Platform.getPlatformByName('CUDA') 
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)
netcdf_reporter = NetCDFReporter('water_traj.nc', trj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(StateDataReporter('water_data.csv', data_freq, step=True, potentialEnergy=True, temperature= True, density=True))

print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy()
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Run the simulation
print("Starting simulation")
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
netcdf_reporter.close()
print("Done!")
