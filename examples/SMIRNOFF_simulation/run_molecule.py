#!/bin/env python

import time
import numpy as np

# Import OpenMM tools
from simtk import openmm, unit
from simtk.openmm import app, Platform

# Use MDTraj to write simulation trajectories
from mdtraj.reporters import NetCDFReporter

# Import the SMIRNOFF forcefield engine and some useful tools
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol

# Define what molecule to work on, and a few simulation parameters
mol_filename = 'AlkEthOH_r51.mol2'
offxml_filename = 'Frosst_AlkEthOH_parmAtFrosst.offxml'  # SMIRNOFF forcefield
time_step = 2*unit.femtoseconds # simulation timestep
temperature = 300*unit.kelvin # simulation temperature
friction = 1/unit.picosecond # collision rate
num_steps = 100000 # number of steps to run
trj_freq = 1000 # number of steps per written trajectory frame
data_freq = 1000 # number of steps per written simulation statistics

# Load molecule
from openforcefield.topology import Molecule
molecule = Molecule.from_file(mol_filename)

# Get positions in OpenMM-compatible format
positions = molecule.positions

# Create an openforcefield Topology
from openforcefield.topology import Topology
topology = Topology.from_molecules(molecule)

# Load a SMIRNOFF small molecule forcefield for alkanes, ethers, and alcohols
forcefield = ForceField(offxml_filename)

# Create the OpenMM system
system = forcefield.create_system(topology)

# Set up an OpenMM simulation
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
platform = openmm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(topology.to_openmm(), system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)
netcdf_reporter = NetCDFReporter('trajectory.nc', trj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(app.StateDataReporter('data.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

# Run the simulation
print("Starting simulation")
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
netcdf_reporter.close()
print("Done!")
