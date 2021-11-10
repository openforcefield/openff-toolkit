#!/usr/bin/env python
# coding: utf-8

#  <div class="alert alert-block alert-info">
#     <b>Note:</b>
#     Interchange is in the process of replacing ParmEd in many workflows, but it still in an alpha testing phase. Our internal tests indicate it is reliable for many small-molecule systems, but it is not yet reliable for complex, multi-component systems and there are likely still rough edges throughout. Feedback is welcome on the <a href=https://github.com/openforcefield/openff-interchange/issues/>Interchange issue tracker.</a></div>

# ## Using OpenFF force fields in Amber and GROMACS
#
# The Open Forcefield Toolkit can create parametrized `openmm.System` objects that can be natively simulated with OpenMM. This example shows the Interchange project can enable parallel workflows using Amber and GROMACS.

# ### Preparing an OpenFF Topology
#
# We start by loading a PDB file containing one copy of ethanol and cyclohexane. Our goal is to create an OpenFF `Topology` object describing this system that we can parametrize with the SMIRNOFF-format "Sage" force field.
#
# The two `Molecule` objects created from the SMILES strings can contain information such as partial charges and stereochemistry that is not included in an OpenMM topology. In this example, partial charges are not explicitly given, and `ForceField` will assign AM1/BCC charges as specified by the "Sage" force field. Note that the OpenFF Toolkit produces partial charges that do not depend on the input conformation of parameterized molecules. See the [FAQ](https://open-forcefield-toolkit.readthedocs.io/en/latest/faq.html#the-partial-charges-generated-by-the-toolkit-don-t-seem-to-depend-on-the-molecule-s-conformation-is-this-a-bug) for more information.

# In[ ]:


try:
    from openmm import app
except ImportError:
    from simtk.openmm import app

from openff.toolkit.utils import get_data_file_path
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField


# In[ ]:


ethanol = Molecule.from_smiles("CCO")
cyclohexane = Molecule.from_smiles("C1CCCCC1")

# Load the PDB file using OpenMM and save the OpenMM Topology
pdbfile = app.PDBFile(
    get_data_file_path("systems/test_systems/1_cyclohexane_1_ethanol.pdb")
)
omm_topology = pdbfile.topology

# Create the OpenFF Topology.
off_topology = Topology.from_openmm(
    omm_topology, unique_molecules=[ethanol, cyclohexane]
)
off_topology


# ### Preparing an OpenFF ForceField
#
# Once the `ForceField` class is imported, the only decision to make is which force field to use. An exhaustive list of force fields released by the Open Force Field Initiative can be found [here](from openff.toolkit.typing.engines.smirnoff import ForceField
# ).
#
# Here we will use force field from the "Sage" line.

# In[ ]:


forcefield = ForceField("openff-2.0.0.offxml")
forcefield


# ### Preparing an OpenMM System
#
# Once a force field and topology have been loaded, an `openmm.System` can be generated natively with the OpenFF Toolkit.

# In[ ]:


omm_system = forcefield.create_openmm_system(off_topology)
omm_system


# ### Preparing an Interchange object
#
# To exports to engines other than OpenMM, we will make use of the [Interchange](https://openff-interchange.readthedocs.io/) project. There is a high-level `Interchange.from_smirnoff` function that consumes OpenFF Toolkit and ForceField objects and produces an `Interchange` object which can then be exported to formats understood by other molecular simulation engines. This extra step is needed to provide a clean interface between _applied_ parameters and engines. Note also that this step does not require an OpenMM System to be generated; `ForceField.create_openmm_system` does not need to be called to use Amber and GROMACS.

# In[ ]:


from openff.interchange.components.interchange import Interchange

interchange = Interchange.from_smirnoff(
    force_field=forcefield,
    topology=off_topology,
)
interchange.positions = pdbfile.positions
interchange


# ### Exporting to Amber and GROMACS files
#
# Once an `Interchange` object has been constructed, its API can be used to export to files understood by GROMACS, Amber, and more.

# In[ ]:


# Export AMBER files.
interchange.to_prmtop("system.prmtop")
interchange.to_inpcrd("system.inpcrd")

# Export GROMACS files.
interchange.to_top("system.top")
interchange.to_gro("system.gro")


# ### Validating the conversion to Amber files
#
# The Interchange project includes functions that take in an `Interchange` object and call out to simulation engines to run single-point energy calculations (with no minimization or dynamics) for the purpose of validating the export layer with each engine. Under the hood, each of these functions calls API points like those used above while converting to files understood by each engine. These rely on having each engine installed and accessible in the current `$PATH`.

# In[ ]:


from openff.interchange.drivers import get_amber_energies, get_openmm_energies


# In[ ]:


openmm_energies = get_openmm_energies(interchange)
openmm_energies.energies


# In[ ]:


# get_ipython().system('cat system.inpcrd')
# get_ipython().system('cat system.prmtop')
# get_ipython().system('cat system.top')
# get_ipython().system('cat system.gro')
amber_energies = get_amber_energies(interchange)
amber_energies.energies


# ### Appendix: Validating the conversion to GROMACS and LAMMPS files
#
# If GROMACS and/or LAMMPS are installed on your machine, the same comparisons can also take place with those engines. They are available via `conda` by running a command like:
#
# ```conda install gromacs lammps -c conda-forge -c bioconda```

# In[ ]:


from distutils.spawn import find_executable
from pprint import pprint

from openff.interchange.drivers import get_gromacs_energies, get_lammps_energies


# In[ ]:


if find_executable("lmp_serial"):
    pprint(get_lammps_energies(interchange).energies)


# In[ ]:


if find_executable("gmx"):
    pprint(get_gromacs_energies(interchange).energies)


# Finally, there is a helper function `get_summary_data` that will attempt to run drivers of each engine. A summary reported is prepared as a Pandas `DataFrame`.

# In[ ]:


from openff.interchange.drivers.all import get_summary_data

get_summary_data(interchange)
