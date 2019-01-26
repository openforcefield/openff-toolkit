#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Utilities for testing.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import collections
import copy
import functools
import itertools
import os
import pprint
import textwrap

import numpy as np
from simtk import unit, openmm

from openforcefield.utils import get_data_filename


#=============================================================================================
# Shortcut functions to get file paths to test data.
#=============================================================================================

def get_amber_filepath(prefix):
    """Get AMBER prmtop and inpcrd test data filepaths.

    Parameters
    ----------
    prefix : str
        The file name without extension of .prmtop and .inpcrd
        files to retrieve from testdata/systems/amber.

    Returns
    -------
    prmtop_filepath : str
        Absolute path to the AMBER prmtop filepath in testdata/systems/amber
    inpcrd_filepath : str
        Absolute path to the AMBER inpcrd filepath in testdata/systems/amber

    """
    prefix = os.path.join('systems', 'amber', prefix)
    prmtop_filepath = get_data_filename(prefix+'.prmtop')
    inpcrd_filepath = get_data_filename(prefix+'.inpcrd')
    return prmtop_filepath, inpcrd_filepath


def get_packmol_pdbfile(prefix='cyclohexane_ethanol_0.4_0.6'):
    """Get PDB filename for a packmol-generated box

    Parameters
    ----------
    prefix : str, optional, default='cyclohexane_ethanol_0.4_0.6'
        The prefix of .pdb file to retrieve from testdata/systems/packmol_boxes

    Returns
    -------
    pdb_filename : str
        Absolute path to the PDB file
    """
    prefix = os.path.join('systems', 'packmol_boxes', prefix)
    pdb_filename = get_data_filename(prefix+'.pdb')
    return pdb_filename


def get_monomer_mol2file(prefix='ethanol'):
    """Get absolute filepath for a mol2 file denoting a small molecule monomer in testdata

    Parameters
    ----------
    prefix : str, optional, default='ethanol'
        The prefix of .mol2 file to retrieve from systems/monomers/

    Returns
    -------
    mol2_filename : str
        Absolute path to the mol2 file
    """
    prefix = os.path.join('systems', 'monomers', prefix)
    mol2_filename = get_data_filename(prefix+'.pdb')
    return mol2_filename


def get_alkethoh_filepath(alkethoh_name, get_amber=False):
    """Retrieve the mol2, top and crd files of a molecule in the AlkEthOH set.

    Parameters
    ----------
    alkethoh_name : str
        The name of the AlkEthOH molecule (e.g. "AlkEthOH_r0", "AlkEthOH_c1266").
    get_amber : bool, optional
        If True, the paths to the top and crd files are returned.

    Returns
    -------
    molecule_file_paths : str or List[str]
        All the requested paths. If ``get_amber`` is False, only a single string
        pointing to the path of the mol2 file is returned, otherwise this is a
        list ``[mol2_path, top_path, crd_path]``.

    """
    import tarfile

    # Determine if this is a ring or a chain molecule and the subfolder name.
    is_ring = alkethoh_name[9] == 'r'
    alkethoh_subdir_name = 'rings' if is_ring else 'chain'
    alkethoh_subdir_name = 'AlkEthOH_' + alkethoh_subdir_name + '_filt1'

    # Determine which paths have to be returned. Paths are
    # relative to the `data/molecules/` folder. We'll re-use
    # these relative paths to extract files from the tar.gz.
    molecule_file_relative_base_path = os.path.join('AlkEthOH_tripos', alkethoh_subdir_name, alkethoh_name)
    # We always return the mol2 file.
    molecule_relative_file_paths = [molecule_file_relative_base_path + '_tripos.mol2']
    # Check if we need to return also Amber files.
    if get_amber:
        molecule_relative_file_paths.append(molecule_file_relative_base_path + '.top')
        molecule_relative_file_paths.append(molecule_file_relative_base_path + '.crd')

    # Build absolute paths.
    molecules_dir_path = get_data_filename('molecules')
    molecule_file_paths = [os.path.join(molecules_dir_path, p) for p in molecule_relative_file_paths]

    # Check if we need to extract some of the files from the tar archive.
    files_to_extract = set()
    for file_idx, molecule_file_path in enumerate(molecule_file_paths):
        if not os.path.isfile(molecule_file_path):
            files_to_extract.add(molecule_relative_file_paths[file_idx])

    # Extract the files.
    if len(files_to_extract) > 0:
        alkethoh_tar_file_path = os.path.join(molecules_dir_path, 'AlkEthOH_tripos.tar.gz')
        with tarfile.open(alkethoh_tar_file_path, 'r:gz') as tar:
            # Find the files to extract.
            members = [m for m in tar.getmembers() if m.name in files_to_extract]
            tar.extractall(path=molecules_dir_path, members=members)

    # Decide whether to return a single path or a list of paths.
    if len(molecule_file_paths) == 1:
        return molecule_file_paths[0]
    return molecule_file_paths


#=============================================================================================
# Shortcut functions to create System objects from system files.
#=============================================================================================

def create_system_from_amber(prmtop_filepath, inpcrd_filepath, *args, **kwargs):
    """Create an OpenMM System and Topology from the AMBER files.

    Parameters
    ----------
    prmtop_filepath : str
        Path to the topology/parameter file in AMBER prmtop format.
    inpcrd_filepath : str
        Path to the coordinates file in AMBER inpcrd or rst7 format.
    *args
    **kwargs
        Other parameters to pass to ``simtk.openmm.app.AmberPrmtopFile.createSystem``.

    Returns
    -------
    system : simtk.openmm.System
        The OpenMM ``System`` object with the parameters.
    topology : simtk.openmm.app.Topology
        The OpenMM ``Topology`` object loaded from the AMBER files.
    positions : simtk.unit.Quantity
        Initial positions loaded from the inpcrd or restart file.

    """
    prmtop_file = openmm.app.AmberPrmtopFile(prmtop_filepath)
    # AmberInpcrdFile parses also rst7 files.
    inpcrd_file = openmm.app.AmberInpcrdFile(inpcrd_filepath)
    # Create system and update box vectors (if needed)
    system = prmtop_file.createSystem(*args, **kwargs)
    if inpcrd_file.boxVectors is not None:
        system.setDefaultPeriodicBoxVectors(*inpcrd_file.boxVectors)
    # Store numpy positions.
    positions = inpcrd_file.getPositions(asNumpy=True)
    return system, prmtop_file.topology, positions


#=============================================================================================
# Utility functions for energy comparisons.
#=============================================================================================

def quantities_allclose(quantity1, quantity2):
    """Check if two Quantity objects are close.

    If the quantities are arrays, all their elements must be close.
    If the given arguments are unitless, this simply calls numpy.allclose().

    Parameters
    ----------
    quantity1 : simtk.unit.Quantity
        The first unit to compare.
    quantity2 : simtk.unit.Quantity
        The second unit to compare.

    Returns
    -------
    is_close : bool
        True if the two quantities are close, False otherwise.

    """
    if isinstance(quantity1, unit.Quantity):
        # Check that the two Quantities have compatible units.
        if not quantity1.unit.is_compatible(quantity2.unit):
            raise ValueError("The two quantities don't have compatible units: "
                             "{} and {}".format(quantity1.unit, quantity2.unit))
        # Compare the values stripped of the units.
        quantity1 = quantity1.value_in_unit_system(unit.md_unit_system)
        quantity2 = quantity2.value_in_unit_system(unit.md_unit_system)
    return np.allclose(quantity1, quantity2)


def get_context_potential_energy(context, positions, box_vectors=None,
                                 by_force_group=True):
    """Compute the potential energy of a System in a Context object.

    This is simply a shortcut that takes care of setting
    box vectors and positions before computing the energy.

    Parameters
    ----------
    context : simtk.openmm.Context
        The Context object containing the system.
    positions : simtk.unit.Quantity
        A n_atoms x 3 arrays of coordinates with units of length.
    box_vectors : simtk.unit.Quantity, optional
        If specified, this should be a 3 x 3 array. Each row should
        be a box vector.
    by_force_group : bool, optional
        If True, it returns a dictionary mapping force groups to
        their potential energies. Default is True.

    Returns
    -------
    potential_energy : simtk.unit.Quantity or Dict[type, Quantity]
        The potential energy of the system at the given coordinates.
        If ``by_force_group`` is True, then this is a dictionary
        mapping force groups to their potential energy.

    """
    # Box vectors must be updated before positions are set.
    if box_vectors is not None:
        context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions)

    # If we don't have energies by force group,
    # just return the energy of the full System.
    if not by_force_group:
        return context.getState(getEnergy=True).getPotentialEnergy()

    # Determine how many force groups we need to compute.
    force_groups = {f.getForceGroup() for f in context.getSystem().getForces()}

    # Compute potential energies of all groups.
    potential_energies = {}
    for force_group in force_groups:
        state = context.getState(getEnergy=True, groups={force_group})
        potential_energies[force_group] = state.getPotentialEnergy()
    return potential_energies


class FailedEnergyComparisonError(AssertionError):
    """Error raised when the energy comparison between two system fails.

    Attributes
    ----------
    potential_energy1
    potential_energy2

    """
    def __init__(self, err_msg, potential_energy1, potential_energy2):
        super().__init__(err_msg)
        self.potential_energy1 = potential_energy1
        self.potential_energy2 = potential_energy2


def compare_context_energies(context1, context2, *args, **kwargs):
    """Compare energies of two Contexts given the same positions.

    Parameters
    ----------
    context1 : simtk.openmm.Context
        The first Context object to compare containing the system.
    context2 : simtk.openmm.Context
        The second Context object to compare containing the system.
    positions : simtk.unit.Quantity
        A n_atoms x 3 arrays of coordinates with units of length.
    box_vectors : simtk.unit.Quantity, optional
        If specified, this should be a 3 x 3 array. Each row should
        be a box vector.
    by_force_group : bool, optional
        If True, the energies of each force groups are compared.
        Default is True.

    Returns
    -------
    potential_energy1 : simtk.unit.Quantity or Dict[int, Quantity]
        The potential energy of context1 at the given coordinates.
        If ``by_force_group`` is True, then this is a dictionary
        mapping force groups to their potential energy.
    potential_energy2 : simtk.unit.Quantity or Dict[int, Quantity]
        The potential energy of context2 at the given coordinates.
        If ``by_force_group`` is True, then this is a dictionary
        mapping force groups to their potential energy.

    Raises
    ------
    FailedEnergyComparisonError
        If the potential energies of the two contexts are different
        or if the two Systems in the contexts have forces divided into
        different force groups. The potential energies of context1 and
        context2 can be accessed in the Exception through the attributes
        potential_energy1 and potential_energy2 respectively.

    """
    potential_energy1 = get_context_potential_energy(context1, *args, **kwargs)
    potential_energy2 = get_context_potential_energy(context2, *args, **kwargs)

    def raise_assert(assertion, err_msg, format_args):
        """Shortcut to raise a custom error and format the error message."""
        if not assertion:
            raise FailedEnergyComparisonError(
                err_msg.format(*format_args), potential_energy1, potential_energy2
            )

    # If by_force_group is True, then the return value will be a
    # dictionary and we need to compare force group by force group.
    if isinstance(potential_energy1, unit.Quantity):
        raise_assert(
            assertion=quantities_allclose(potential_energy1, potential_energy2),
            err_msg='potential energy 1 {}, potential energy 2: {}',
            format_args=(potential_energy1, potential_energy2)
        )
    else:  # potential_energy1 is a dict.
        # If the two context expose different force groups raise an error.
        raise_assert(
            assertion=set(potential_energy1) == set(potential_energy2),
            err_msg='The two contexts have different force groups: context1 {}, context2 {}',
            format_args=(sorted(potential_energy1), sorted(potential_energy2))
        )

        # Check the potential energies of all force groups.
        for force_group in potential_energy1:
            energy1 = potential_energy1[force_group]
            energy2 = potential_energy2[force_group]
            raise_assert(
                assertion=quantities_allclose(energy1, energy2),
                err_msg='Force group {} do not have the same energies: context1 {}, context2 {}',
                format_args=(force_group, energy1, energy2)
            )

    return potential_energy1, potential_energy2


def compare_system_energies(system1, system2, positions, box_vectors=None,
                            by_force_type=True, modify_system=False):
    """Compare energies of two Systems given the same positions.

    Parameters
    ----------
    system1 : simtk.openmm.System
        The first Context object to compare containing the system.
    system2 : simtk.openmm.System
        The second Context object to compare containing the system.
    positions : simtk.unit.Quantity
        A n_atoms x 3 arrays of coordinates with units of length.
    box_vectors : simtk.unit.Quantity, optional
        If specified, this should be a 3 x 3 array. Each row should
        be a box vector.
    by_force_type : bool, optional
        If True, the contribution to the energy of each force type
        is compared instead of just the overall potential energy.
        Default is True.
    modify_system : bool, optional
        If ``by_force_type`` is True, the function will change the
        force groups of two ``Systems`` to compute the separate
        contributions. Set this to False to work on a copy and avoid
        modifying the original objects. Default is False.

    Returns
    -------
    potential_energy1 : simtk.unit.Quantity or Dict[str, Quantity]
        The potential energy of context1 at the given coordinates.
        If ``by_force_type`` is True, then this is a dictionary
        mapping force names to their potential energy.
    potential_energy2 : simtk.unit.Quantity or Dict[str, Quantity]
        The potential energy of context2 at the given coordinates.
        If ``by_force_type`` is True, then this is a dictionary
        mapping force names to their potential energy.

    Raises
    ------
    FailedEnergyComparisonError
        If the potential energies of the two contexts are different
        or if the two Systems in the contexts have forces divided into
        different force groups. The potential energies of context1 and
        context2 can be accessed in the Exception through the attributes
        potential_energy1 and potential_energy2 respectively.

    """
    if by_force_type:
        if not modify_system:
            system1 = copy.deepcopy(system1)
            system2 = copy.deepcopy(system2)

        # First check that the two systems have the same force types.
        force_names1 = {f.__class__.__name__ for f in system1.getForces()}
        force_names2 = {f.__class__.__name__ for f in system2.getForces()}
        err_msg = 'The two systems have difference force types: system1 {}, system2 {}'
        assert force_names1 == force_names2, err_msg.format(force_names1, force_names2)

        # Create a map from force group to force class and viceversa.
        group_to_force = {i: force_name for i, force_name in enumerate(force_names1)}
        force_to_group = {force_name: i for i, force_name in group_to_force.items()}

        # Assign force groups.
        for system in [system1, system2]:
            for force in system.getForces():
                force.setForceGroup(force_to_group[force.__class__.__name__])

    # Create Contexts and compare the energies.
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    context1 = openmm.Context(system1, integrator)
    context2 = openmm.Context(system2, copy.deepcopy(integrator))

    def map_energies_by_force_type(potential_energy1, potential_energy2):
        """Convert dictionary force_group -> energy to force_type -> energy."""
        potential_energy1 = {group_to_force[group]: energy
                             for group, energy in potential_energy1.items()}
        potential_energy2 = {group_to_force[group]: energy
                             for group, energy in potential_energy2.items()}
        return potential_energy1, potential_energy2

    # Catch the exception and log table force type -> energy.
    try:
        potential_energy1, potential_energy2 = compare_context_energies(
            context1, context2, positions, box_vectors,
            by_force_group=by_force_type
        )
    except FailedEnergyComparisonError as e:
        # We don't need to convert force groups into force types
        if not by_force_type:
            raise
        # Add to the error message the table of energies by force type.
        potential_energy1, potential_energy2 = map_energies_by_force_type(
            e.potential_energy1, e.potential_energy2
        )
        # Pretty-print the dictionaries.
        table = '\n\npotential energy system1:\n' + pprint.pformat(potential_energy1)
        table += '\n\npotential energy system2:\n' + pprint.pformat(potential_energy2)
        raise type(e)(str(e) + table, e.potential_energy1, e.potential_energy2)

    return potential_energy1, potential_energy2


#=============================================================================================
# Utility functions for parameters comparisons.
#=============================================================================================

class _ParametersComparer:
    """This is just a convenience class to compare and print parameters.

    Parameters
    ----------
    **kwargs
        All the parameters to compare.

    Attributes
    ----------
    parameters : dict
        The dictionary of parameters.

    Examples
    --------
    >>> import copy
    >>> from simtk import unit
    >>> par1 = _ParametersComparer(charge=1.0*unit.elementary_charge,
    ...                            rmin_half=1.5*unit.angstrom)
    >>> par2 = copy.deepcopy(par1)
    >>> par3 = _ParametersComparer(charge=-1.0*unit.elementary_charge,
    ...                            rmin_half=0.15*unit.nanometers)
    >>> par1 == par2
    True
    >>> par1 == par3
    False
    >>> str(par1)
    '(charge: 1 e**2, rmin_half: 0.15 nm)'

    """

    def __init__(self, **parameters):
        self.parameters = parameters

    def pretty_format_diff(self, other, new_line=True, indent=True):
        """Return a pretty-formatted string describing the differences between parameters.

        Parameters
        ----------
        other : _ParametersComparer
        new_line : bool, optional
            Separate different parameters with new lines.
        indent : bool, optional
            Indent the formatting.

        Returns
        -------
        diff_str : str

        """
        diff = self._get_diff(other)
        diff_list = []
        for par_name in sorted(diff.keys()):
            par1, par2 = diff[par_name]
            diff_list.append('{par_name}: {par1} != {par2}'.format(
                par_name=par_name, par1=par1, par2=par2))
        separator = '\n' if new_line else ''
        diff_str = separator.join(diff_list)
        if indent:
            diff_str = textwrap.indent(diff_str, prefix='    ')
        return diff_str

    def _get_diff(self, other):
        """Build a 'diff' including only the parameters that are different.
        """
        assert set(self.parameters.keys()) == set(other.parameters.keys())
        diff = {}
        for par_name, par1_value in self.parameters.items():
            par2_value = other.parameters[par_name]
            # Determine whether this parameter is close or not.
            if not quantities_allclose(par1_value, par2_value):
                diff[par_name] = (par1_value, par2_value)
        return diff

    def __eq__(self, other):
        return len(self._get_diff(other)) == 0

    def __str__(self):
        # Reorder the parameters to print in deterministic order.
        parameter_order = sorted(self.parameters)
        # Quantities in a dictionary are normally printed in the
        # format "Quantity(value, unit=unit)" so we make it prettier.
        par_str = ', '.join('{}: {}'.format(p, self.parameters[p]) for p in parameter_order)
        return '(' + par_str + ')'


class _TorsionParametersComparer:
    """A ParametersComparer class for torsions.

    Torsions require a different comparison because multiple set of
    (periodicity, phase, k) parameters can be associated to the same
    four atoms.

    Parameters
    ----------
    *parameters
        Variable length arguments. Each is a _ParametersComparer instance.

    Attributes
    ----------
    parameters : List[_ParametersComparer]
        The list of _ParametersComparer instances.

    Examples
    --------
    >>> from simtk import unit
    >>> par1 = _ParametersComparer(periodicity=2, phase=0.0*unit.degrees)
    >>> par2 = _ParametersComparer(periodicity=4, phase=180.0*unit.degrees)
    >>> torsion_par = _TorsionParametersComparer(par1, par2)

    """

    def __init__(self, *parameters):
        # *args is a tuple. Convert it to a list to make it appendable.
        self.parameters = list(parameters)

    def pretty_format_diff(self, other, new_line=True, indent=True):
        """Return a pretty-formatted string describing the differences between parameters.

        Parameters
        ----------
        other : _TorsionParametersComparer
        new_line : bool, optional
            Separate different parameters with new lines.
        indent : bool, optional
            Indent the formatting.

        Returns
        -------
        diff_str : str

        """
        diff_str = 'Parameters in first system:\n'
        diff_str += self._pretty_format_parameters(self.parameters, new_line=new_line)
        diff_str = '\nParameters in second system:\n'
        diff_str += self._pretty_format_parameters(other.parameters, new_line=new_line)
        if indent:
            diff_str = textwrap.indent(diff_str, '    ')
        return diff_str

    def _pretty_format_parameters(self, parameters, new_line=True, indent=True):
        # Reorder the parameters by periodicity and then phase to print in deterministic order.
        sort_key = lambda x: (x.parameters['periodicity'], x.parameters['phase'])
        parameters = sorted(parameters, key=sort_key)
        # Quantities in a dictionary are normally printed in the
        # format "Quantity(value, unit=unit)" so we make it prettier.
        separator = '\n' if new_line else ', '
        prefix = '['
        if indent and new_line:
            separator += '    '
        elif indent:
            prefix += '    '
        return prefix + separator.join(str(pars) for pars in parameters) + ']'

    def __eq__(self, other):
        # Look for self.parameters that don't match any other.parameters.
        for parameters1 in self.parameters:
            # Find at least 1 set of parameters in the list that are equal.
            if all(parameters1 != parameters2 for parameters2 in other.parameters):
                return False

        # Look for other.parameters that don't match any other self.parameters.
        for parameters2 in other.parameters:
            # Find at least 1 set of parameters in the list that are equal.
            if all(parameters1 != parameters2 for parameters1 in self.parameters):
                return False
        return True

    def __str__(self):
        return self._pretty_format_parameters(self.parameters, new_line=False, indent=False)


class FailedParameterComparisonError(AssertionError):
    """Error raised when the parameter comparison between two forces fails.

    Attributes
    ----------
    different_parameters : Dict[Hashable, Tuple[_ParameterComparer]]
        All the parameters for which a differences was detected.
        different_parameters[atom_indices] is a tuple
        (parameter_force1, parameter_force2).

    """
    def __init__(self, err_msg, different_parameters):
        super().__init__(err_msg)
        self.different_parameters = different_parameters


def _compare_parameters(parameters_force1, parameters_force2, interaction_type,
                        force_name='', systems_labels=None):
    """Compare the parameters of 2 forces and raises an exception if they are different.

    Parameters
    ----------
    parameters_force1 : Dict[Hashable, _ParametersComparer]
        A dictionary associating a _ParametersComparer entry to a key
        that is usually a set of atom indices (e.g., the atom index
        for nonbonded interactions, two atom indices for bonds and
        1-4 exceptions, 3 indices for angles, and 4 for torsions).
    parameters_force2 : Dict[Hashable, _ParametersComparer]
        The parameters of the force to compare.
    interaction_type : str
        A string describing the type of interactions (e.g., "bond",
        "particle exception", "proper torsion"). This is only used to
        improve the error message where differences between the two
        forces are detected.
    force_name : str, optional
        The name of the force to optionally include in the eventual
        error message.
    systems_labels : Tuple[str], optional
        A pair of strings with a meaningful name for the system. If
        specified, this will be included in the error message to
        improve its readability.

    Raises
    ------
    FailedParameterComparisonError
        If there are differences in the parameters of the two forces.
        The exceptions exposes an attribute ``different_parameters``
        with the different parameters.

    """
    diff_msg = ''

    # Handle force and systems labels default arguments.
    if force_name != '':
        force_name += ' '  # Add space after.
    if systems_labels is not None:
        systems_labels = ' for {} and {} respectively'.format(*systems_labels)
    else:
        systems_labels = ''

    # First check the parameters that are unique to only one of the forces.
    unique_keys1 = set(parameters_force1) - set(parameters_force2)
    unique_keys2 = set(parameters_force2) - set(parameters_force1)
    err_msg = '\n\nForce{} has the following unique ' + interaction_type + 's: {}\n'
    if len(unique_keys1) != 0:
        diff_msg += err_msg.format(1, sorted(unique_keys1))
    if len(unique_keys2) != 0:
        diff_msg += err_msg.format(2, sorted(unique_keys2))

    # Create a diff for entries that have same keys but different parameters.
    different_parameters = {}
    for key in set(parameters_force1).intersection(set(parameters_force2)):
        if parameters_force1[key] != parameters_force2[key]:
            different_parameters[key] = (parameters_force1[key], parameters_force2[key])

    # Print error.
    if len(different_parameters) > 0:
        diff_msg += ('\n\nThe following {interaction}s have different parameters '
                     'in the two {force_name}forces{systems_labels}:'.format(
            interaction=interaction_type, force_name=force_name,
            systems_labels=systems_labels)
        )
        for key in sorted(different_parameters.keys()):
            param1, param2 = different_parameters[key]
            parameters_diff = param1.pretty_format_diff(param2)
            diff_msg += '\n{} {}:\n{}'.format(interaction_type, key, parameters_diff)

    if diff_msg != '':
        diff_msg = ('A difference between {} was detected. '
                    'Details follow.').format(interaction_type) + diff_msg
        raise FailedParameterComparisonError(diff_msg + '\n', different_parameters)


@functools.singledispatch
def _get_force_parameters(force, system):
    """This function builds a _ParameterComparer representation of the
    force parameters that can be used for comparison to other forces
    with _compare_parameters.

    Each _ParameterComparer must be associated to a key, which is usually
    a tuple of atom indices (e.g., the atom index for nonbonded interactions,
    two atom indices for bonds and 1-4 exceptions, 3 indices for angles,
    and 4 for torsions).

    A single force can generate multiple representations (e.g., a NonbondedForce
    can generate one representations for particle parameters and one for
    parameter exceptions.

    Parameters
    ----------
    force
        The OpenMM Force object.
    system : simtk.openmm.System
        The System to which this force belongs to.

    Returns
    -------
    force_parameters : Dict[str, Dict[Hashable, _ParameterComparer]]
        The parameter representation. force_parameters[interaction_type]
        is a dictionary mapping a tuple of atom indices to the
        associated parameters. interaction_type can be any string
        identifying the type of parameters (e.g. "particle exception",
        "proper torsion", "bond").

    """
    raise NotImplementedError('Comparison between {}s is not currently '
                              'supported.'.format(type(force)))


@_get_force_parameters.register(openmm.HarmonicBondForce)
def _get_bond_force_parameters(force, _):
    """Implementation of _get_force_parameters for HarmonicBondForces."""
    # Build the dictionary of the parameters of a single force.
    force_parameters = {}
    for bond_idx in range(force.getNumBonds()):
        atom1, atom2, r0, k = force.getBondParameters(bond_idx)

        # Ignore bonds with 0.0 spring constant that don't affect the energy.
        if (k / k.unit) == 0.0:
            continue

        # Reorder the bond to have a canonical key.
        bond_key = tuple(sorted([atom1, atom2]))
        force_parameters[bond_key] = _ParametersComparer(r0=r0, k=k)

    return {'bond': force_parameters}


@_get_force_parameters.register(openmm.HarmonicAngleForce)
def _get_angle_force_parameters(force, _):
    """Implementation of _get_force_parameters for HarmonicAngleForces."""
    # Build the dictionary of the parameters of a single force.
    force_parameters = {}
    for angle_idx in range(force.getNumAngles()):
        atom1, atom2, atom3, theta0, k = force.getAngleParameters(angle_idx)

        # Ignore angles with 0.0 spring constant that don't affect the energy.
        if (k / k.unit) == 0.0:
            continue

        # Reorder the bond to have a canonical key.
        angle_key = min(atom1, atom3), atom2, max(atom1, atom3)
        force_parameters[angle_key] = _ParametersComparer(theta0=theta0, k=k)

    return {'angle': force_parameters}


@_get_force_parameters.register(openmm.NonbondedForce)
def _get_nonbonded_force_parameters(force, _):
    """Implementation of _get_force_parameters for NonbondedForces."""
    # Build the dictionary of the particle parameters of the force.
    particle_parameters = {}
    for particle_idx in range(force.getNumParticles()):
        charge, sigma, epsilon = force.getParticleParameters(particle_idx)

        # Ignore sigma parameter if epsilon is 0.0.
        particle_parameters[particle_idx] = _ParametersComparer(
            charge=charge, epsilon=epsilon)
        if (epsilon / epsilon.unit) != 0.0:
            particle_parameters[particle_idx].parameters['sigma'] = sigma

    # Build the dictionary representation of the particle exceptions.
    exception_parameters = {}
    for exception_idx in range(force.getNumExceptions()):
        atom1, atom2, chargeprod, sigma, epsilon = force.getExceptionParameters(exception_idx)

        # Reorder the atom indices to have a canonical key.
        exception_key = tuple(sorted([atom1, atom2]))
        # Ignore sigma parameter if epsilon is 0.0.
        exception_parameters[exception_key] = _ParametersComparer(
            charge=chargeprod, epsilon=epsilon)
        if (epsilon / epsilon.unit) != 0.0:
            exception_parameters[exception_key].parameters['sigma'] = sigma

    return {'particle': particle_parameters, 'particle exception': exception_parameters}


@_get_force_parameters.register(openmm.PeriodicTorsionForce)
def _get_torsion_force_parameters(force, system):
    """Implementation of _get_force_parameters for NonbondedForces."""
    # Find all bonds. We'll use this to distinguish
    # between proper and improper torsions.
    bond_set = _find_all_bonds(system)

    proper_parameters = {}
    improper_parameters = {}
    for torsion_idx in range(force.getNumTorsions()):
        atom1, atom2, atom3, atom4, periodicity, phase, k = force.getTorsionParameters(torsion_idx)

        # Ignore torsions that don't contribute to the energy.
        if (k / k.unit) == 0.0:
            continue

        torsion_key = [atom1, atom2, atom3, atom4]
        if len(set(torsion_key)) != 4:
            raise ValueError('Torsion {} is defined on less than 4 atoms: {}'.format(torsion_key))

        # Check if this is proper or not.
        torsion_bonds = [(torsion_key[i], torsion_key[i+1]) for i in range(3)]
        is_proper = all(bond in bond_set for bond in torsion_bonds)

        # Determine the canonical order of the torsion key.
        if is_proper:
            torsion_key = _get_proper_torsion_canonical_order(*torsion_key)
            force_parameters = proper_parameters
        else:
            torsion_key = _get_improper_torsion_canonical_order(bond_set, *torsion_key)
            force_parameters = improper_parameters

        # Update the dictionary representation.
        parameters = _ParametersComparer(periodicity=periodicity, phase=phase, k=k)
        # Each set of 4 atoms can have multiple torsion terms.
        try:
            force_parameters[torsion_key].parameters.append(parameters)
        except KeyError:
            force_parameters[torsion_key] = _TorsionParametersComparer(parameters)

    return {'proper torsion': proper_parameters, 'improper torsion': improper_parameters}


def _find_all_bonds(system):
    """Find all bonds in the given system

    Parameters
    ----------
    system : simtk.openmm.System
        A System object containing a HarmonicBondForce from which the
        bonds are inferred.

    Returns
    -------
    bond_set : Set[Tuple[int]]
        A set of pairs (atom_index1, atom_index2) associated to bonds.
        For each bond two pairs are created with the atom indices in
        inverse order.

    """
    # Find the force with information on the bonds.
    bond_force = [f for f in system.getForces() if isinstance(f, openmm.HarmonicBondForce)]
    assert len(bond_force) == 1
    bond_force = bond_force[0]

    # Create a set of all bonds.
    bond_set = set()
    for bond_idx in range(bond_force.getNumBonds()):
        atom1, atom2, _, _ = bond_force.getBondParameters(bond_idx)
        bond_set.add((atom1, atom2))
        bond_set.add((atom2, atom1))

    return bond_set


def _get_proper_torsion_canonical_order(i0, i1, i2, i3):
    """Create a unique order of the 4 atom indices of a proper torsion.

    The atom indices of a proper torsion are reordered so that the
    first atom is the smallest index.

    Parameters
    ----------
    i0, i1, i2, i3 : int
        Atom indices of the proper torsion.

    Returns
    -------
    j0, j1, j2, j3 : int
        Reordered atom indices of the proper torsion.

    """
    if i0 < i3:
        return i0, i1, i2, i3
    else:
        return i3, i2, i1, i0


def _get_improper_torsion_canonical_order(bond_set, i0, i1, i2, i3):
    """Create a unique order of the 4 atom indices of an improper torsion.

    Return j0, j1, j2, j3, where j0 is the central index and j1, j2, j3
    are in sorted() order. The centrality is determined by the maximum
    counts in the adjacency matrix.

    Parameters
    ----------
    bond_set : set
        The set of all bonds as determined by _find_all_bonds.
    i0, i1, i2, i3 : int,
        Atom indices of the improper torsion.

    Returns
    -------
    j0, j1, j2, j3 : int
        Reordered atom indices of the improper torsion, with j0 being
        the central index.

    """
    connections = np.zeros((4, 4))

    mapping = {i0: 0, i1: 1, i2: 2, i3: 3}
    inv_mapping = dict([(val, key) for key, val in mapping.items()])

    for (a, b) in itertools.combinations([i0, i1, i2, i3], 2):
        if (a, b) in bond_set:
            i, j = mapping[a], mapping[b]
            connections[i, j] += 1.
            connections[j, i] += 1.

    central_ind = connections.sum(0).argmax()
    central_ind = inv_mapping[central_ind]
    other_ind = sorted([i0, i1, i2, i3])
    other_ind.remove(central_ind)

    return central_ind, other_ind[0], other_ind[1], other_ind[2]


def compare_system_parameters(system1, system2, systems_labels=None):
    """Check that two OpenMM systems have the same parameters.

    Parameters
    ----------
    system1 : simtk.openmm.System
        The first system to compare.
    system2 : simtk.openmm.System
        The second system to compare.
    systems_labels : Tuple[str], optional
        A pair of strings with a meaningful name for the system. If
        specified, this will be included in the error message to
        improve its readability.

    Raises
    ------
    FailedParameterComparisonError
        If there are differences in the parameters of the two forces.
        The exceptions exposes an attribute ``different_parameters``
        with the different parameters.

    """
    # We need to perform some checks on the type and number of forces in the Systems.
    force_names1 = collections.Counter(f.__class__.__name__ for f in system1.getForces())
    force_names2 = collections.Counter(f.__class__.__name__ for f in system2.getForces())
    err_msg = 'Only systems having 1 force per type are supported.'
    assert set(force_names1.values()) == {1}, err_msg
    assert set(force_names2.values()) == {1}, err_msg

    # Check that the two systems have the same forces.
    err_msg = 'The two Systems have different forces: system1 {}, system2 {}'
    assert set(force_names1) == set(force_names2), err_msg.format(force_names1, force_names2)

    # Find all the pair of forces to compare.
    force_pairs = {force_name: [] for force_name in force_names1}
    for system in [system1, system2]:
        for force in system.getForces():
            force_pairs[force.__class__.__name__].append(force)

    # Compare all pairs of forces
    for force_name, (force1, force2) in force_pairs.items():
        parameters_force1 = _get_force_parameters(force1, system1)
        parameters_force2 = _get_force_parameters(force2, system2)
        for parameter_type, parameters1 in parameters_force1.items():
            parameters2 = parameters_force2[parameter_type]
            _compare_parameters(parameters1, parameters2,
                                interaction_type=parameter_type,
                                force_name=force_name,
                                systems_labels=systems_labels)


#=============================================================================================
# Utility functions to compare SMIRNOFF and AMBER force fields.
#=============================================================================================

def compare_amber_smirnoff(prmtop_filepath, inpcrd_filepath, forcefield, molecule):
    """
    Compare energies and parameters for OpenMM Systems/topologies created
    from an AMBER prmtop and crd versus from a SMIRNOFF forcefield file which
    should parametrize the same system with same parameters.

    Parameters
    ----------
    prmtop_filepath : str
        Path to the topology/parameter file in AMBER prmtop format
    inpcrd_filepath : str
        Path to the coordinates file in AMBER inpcrd or rst7 format
    forcefield : ForceField
        Force field instance used to create the system to compare.
    molecule : topology.molecule.Molecule
        The molecule object to test.

    Returns
    -------
    amber_energies : Dict[str, simtk.unit.Quantity]
        The potential energy of the AMBER system for each force type.
    forcefield_energies : Dict[str, simtk.unit.Quantity]
        The potential energy of the ForceField system for each force type.

    Raises
    ------
    FailedEnergyComparisonError
        If the potential energies of the two systems are different
        or if the two Systems in the contexts have forces divided into
        different force groups. The potential energies of the AMBER and
        ForceField systems can be accessed in the Exception through the
        attributes potential_energy1 and potential_energy2 respectively.
    FailedParameterComparisonError
        If there are differences in the parameters of the two systems.
        The exceptions exposes an attribute ``different_parameters``
        with the different parameters.

    """
    from openforcefield.topology import Topology

    # Create System from AMBER files. By default, contrarily to ForceField,
    # systems from AMBER files are created with removeCMMotion=True
    amber_system, openmm_topology, positions = create_system_from_amber(
        prmtop_filepath, inpcrd_filepath, removeCMMotion=False)
    box_vectors = amber_system.getDefaultPeriodicBoxVectors()

    # Create System from forcefield.
    openff_topology = Topology.from_openmm(openmm_topology, unique_molecules=[molecule])
    ff_system = forcefield.create_openmm_system(openff_topology)

    # Test energies and parameters.
    compare_system_parameters(amber_system, ff_system, systems_labels=('AMBER', 'SMIRNOFF'))
    amber_energies, forcefield_energies = compare_system_energies(
        amber_system, ff_system, positions, box_vectors)

    return amber_energies, forcefield_energies
