"""
Helper functions for producing `openmm.Force` objects for non-bonded terms.
"""
from collections import defaultdict
from typing import DefaultDict, Optional, Union

import openmm
from openff.toolkit import Molecule
from openff.units import unit as off_unit
from openff.units.openmm import to_openmm as to_openmm_quantity
from openmm import unit
from typing_extensions import TypeAlias

from openff.interchange import Interchange
from openff.interchange.common._nonbonded import ElectrostaticsCollection, vdWCollection
from openff.interchange.components.potentials import Collection
from openff.interchange.constants import _PME
from openff.interchange.exceptions import (
    CannotSetSwitchingFunctionError,
    InternalInconsistencyError,
    UnsupportedCutoffMethodError,
    UnsupportedExportError,
)
from openff.interchange.models import TopologyKey, VirtualSiteKey

_DATA_DICT: TypeAlias = dict[str, Union[None, str, bool, "Collection"]]

_MIXING_RULE_EXPRESSIONS: dict[str, str] = {
    "lorentz-berthelot": "sigma=(sigma1+sigma2)/2; epsilon=sqrt(epsilon1*epsilon2); ",
}


def _process_nonbonded_forces(
    interchange: "Interchange",
    system: openmm.System,
    combine_nonbonded_forces: bool = False,
) -> dict[Union[int, VirtualSiteKey], int]:
    """
    Process the non-bonded collections in an Interchange into corresponding openmm objects.

    This typically involves processing the vdW and Electrostatics sections of an Interchange object
    into a corresponding openmm.NonbondedForce (if `combine_nonbonded_forces=True`) or a
    collection of other forces (NonbondedForce, CustomNonbondedForce, CustomBondForce) if
    `combine_nonbondoed_forces=False`.

    """
    from openff.interchange.common._nonbonded import _NonbondedCollection

    for collection in interchange.collections.values():
        if isinstance(collection, _NonbondedCollection):
            break
    else:
        # If there are no non-bonded collections, assume here that there can be no virtual sites,
        # so just return an i-i mapping between OpenFF and OpenMM indices
        return {i: i for i in range(interchange.topology.n_atoms)}

    has_virtual_sites = "VirtualSites" in interchange.collections

    if has_virtual_sites:
        from openff.interchange.interop._virtual_sites import (
            _virtual_site_parent_molecule_mapping,
        )
        from openff.interchange.interop.openmm._virtual_sites import (
            _check_virtual_site_exclusion_policy,
        )

        virtual_sites = interchange["VirtualSites"]

        _check_virtual_site_exclusion_policy(virtual_sites)

        virtual_site_molecule_map: dict[
            VirtualSiteKey,
            int,
        ] = _virtual_site_parent_molecule_mapping(interchange)

        molecule_virtual_site_map: dict[int, list[VirtualSiteKey]] = defaultdict(list)

        for virtual_site_key, molecule_index in virtual_site_molecule_map.items():
            molecule_virtual_site_map[molecule_index].append(virtual_site_key)

    else:
        molecule_virtual_site_map = defaultdict(list)

    # Mapping between OpenFF "particles" and OpenMM particles (via index). OpenFF objects
    # (keys) are either atom indices (if atoms) or `VirtualSitesKey`s if virtual sites
    # openff_openmm_particle_map: Dict[Union[int, VirtualSiteKey], int] = dict()
    openff_openmm_particle_map = _add_particles_to_system(
        interchange,
        system,
        molecule_virtual_site_map,
    )

    # TODO: Process ElectrostaticsHandler.exception_potential
    has_vdw = False

    for name, collection in interchange.collections.items():
        if name == "vdW":
            has_vdw = True
            break
        if collection.is_plugin:
            # TODO: Here is where to detect an electrostatics plugin, if one ever exists
            if collection.acts_as == "vdW":
                has_vdw = True
                break

    if not has_vdw:
        if has_virtual_sites:
            raise UnsupportedExportError(
                "Virtual sites with no vdW handler not currently supported. If this use case is "
                "important to you, please raise an issue describing the functionality you wish to "
                "see.",
            )

        try:
            interchange["Electrostatics"]
        except LookupError:
            raise InternalInconsistencyError(
                "In a confused state, could not find any vdW interactions but also failed to find "
                "any electrostatics collection. This is a supported use case but should have been caught "
                "earlier in this function. Please file an issue with a minimal reproducing example.",
            )

    _data = _prepare_input_data(interchange)

    if combine_nonbonded_forces:
        _func = _create_single_nonbonded_force
    else:
        _func = _create_multiple_nonbonded_forces

    _func(
        _data,
        interchange,
        system,
        molecule_virtual_site_map,
        openff_openmm_particle_map,
    )

    return openff_openmm_particle_map


def _add_particles_to_system(
    interchange: "Interchange",
    system: openmm.System,
    molecule_virtual_site_map,
) -> dict[Union[int, VirtualSiteKey], int]:
    openff_openmm_particle_map: dict[Union[int, VirtualSiteKey], int] = dict()

    for molecule in interchange.topology.molecules:
        for atom in molecule.atoms:
            atom_index = interchange.topology.atom_index(atom)

            # Skip unit check for speed, trust that the toolkit reports mass in Dalton
            system_index = system.addParticle(mass=atom.mass.m)

            openff_openmm_particle_map[atom_index] = system_index

        for virtual_site_key in molecule_virtual_site_map[
            interchange.topology.molecule_index(molecule)
        ]:
            from openff.interchange.interop.openmm._virtual_sites import (
                _create_openmm_virtual_site,
                _create_virtual_site_object,
            )

            system_index = system.addParticle(mass=0.0)

            openff_openmm_particle_map[virtual_site_key] = system_index

            virtual_site_potential = interchange["VirtualSites"].potentials[
                interchange["VirtualSites"].key_map[virtual_site_key]
            ]

            virtual_site_object = _create_virtual_site_object(
                virtual_site_key,
                virtual_site_potential,
            )

            openmm_particle: openmm.VirtualSite = _create_openmm_virtual_site(
                virtual_site_object,
                openff_openmm_particle_map,
            )

            system.setVirtualSite(system_index, openmm_particle)

    return openff_openmm_particle_map


def _prepare_input_data(interchange: "Interchange") -> _DATA_DICT:
    try:
        vdw: "vdWCollection" = interchange["vdW"]
    except LookupError:
        for collection in interchange.collections.values():
            if collection.is_plugin:
                if collection.acts_as == "vdW":
                    # We can't be completely sure all plugins subclass out of vdWCollection here
                    vdw = collection  # type: ignore[assignment]
                    break
        else:
            vdw = None  # type: ignore[assignment]

    if vdw:
        vdw_cutoff: Optional[unit.Quanaity] = vdw.cutoff
        vdw_method: Optional[str] = getattr(vdw, "method", "cutoff").lower()
        mixing_rule: Optional[str] = getattr(vdw, "mixing_rule", None)
        vdw_expression: Optional[str] = vdw.expression.replace("**", "^")
    else:
        vdw_cutoff = None
        vdw_method = None
        mixing_rule = None
        vdw_expression = None

    try:
        electrostatics: "ElectrostaticsCollection" = interchange["Electrostatics"]
    except LookupError:
        electrostatics = None  # type: ignore[assignment]

    if electrostatics is None:
        electrostatics_method: Optional[str] = None
    else:
        if interchange.box is None:
            electrostatics_method = getattr(
                electrostatics,
                "nonperiodic_potential",
                "Coulomb",
            )
        else:
            electrostatics_method = getattr(electrostatics, "periodic_potential", _PME)

    return {
        "vdw_collection": vdw,
        "vdw_cutoff": vdw_cutoff,
        "vdw_method": vdw_method,
        "vdw_expression": vdw_expression,
        "mixing_rule": mixing_rule,
        "mixing_rule_expression": _MIXING_RULE_EXPRESSIONS.get(mixing_rule, "")
        if isinstance(mixing_rule, str)
        else None,
        "electrostatics_collection": electrostatics,
        "electrostatics_method": electrostatics_method,
        "periodic": interchange.box is None,
    }


def _create_single_nonbonded_force(
    data: dict,
    interchange: "Interchange",
    system: openmm.System,
    molecule_virtual_site_map: dict["Molecule", list[VirtualSiteKey]],
    openff_openmm_particle_map: dict[Union[int, VirtualSiteKey], int],
):
    """Create a single openmm.NonbondedForce from vdW/electrostatics/virtual site collections."""
    if data["mixing_rule"] not in ("lorentz-berthelot", None):
        raise UnsupportedExportError(
            "OpenMM's default NonbondedForce only supports Lorentz-Berthelot mixing rules."
            "Try setting `combine_nonbonded_forces=False`.",
        )

    if molecule_virtual_site_map in (None, dict()):
        has_virtual_sites = False
    elif all([len(v) == 0 for v in molecule_virtual_site_map.values()]):
        has_virtual_sites = False
    else:
        has_virtual_sites = True

    non_bonded_force = openmm.NonbondedForce()
    system.addForce(non_bonded_force)

    if interchange.box is None:
        if (data["vdw_method"] in ("cutoff", None)) and (
            data["electrostatics_method"] in ("Coulomb", None)
        ):
            non_bonded_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
            non_bonded_force.setUseDispersionCorrection(True)
            if data["vdw_cutoff"]:
                non_bonded_force.setCutoffDistance(
                    to_openmm_quantity(data["vdw_cutoff"]),
                )
        else:
            raise UnsupportedCutoffMethodError(
                f"Combination of non-bonded cutoff methods {data['vdw_method']} (vdW) and "
                f"{data['electrostatics_method']} (Electrostatics) not currently supported or "
                f"invalid with `combine_nonbonded_forces=True` and `.box={interchange.box}`.",
            )

    else:
        if data["vdw_method"] in ("cutoff", None) and data["electrostatics_method"] in (
            _PME,
            None,
        ):
            non_bonded_force.setNonbondedMethod(openmm.NonbondedForce.PME)
            non_bonded_force.setEwaldErrorTolerance(1.0e-4)
            non_bonded_force.setUseDispersionCorrection(True)
            if not data["vdw_cutoff"]:
                # With no vdW handler and/or ambiguous cutoff, cannot set it,
                # thereforce silently fall back to OpenMM's default. It's not
                # clear if this value matters with only (PME) charges and no
                # vdW interactions in the system.
                pass
            else:
                non_bonded_force.setCutoffDistance(
                    to_openmm_quantity(data["vdw_cutoff"]),
                )
        elif data["vdw_method"] == "pme" and data["electrostatics_method"] == _PME:
            non_bonded_force.setNonbondedMethod(openmm.NonbondedForce.LJPME)
            non_bonded_force.setEwaldErrorTolerance(1.0e-4)
        else:
            raise UnsupportedCutoffMethodError(
                f"Combination of non-bonded cutoff methods {data['vdw_method']} (vdW) and "
                "{data['electrostatics_method']} (Electrostatics) not currently supported or invalid with "
                f"`combine_nonbonded_forces=True` and `.box={interchange.box}`.",
            )

    if data["electrostatics_collection"] is not None:
        if has_virtual_sites:
            partial_charges = data[
                "electrostatics_collection"
            ].charges_with_virtual_sites
        else:
            partial_charges = data["electrostatics_collection"].charges

    # mapping between (openmm) index of each atom and the (openmm) index of each virtual particle
    #   of that parent atom (if any)
    # if no virtual sites at all, this remains an empty dict
    parent_virtual_particle_mapping: DefaultDict[int, list[int]] = defaultdict(list)

    vdw = data["vdw_collection"]

    for molecule in interchange.topology.molecules:
        for atom in molecule.atoms:
            non_bonded_force.addParticle(0.0, 1.0, 0.0)

            atom_index = interchange.topology.atom_index(atom)

            top_key = TopologyKey(atom_indices=(atom_index,))

            if data["electrostatics_collection"] is not None:
                partial_charge = partial_charges[top_key].m_as(off_unit.e)
            else:
                partial_charge = 0.0

            if vdw is not None:
                pot_key = vdw.key_map[top_key]

                sigma = vdw.potentials[pot_key].parameters["sigma"]
                epsilon = vdw.potentials[pot_key].parameters["epsilon"]

                sigma = sigma.m_as(off_unit.nanometer)
                epsilon = epsilon.m_as(off_unit.kilojoule / off_unit.mol)
            else:
                sigma = unit.Quantity(0.0, unit.nanometer)
                epsilon = unit.Quantity(0.0, unit.kilojoules_per_mole)

            openmm_atom_index = openff_openmm_particle_map[atom_index]

            non_bonded_force.setParticleParameters(
                openmm_atom_index,
                partial_charge,
                sigma,
                epsilon,
            )

        if has_virtual_sites:
            molecule_index = interchange.topology.molecule_index(molecule)
        else:
            continue

        for virtual_site_key in molecule_virtual_site_map[molecule_index]:
            # TODO: Move this function to openff/interchange/interop/_particles.py ?
            from openff.interchange.interop.openmm._virtual_sites import (
                _create_openmm_virtual_site,
                _create_virtual_site_object,
            )

            _potential_key = interchange["VirtualSites"].key_map[virtual_site_key]
            virtual_site_potential = interchange["VirtualSites"].potentials[
                _potential_key
            ]
            virtual_site_object = _create_virtual_site_object(
                virtual_site_key,
                virtual_site_potential,
            )

            openmm_particle: openmm.VirtualSite = _create_openmm_virtual_site(
                virtual_site_object,
                openff_openmm_particle_map,
            )

            vdw = interchange["vdW"]
            coul = interchange["Electrostatics"]

            vdw_key = vdw.key_map.get(virtual_site_key)
            coul_key = coul.key_map.get(virtual_site_key)
            if vdw_key is None or coul_key is None:
                raise InternalInconsistencyError(
                    f"Virtual site {virtual_site_key} is not associated with any "
                    "vdW and/or electrostatics interactions",
                )

            charge_increments = coul.potentials[coul_key].parameters[
                "charge_increments"
            ]
            charge = to_openmm_quantity(-sum(charge_increments))

            vdw_parameters = vdw.potentials[vdw_key].parameters
            sigma = to_openmm_quantity(vdw_parameters["sigma"])
            epsilon = to_openmm_quantity(vdw_parameters["epsilon"])

            system_index: int = openff_openmm_particle_map[virtual_site_key]
            force_index = non_bonded_force.addParticle(charge, sigma, epsilon)

            if system_index != force_index:
                raise InternalInconsistencyError(
                    "Mismatch in system and force indexing",
                )

            parent_atom_index = openff_openmm_particle_map[
                virtual_site_object.orientations[0]
            ]

            parent_virtual_particle_mapping[parent_atom_index].append(force_index)

            system.setVirtualSite(system_index, openmm_particle)

    _create_exceptions(
        data,
        non_bonded_force,
        interchange,
        openff_openmm_particle_map,
        parent_virtual_particle_mapping,
    )

    _apply_switching_function(data["vdw_collection"], non_bonded_force)


def _create_exceptions(
    data: dict,
    non_bonded_force: openmm.NonbondedForce,
    interchange: "Interchange",
    openff_openmm_particle_map: dict,
    parent_virtual_particle_mapping: DefaultDict[int, list[int]],
):
    # The topology indices reported by toolkit methods must be converted to openmm indices
    bonds = [
        sorted(
            openff_openmm_particle_map[interchange.topology.atom_index(a)]
            for a in bond.atoms
        )
        for bond in interchange.topology.bonds
    ]

    coul_14 = getattr(data["electrostatics_collection"], "scale_14", 1.0)
    vdw_14 = getattr(data["vdw_collection"], "scale_14", 1.0)

    # First, create all atom-atom exceptions according to the conventional pattern
    non_bonded_force.createExceptionsFromBonds(
        bonds=bonds,
        coulomb14Scale=coul_14,
        lj14Scale=vdw_14,
    )

    # Faster to loop through exceptions and look up parents than opposite
    if parent_virtual_particle_mapping not in (None, dict()):
        # First add exceptions between each virtual particle and parent atom
        for (
            parent,
            virtual_particles_of_this_parent,
        ) in parent_virtual_particle_mapping.items():
            for virtual_particle in virtual_particles_of_this_parent:
                # These indices are of the OpenMM particles, so no need to use map
                non_bonded_force.addException(
                    parent,
                    virtual_particle,
                    0.0,
                    1.0,
                    0.0,
                    True,
                )

        for exception_index in range(non_bonded_force.getNumExceptions()):
            # These particles should only be atoms in this loop
            (
                p1,
                p2,
                charge_prod,
                _,
                epsilon,
            ) = non_bonded_force.getExceptionParameters(exception_index)
            for virtual_particle_of_p1 in parent_virtual_particle_mapping[p1]:
                # If this iterable is not empty, add an exception between p1's virtual
                # particle and the "other" atom in p1's exception
                if virtual_particle_of_p1 == p2:
                    continue

                if charge_prod._value == epsilon._value == 0.0:
                    non_bonded_force.addException(
                        particle1=virtual_particle_of_p1,
                        particle2=p2,
                        chargeProd=0.0,
                        sigma=1.0,
                        epsilon=0.0,
                        replace=True,
                    )
                else:
                    # TODO: Pass mixing rule into Decide on best logic for inheriting scaled 1-4 interactions
                    v1_parameters = non_bonded_force.getParticleParameters(
                        virtual_particle_of_p1,
                    )
                    p2_parameters = non_bonded_force.getParticleParameters(p2)
                    non_bonded_force.addException(
                        particle1=virtual_particle_of_p1,
                        particle2=p2,
                        chargeProd=v1_parameters[0] * p2_parameters[0],
                        sigma=(v1_parameters[1] + p2_parameters[1]) * 0.5,
                        epsilon=(v1_parameters[2] * p2_parameters[2]) ** 0.5,
                    )
            for virtual_particle_of_p2 in parent_virtual_particle_mapping[p2]:
                # If this iterable is not empty, add an exception between p1's virtual
                # particle and the "other" atom in p1's exception
                if virtual_particle_of_p2 == p1:
                    continue

                if charge_prod._value == epsilon._value == 0.0:
                    non_bonded_force.addException(
                        particle1=virtual_particle_of_p2,
                        particle2=p1,
                        chargeProd=0.0,
                        sigma=1.0,
                        epsilon=0.0,
                        replace=True,
                    )
                else:
                    # TODO: Pass mixing rule into Decide on best logic for inheriting scaled 1-4 interactions
                    v2_parameters = non_bonded_force.getParticleParameters(
                        virtual_particle_of_p2,
                    )
                    p1_parameters = non_bonded_force.getParticleParameters(p1)
                    non_bonded_force.addException(
                        particle1=virtual_particle_of_p2,
                        particle2=p1,
                        chargeProd=v2_parameters[0] * p1_parameters[0],
                        sigma=(v2_parameters[1] + p1_parameters[1]) * 0.5,
                        epsilon=(v2_parameters[2] * p1_parameters[2]) ** 0.5,
                    )


def _create_multiple_nonbonded_forces(
    data: dict,
    interchange: "Interchange",
    system: openmm.System,
    molecule_virtual_site_map: dict,
    openff_openmm_particle_map: dict[Union[int, VirtualSiteKey], int],
):
    from openff.interchange.components.toolkit import _get_14_pairs

    if molecule_virtual_site_map in (None, dict()):
        has_virtual_sites = False
    elif all([len(v) == 0 for v in molecule_virtual_site_map.values()]):
        has_virtual_sites = False
    else:
        has_virtual_sites = True

    vdw_force = _create_vdw_force(
        data,
        interchange,
        molecule_virtual_site_map,
        has_virtual_sites,
    )

    electrostatics_force: openmm.NonbondedForce = _create_electrostatics_force(
        data,
        interchange,
        molecule_virtual_site_map,
        has_virtual_sites,
        openff_openmm_particle_map,
    )

    _set_particle_parameters(
        data,
        vdw_force,
        electrostatics_force,
        interchange,
        has_virtual_sites,
        molecule_virtual_site_map,
        openff_openmm_particle_map,
    )

    # Attempting to match the value used internally by OpenMM; The source of this value is likely
    # https://github.com/openmm/openmm/issues/1149#issuecomment-250299854
    # 1 / * (4pi * eps0) * elementary_charge ** 2 / nanometer ** 2
    coul_const = 138.935456  # kJ/nm

    if vdw_force is not None:
        vdw = data["vdw_collection"]

        if vdw.is_plugin:
            vdw_14_force = openmm.CustomBondForce(
                _get_scaled_potential_function(data["vdw_expression"]),
            )

            # feed in r_min1, epsilon1, ..., r_min2, epsilon2, ... each as individual parameters
            for index in [1, 2]:
                for parameter in vdw.potential_parameters():
                    vdw_14_force.addPerBondParameter(f"{parameter}{str(index)}")

            vdw_14_force.addGlobalParameter("scale14", vdw.scale_14)

            for global_parameter in vdw.global_parameters():
                vdw_14_force.addGlobalParameter(
                    global_parameter,
                    getattr(vdw, global_parameter).m,
                )

            for term, value in data["vdw_collection"].pre_computed_terms().items():
                vdw_14_force.addGlobalParameter(term, value)

        else:
            vdw_14_force = openmm.CustomBondForce(data["vdw_expression"])

            for parameter in data["vdw_collection"].potential_parameters():
                vdw_14_force.addPerBondParameter(parameter)

        vdw_14_force.setUsesPeriodicBoundaryConditions(interchange.box is not None)

    else:
        vdw_14_force = None

    coul_14_force = openmm.CustomBondForce(f"{coul_const}*qq/r")
    coul_14_force.addPerBondParameter("qq")
    coul_14_force.setUsesPeriodicBoundaryConditions(interchange.box is not None)

    coul_14, vdw_14 = _get_14_scaling_factors(data)

    openmm_pairs = list()

    for atom1, atom2 in _get_14_pairs(interchange.topology):
        openff_indices = (
            interchange.topology.atom_index(atom1),
            interchange.topology.atom_index(atom2),
        )

        openmm_indices = (
            openff_openmm_particle_map[openff_indices[0]],
            openff_openmm_particle_map[openff_indices[1]],
        )

        openmm_pairs.append(openmm_indices)

    if electrostatics_force is not None:
        for i in range(electrostatics_force.getNumExceptions()):
            (p1, p2, _, _, _) = electrostatics_force.getExceptionParameters(i)

            if (p1, p2) in openmm_pairs or (p2, p1) in openmm_pairs:
                if vdw_force is not None:
                    if data["vdw_collection"].is_plugin:
                        # Since we fed in in r_min1, epsilon1, ..., r_min2, epsilon2, ...
                        # each as individual parameters, we need to prepare a list of
                        # length 2 * len(potential_parameters) to this constructor
                        parameters1 = vdw_force.getParticleParameters(p1)
                        parameters2 = vdw_force.getParticleParameters(p2)
                        vdw_14_force.addBond(p1, p2, [*parameters1, *parameters2])

                    else:
                        # Look up the vdW parameters for each particle
                        sig1, eps1 = vdw_force.getParticleParameters(p1)
                        sig2, eps2 = vdw_force.getParticleParameters(p2)

                        # manually compute ...
                        sig_14 = (sig1 + sig2) * 0.5
                        eps_14 = (eps1 * eps2) ** 0.5 * vdw_14

                        # ... and set the 1-4 interactions
                        vdw_14_force.addBond(p1, p2, [sig_14, eps_14])

                # Look up the partial charges for each particle
                q1 = electrostatics_force.getParticleParameters(p1)[0]
                q2 = electrostatics_force.getParticleParameters(p2)[0]

                # manually compute ...
                qq = q1 * q2 * coul_14

                # ... and set the 1-4 interactions
                coul_14_force.addBond(p1, p2, [qq])

            if vdw_force is not None:
                vdw_force.addExclusion(p1, p2)

            if electrostatics_force is not None:
                electrostatics_force.setExceptionParameters(i, p1, p2, 0.0, 0.0, 0.0)

    for force in [vdw_force, electrostatics_force, vdw_14_force, coul_14_force]:
        if force is not None:
            system.addForce(force)

    if vdw_force is not None and electrostatics_force is not None:
        if (vdw_force.getNonbondedMethod() > 0) ^ (
            electrostatics_force.getNonbondedMethod() > 0
        ):
            raise UnsupportedCutoffMethodError(
                "When using `openmm.CustomNonbondedForce`, vdW and electrostatics cutoff methods "
                "must agree on whether or not periodic boundary conditions should be used."
                f"OpenMM will throw an error. Found vdw method {vdw_force.getNonbondedMethod()}, "
                f"and electrostatics method {electrostatics_force.getNonbondedMethod()}, ",
            )


def _create_vdw_force(
    data: _DATA_DICT,
    interchange: "Interchange",
    molecule_virtual_site_map: dict[int, list[VirtualSiteKey]],
    has_virtual_sites: bool,
) -> Optional[openmm.CustomNonbondedForce]:
    vdw_collection: Optional["vdWCollection"] = data["vdw_collection"]  # type: ignore[assignment]

    if vdw_collection is None:
        return None

    vdw_expression: str = data["vdw_expression"]  # type: ignore[assignment]
    mixing_rule_expression: str = data["mixing_rule_expression"]  # type: ignore[assignment]

    vdw_force = openmm.CustomNonbondedForce(
        f"{vdw_expression}"
        if mixing_rule_expression in (None, "")
        else f"{vdw_expression}; {mixing_rule_expression}",
    )

    for potential_parameter in vdw_collection.potential_parameters():
        vdw_force.addPerParticleParameter(potential_parameter)

    if vdw_collection.is_plugin:
        # TODO: Move this block outside of the plugin conditional
        for global_parameter in vdw_collection.global_parameters():
            vdw_force.addGlobalParameter(
                global_parameter,
                getattr(vdw_collection, global_parameter).m,
            )

        for term, value in vdw_collection.pre_computed_terms().items():
            vdw_force.addGlobalParameter(term, value)

    for molecule in interchange.topology.molecules:
        for _ in molecule.atoms:
            vdw_force.addParticle(vdw_collection.default_parameter_values())

        if has_virtual_sites:
            molecule_index = interchange.topology.molecule_index(molecule)
            for _ in molecule_virtual_site_map[molecule_index]:
                vdw_force.addParticle(vdw_collection.default_parameter_values())

    if data["vdw_method"] == "cutoff":
        if interchange.box is None:
            vdw_force.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)
        else:
            vdw_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        vdw_force.setUseLongRangeCorrection(True)
        vdw_force.setCutoffDistance(to_openmm_quantity(data["vdw_cutoff"]))

        _apply_switching_function(vdw_collection, vdw_force)

    elif data["vdw_method"] == "pme":
        if interchange.box is None:
            raise UnsupportedCutoffMethodError(
                "vdW method pme/ljpme is not valid for non-periodic systems.",
            )
        else:
            raise UnsupportedCutoffMethodError(
                "LJ-PME with split non-bonded forces is not supported due to openmm.CustomNonbondedForce "
                "not supporting PME. If also using PME electrostatics, try `combine_nonbonded_forces=True`,  "
                "which should produce a single force with NonbondedForce.LJPME, which uses PME for both "
                "electrostatics and LJ forces tersm. If your use case would benenfit from split non-bonded "
                "forces with LJPME, please file an feature request.",
            )

    return vdw_force


def _create_electrostatics_force(
    data: _DATA_DICT,
    interchange: "Interchange",
    molecule_virtual_site_map: dict[int, list[VirtualSiteKey]],
    has_virtual_sites: bool,
    openff_openmm_particle_map,
) -> Optional[openmm.NonbondedForce]:
    if data["electrostatics_collection"] is None:
        return None

    electrostatics_force = openmm.NonbondedForce()

    # mapping between (openmm) index of each atom and the (openmm) index of each virtual particle
    #   of that parent atom (if any)
    # if no virtual sites at all, this remains an empty dict
    parent_virtual_particle_mapping: DefaultDict[int, list[int]] = defaultdict(list)

    for molecule in interchange.topology.molecules:
        for _ in molecule.atoms:
            electrostatics_force.addParticle(0.0, 1.0, 0.0)

        if has_virtual_sites:
            molecule_index = interchange.topology.molecule_index(molecule)
            for virtual_site_key in molecule_virtual_site_map[molecule_index]:
                force_index = electrostatics_force.addParticle(0.0, 1.0, 0.0)

                parent_atom_index = virtual_site_key.orientation_atom_indices[0]
                parent_virtual_particle_mapping[
                    openff_openmm_particle_map[parent_atom_index]
                ].append(force_index)

    if data["electrostatics_method"] == "reaction-field":
        raise UnsupportedExportError(
            "Reaction field electrostatics not supported. If this use case is important to you, "
            "please raise an issue describing the scope of functionality you would like to use.",
        )

    elif data["electrostatics_method"] == _PME:
        electrostatics_force.setNonbondedMethod(openmm.NonbondedForce.PME)
        electrostatics_force.setEwaldErrorTolerance(1.0e-4)
        electrostatics_force.setUseDispersionCorrection(True)
        if data["vdw_cutoff"] is not None:
            # All nonbonded forces must use the same cutoff, even though PME doesn't have a cutoff
            electrostatics_force.setCutoffDistance(
                to_openmm_quantity(data["vdw_cutoff"]),
            )
    elif data["electrostatics_method"] == "Coulomb":
        if interchange.box is None:
            electrostatics_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
        else:
            raise UnsupportedCutoffMethodError(
                f"Electrostatics method {data['electrostatics_method']} ambiguous with a periodic system.",
            )

    elif data["electrostatics_method"] == "cutoff":
        if interchange.box is None:
            electrostatics_force.setNonbondedMethod(
                openmm.NonbondedForce.CutoffNonPeriodic,
            )
            electrostatics_force.setCutoffDistance(2.0)
        else:
            electrostatics_force.setNonbondedMethod(
                openmm.NonbondedForce.CutoffPeriodic,
            )
            electrostatics_force.setCutoffDistance(0.9)

        # TODO: Wire through electrostatics cutoff?
    else:
        raise UnsupportedCutoffMethodError(
            f"Electrostatics method {data['electrostatics_method']} not supported",
        )

    _create_exceptions(
        data,
        electrostatics_force,
        interchange,
        openff_openmm_particle_map,
        parent_virtual_particle_mapping,
    )
    return electrostatics_force


def _set_particle_parameters(
    data: _DATA_DICT,
    vdw_force: openmm.CustomNonbondedForce,
    electrostatics_force: openmm.NonbondedForce,
    interchange: "Interchange",
    has_virtual_sites: bool,
    molecule_virtual_site_map: dict[int, list[VirtualSiteKey]],
    openff_openmm_particle_map: dict[Union[int, VirtualSiteKey], int],
):
    if electrostatics_force is not None:
        electrostatics: ElectrostaticsCollection = data[  # type: ignore[assignment]
            "electrostatics_collection"
        ]

        partial_charges = getattr(
            electrostatics,
            "charges_with_virtual_sites" if has_virtual_sites else "charges",
        )

    else:
        partial_charges = None

    vdw: "vdWCollection" = data["vdw_collection"]  # type: ignore[assignment]

    for molecule in interchange.topology.molecules:
        for atom in molecule.atoms:
            atom_index = interchange.topology.atom_index(atom)
            particle_index = openff_openmm_particle_map[atom_index]
            # TODO: Actually process virtual site vdW parameters here

            top_key = TopologyKey(atom_indices=(atom_index,))

            if partial_charges is not None:
                partial_charge = partial_charges[top_key].m_as(off_unit.e)
            else:
                partial_charge = 0.0

            if vdw is not None:
                pot_key = vdw.key_map[top_key]

                if vdw.is_plugin:
                    parameters = vdw.potentials[pot_key].parameters

                    if hasattr(vdw, "modify_parameters"):
                        # This method strips units ..
                        parameters = vdw.modify_parameters(parameters)
                    else:
                        # so manually strip them if the method is not present
                        parameters = {key: val.m for key, val in parameters.items()}

                else:
                    sigma = vdw.potentials[pot_key].parameters["sigma"]
                    epsilon = vdw.potentials[pot_key].parameters["epsilon"]

                    sigma = sigma.m_as(off_unit.nanometer)
                    epsilon = epsilon.m_as(off_unit.kilojoule / off_unit.mol)
            else:
                sigma = unit.Quantity(0.0, unit.nanometer)
                epsilon = unit.Quantity(0.0, unit.kilojoules_per_mole)

            if vdw_force is not None:
                if vdw.is_plugin:
                    vdw_force.setParticleParameters(
                        particle_index,
                        [*parameters.values()],
                    )
                else:
                    vdw_force.setParticleParameters(particle_index, [sigma, epsilon])

            if electrostatics_force is not None:
                electrostatics_force.setParticleParameters(
                    particle_index,
                    partial_charge,
                    0.0,
                    0.0,
                )

        for virtual_site_key in molecule_virtual_site_map[
            interchange.topology.molecule_index(molecule)
        ]:
            # TODO: Also set the vdW parameters, including figuring out how to wire up
            #       custom non-bonded functional forms
            particle_index = openff_openmm_particle_map[virtual_site_key]

            partial_charge = partial_charges[virtual_site_key].m_as(off_unit.e)

            if electrostatics_force is not None:
                electrostatics_force.setParticleParameters(
                    particle_index,
                    partial_charges[virtual_site_key].m_as(off_unit.e),
                    0.0,
                    0.0,
                )


def _get_14_scaling_factors(data: _DATA_DICT) -> tuple[float, float]:
    if data.get("electrostatics_collection", None) is not None:
        coul_14 = data["electrostatics_collection"].scale_14  # type: ignore[union-attr]
    else:
        coul_14 = 1.0

    if data.get("vdw_collection", None) is not None:
        vdw_14 = data["vdw_collection"].scale_14  # type: ignore[union-attr]
    else:
        vdw_14 = 1.0

    return coul_14, vdw_14


def _apply_switching_function(
    vdw_collection: "vdWCollection",
    force: openmm.NonbondedForce,
):
    if not hasattr(force, "setUseSwitchingFunction"):
        raise CannotSetSwitchingFunctionError(
            "Attempting to set switching function on an OpenMM force that does nont support it."
            f"Passed force of type {type(force)}.",
        )
    if getattr(vdw_collection, "switch_width", None) is None:
        force.setUseSwitchingFunction(False)
    elif vdw_collection.switch_width.m == 0.0:
        force.setUseSwitchingFunction(False)
    else:
        switching_distance = to_openmm_quantity(
            vdw_collection.cutoff - vdw_collection.switch_width,
        )

        if switching_distance._value < 0:
            raise UnsupportedCutoffMethodError(
                "Found a `switch_width` greater than the cutoff distance. It's not clear "
                "what this means and it's probably invalid. Found "
                f"switch_width{vdw_collection.switch_width} and cutoff {vdw_collection.cutoff}",
            )

        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(switching_distance)


def _get_scaled_potential_function(potential: str) -> str:
    # https://github.com/openforcefield/smirnoff-plugins/blob/979907fc93e8b6217d40e6f469809c230b86b012/smirnoff_plugins/handlers/nonbonded.py#L63-L83
    split_potential = potential.split(";")

    for i, expression in enumerate(split_potential):
        if "=" not in expression:
            # This is the actual expression, final energy so modify
            split_potential[i] = f"({expression})*scale14"
            break

    return ";".join(split_potential)
