import openmm

from openff.interchange import Interchange
from openff.interchange.exceptions import UnsupportedExportError


def _process_gbsa(
    interchange: "Interchange",
    system: "openmm.System",
):
    import openmm.app
    import openmm.unit
    from openff.units import unit

    try:
        collection = interchange.collections["GBSA"]
    except KeyError:
        return

    existing_forces = [
        force
        for force in system.getForces()
        if isinstance(force, (openmm.NonbondedForce, openmm.CustomNonbondedForce))
    ]

    if len(existing_forces) != 1:
        raise UnsupportedExportError(
            "GBSA implementation assumes exactly one NonbondedForce or CustomNonbondedForce is present."
            f"Found these ({len(existing_forces)=}) forces: {existing_forces=}",
        )

    non_bonded_force = existing_forces[0]

    if non_bonded_force.getNonbondedMethod() == openmm.NonbondedForce.NoCutoff:
        amber_cutoff = None
    else:
        amber_cutoff = non_bonded_force.getCutoffDistance()

    if collection.gb_model == "OBC2":
        force = openmm.GBSAOBCForce()
    elif collection.gb_model in ["OBC1", "HCT"]:
        if collection.gb_model == "HCT":
            force_type = openmm.app.internal.customgbforces.GBSAHCTForce
        elif collection.gb_model == "OBC1":
            force_type = openmm.app.internal.customgbforces.GBSAOBC1Force

        force = force_type(
            solventDielectric=collection.solvent_dielectric,
            soluteDielectric=collection.solute_dielectric,
            SA=collection.sa_model,
            cutoff=amber_cutoff,
            kappa=0,
        )

    system.addForce(force)

    if amber_cutoff is not None:
        force.setCutoffDistance(amber_cutoff)

    if non_bonded_force.usesPeriodicBoundaryConditions():
        force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    else:
        force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

    if collection.gb_model == "OBC2":
        force.setSolventDielectric(collection.solvent_dielectric)
        force.setSoluteDielectric(collection.solute_dielectric)
        force.setSurfaceAreaEnergy(
            0
            if collection.sa_model is None
            else collection.surface_area_penalty.to_openmm(),
        )

    for topology_key, potential_key in collection.key_map.items():
        charge, *_ = non_bonded_force.getParticleParameters(
            topology_key.atom_indices[0],
        )
        _parameters = collection.potentials[potential_key].parameters

        parameters = [
            charge,
            _parameters["radius"].m_as(unit.nanometer),
            _parameters["scale"],
        ]

        if collection.gb_model == "OBC2":
            force.addParticle(*parameters)
        else:
            force.addParticle(parameters)

    if collection.gb_model != "OBC2":
        force.finalize()
