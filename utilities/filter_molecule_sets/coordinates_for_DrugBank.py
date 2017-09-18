from openforcefield.utils import utils
from openeye import oechem
from openeye import oeomega
import openmoltools
import multiprocessing


def genConfs(c_mol, ofsff, ofsTri, index):
    # set omega settings
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(False)
    omega.SetEnergyWindow(15.0)
    strict_stereo = True 
    omega.SetStrictStereo(strict_stereo)
    omega.SetSampleHydrogens(True)
    omega.SetStrictAtomTypes(True)

    mol = oechem.OEMol(c_mol)
    status = omega(mol)

    if status:
        # change title
        mol.SetTitle('DrugBank_%i' % index)
        # save forcefield type
        mol1 = oechem.OEMol(mol)
        oechem.OETriposAtomNames(mol1)
        oechem.OEWriteConstMolecule(ofsff, mol1)

        # save Tripos atom types
        mol2 = oechem.OEMol(mol)
        oechem.OETriposAtomTypeNames(mol2)
        oechem.OEWriteConstMolecule(ofsTri, mol2)

    return status

flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield

in_file = utils.get_data_filename('molecules/DrugBank_atyped.oeb')
ff_out = 'DrugBank_ff.mol2'
tripos_out = 'DrugBank_tripos.mol2'
failed_out = 'DrugBank_no3D.mol2'

# open files
ofsff = oechem.oemolostream()
ofsff.SetFlavor(oechem.OEFormat_MOL2, flavor)
ofsff.open(ff_out)

ofsTri = oechem.oemolostream()
ofsTri.SetFlavor(oechem.OEFormat_MOL2, flavor)
ofsTri.open(tripos_out)

ofsFail = oechem.oemolostream()
ofsFail.SetFlavor(oechem.OEFormat_MOL2, flavor)
ofsFail.open(failed_out)

success = 0
time_out = 0
conf_fail = 0
index = 0

ifs = oechem.oemolistream(in_file) 
ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)

c_mol = oechem.OECreateOEGraphMol()
while oechem.OEReadMolecule(ifs, c_mol):
    index += 1
    # process molecules individually, storing less
    p = multiprocessing.Process(target=genConfs, args=(c_mol,ofsff, ofsTri, index,))
    p.start()
    p.join(24)
    if p.is_alive():
        print("TIMED OUT %s" % oechem.OECreateIsoSmiString(c_mol))
        oechem.OEWriteConstMolecule(ofsFail, oechem.OEMol(c_mol))
        time_out += 1
        p.terminate()
        p.join()
    elif p.exitcode:
        success += 1
        p.terminate()
        p.join()
    else:
        print("CONF FAIL %s" % oechem.OECreateIsoSmiString(c_mol))
        oechem.OEWriteConstMolecule(ofsFail, oechem.OEMol(c_mol))
        conf_fail += 1
        p.terminate()
        p.join()

# Print data
print("Success %i out of %i" % (success, index))
print("%i timed out" % time_out)
print("%i failed during conformation generation" % conf_fail)

# close files
ofsff.close()
ofsTri.close()
ofsFail.close()
