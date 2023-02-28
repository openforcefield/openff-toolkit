"""
Tests for forcefield class

"""
import copy
import itertools
import os
from tempfile import NamedTemporaryFile

import numpy as np
import openmm
import pytest
from numpy.testing import assert_almost_equal
from openff.units import unit
from openff.units.openmm import from_openmm, to_openmm
from openmm import NonbondedForce, Platform, XmlSerializer, app
from openmm import unit as openmm_unit
from pydantic import ValidationError

from openff.toolkit.tests.create_molecules import (
    create_acetate,
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
    create_water,
)
from openff.toolkit.tests.utils import (
    compare_partial_charges,
    get_14_scaling_factors,
    requires_openeye,
    requires_openeye_mol2,
    requires_rdkit,
    unimplemented_interchange,
)
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import (
    ElectrostaticsHandler,
    ForceField,
    LibraryChargeHandler,
    ParameterHandler,
    ToolkitAM1BCCHandler,
    XMLParameterIOHandler,
    get_available_force_fields,
    vdWHandler,
)
from openff.toolkit.utils import (
    AmberToolsToolkitWrapper,
    BuiltInToolkitWrapper,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
    get_data_file_path,
)
from openff.toolkit.utils.exceptions import (
    ChargeMethodUnavailableError,
    IncompatibleParameterError,
    ParameterLookupError,
    SMIRNOFFAromaticityError,
    SMIRNOFFSpecError,
    SMIRNOFFSpecUnimplementedError,
    SMIRNOFFVersionError,
)

XML_FF_GENERICS = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Bonds version="0.4" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
    <Bond smirks="[*:1]~[*:2]" id="b1" k="680.0 * kilocalories_per_mole/angstrom**2" length="1.09 * angstrom"/>
  </Bonds>
  <Angles version="0.3">
    <Angle smirks="[*:1]~[*:2]~[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]~[*:2]~[*:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[*:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1"  switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[*:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
"""

xml_without_section_version = """<?xml version="1.0" encoding="ASCII"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ToolkitAM1BCC/>
</SMIRNOFF>
"""

xml_simple_ff = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Bonds version="0.3">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0 * kilocalories_per_mole/angstrom**2" length="1.526 * angstrom"/>
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0 * kilocalories_per_mole/angstrom**2" length="1.51 * angstrom"/>
    <Bond smirks="[#6X3:1]-[#6X3:2]" id="b3" k_bondorder1="100.0*kilocalories_per_mole/angstrom**2" k_bondorder2="4000.0*kilocalories_per_mole/angstrom**2" length="1.5 * angstrom"/>
  </Bonds>
  <Angles version="0.3">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5 * degree" id="a2" k="70.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree" periodicity2="2" phase2="180.0 * degree" idivf2="1" k2="0.250 * kilocalories_per_mole" periodicity3="1" phase3="180.0 * degree" idivf3="1" k3="0.200 * kilocalories_per_mole"/>
    <Proper smirks="[*:1]:[#6X4:2]~[#6X4:3]:[*:4]" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157 * kilocalories_per_mole" id="n2" rmin_half="1.4870 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
"""

xml_ff_w_comments = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <!-- SMIRNOFF (SMIRKS Native Open Force Field) template file -->
  <Date>2018-07-14</Date>
  <Author>C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, UC Irvine; D. L. Mobley, UC Irvine</Author>
  <!-- This file is meant for processing via openff.toolkit.typing.engines.smirnoff -->
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the bond energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Bonds version="0.3">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0 * kilocalories_per_mole/angstrom**2" length="1.526 * angstrom" />
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0 * kilocalories_per_mole/angstrom**2" length="1.51 * angstrom"/>
  </Bonds>
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the angle energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Angles version="0.3">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5 * degree" id="a2" k="70.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree" periodicity2="2" phase2="180.0 * degree" idivf2="1" k2="0.250 * kilocalories_per_mole" periodicity3="1" phase3="180.0 * degree" idivf3="1" k3="0.200 * kilocalories_per_mole"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157 * kilocalories_per_mole" id="n2" rmin_half="1.4870 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" scale15="1" cutoff="9.0 * angstrom" pme_tolerance="0.00001"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
"""

xml_missing_torsion = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[#99:1]-[#99X4:2]-[#99:3]-[#99:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
  </ProperTorsions>
</SMIRNOFF>
"""

xml_ff_w_cosmetic_elements = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <!-- SMIRNOFF (SMIRKS Native Open Force Field) template file -->
  <Date>MMXVIII-VII-XIV</Date>
  <Author>Alice and Bob</Author>
  <!-- This file is meant for processing via openff.toolkit.typing.engines.smirnoff -->
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the bond energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Bonds version="0.3">
    <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" k="620.0 * kilocalories_per_mole/angstrom**2" length="1.526 * angstrom" parameters="k, length" parameterize_eval="blah=blah2"/>
    <Bond smirks="[#6X4:1]-[#6X3:2]" id="b2" k="634.0 * kilocalories_per_mole/angstrom**2" length="1.51 * angstrom"/>
  </Bonds>
  <!-- WARNING: AMBER functional forms drop the factor of 2 in the angle energy term, so cross-comparing this file with a corresponding .frcmod file, it will appear that the values here are twice as large as they should be. -->
  <Angles version="0.3" cosmetic_element="why not?">
    <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="109.5 * degree" id="a1" k="100.0 * kilocalories_per_mole/radian**2"/>
    <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.5 * degree" id="a2" k="70.0 * kilocalories_per_mole/radian**2"/>
  </Angles>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
    <Proper smirks="[*:1]-[#6X4:2]-[#6X4:3]-[*:4]" id="t1" idivf1="1" k1="0.156 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree"/>
    <Proper smirks="[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]" id="t2" idivf1="1" k1="0.180 * kilocalories_per_mole" periodicity1="3" phase1="0.0 * degree" periodicity2="2" phase2="180.0 * degree" idivf2="1" k2="0.250 * kilocalories_per_mole" periodicity3="1" phase3="180.0 * degree" idivf3="1" k3="0.200 * kilocalories_per_mole"/>
    <Proper smirks="[*:1]:[#6X4:2]~[#6X4:3]:[*:4]" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
  <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i1" k1="1.1 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
    <Improper smirks="[*:1]~[#6X3:2](~[#8X1:3])~[#8:4]" id="i2" k1="10.5 * kilocalories_per_mole" periodicity1="2" phase1="180. * degree"/>
  </ImproperTorsions>
  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
    <Atom smirks="[#1:1]" epsilon="0.0157 * kilocalories_per_mole" id="n1" rmin_half="0.6000 * angstrom"/>
    <Atom smirks="[#1:1]-[#6X4]" epsilon="0.0157 * kilocalories_per_mole" id="n2" rmin_half="1.4870 * angstrom"/>
  </vdW>
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" scale15="1" cutoff="9.0 * angstrom" pme_tolerance="0.00001"/>
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
"""

xml_toolkitam1bcc_ff = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ToolkitAM1BCC version="0.3"/>
</SMIRNOFF>
"""

xml_ethanol_library_charges_ff = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]-[#6:2](-[#1:3])(-[#1:4])-[#6:5](-[#1:6])(-[#1:7])-[#8:8]-[#1:9]" charge1="-0.02*elementary_charge" charge2="-0.2*elementary_charge" charge3="-0.02*elementary_charge" charge4="-0.02*elementary_charge" charge5="-0.1*elementary_charge" charge6="-0.01*elementary_charge" charge7="-0.01*elementary_charge" charge8="0.3*elementary_charge" charge9="0.08*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
"""

xml_ethanol_library_charges_in_parts_ff = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <!-- Note that the oxygen is covered twice here. The correct behavior should be to take the charge from the SECOND LibraryCharge, as it should overwrite the first -->
       <LibraryCharge smirks="[#1:1]-[#6:2](-[#1:3])(-[#1:4])-[#6:5](-[#1:6])(-[#1:7])-[#8:8]" charge1="-0.02*elementary_charge" charge2="-0.2*elementary_charge" charge3="-0.02*elementary_charge" charge4="-0.02*elementary_charge" charge5="-0.1*elementary_charge" charge6="-0.01*elementary_charge" charge7="-0.01*elementary_charge" charge8="-999*elementary_charge" />
       <LibraryCharge smirks="[#8:1]-[#1:2]" charge1="0.3*elementary_charge" charge2="0.08*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
"""

xml_ethanol_library_charges_by_atom_ff = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]-[#6]" charge1="-0.02*elementary_charge" />
       <LibraryCharge smirks="[#6X4:1]" charge1="-0.2*elementary_charge" />
       <LibraryCharge smirks="[#1:1]-[#6]-[#8]" charge1="-0.01*elementary_charge" />
       <LibraryCharge smirks="[#6X4:1]-[#8]" charge1="-0.1*elementary_charge" />
       <LibraryCharge smirks="[#8X2:1]" charge1="0.3*elementary_charge" />
       <LibraryCharge smirks="[#1:1]-[#8]" charge1="0.08*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
"""

xml_OH_library_charges_xml = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]" charge1="1.*elementary_charge" />
       <LibraryCharge smirks="[#8:1]" charge1="-2.*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
"""

xml_CH_zeroes_library_charges_xml = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge smirks="[#1:1]" charge1="0.*elementary_charge" />
       <LibraryCharge smirks="[#6:1]" charge1="0.*elementary_charge" />
    </LibraryCharges>
</SMIRNOFF>
"""

xml_spec_docs_ala_library_charges_xml = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge name="ALA" smirks="[NX3:1]([#1:2])([#6])[#6H1:3]([#1:4])([#6:5]([#1:6])([#1:7])[#1:8])[#6:9](=[#8:10])[#7]" charge1="-0.4157*elementary_charge" charge2="0.2719*elementary_charge" charge3="0.0337*elementary_charge" charge4="0.0823*elementary_charge" charge5="-0.1825*elementary_charge" charge6="0.0603*elementary_charge" charge7="0.0603*elementary_charge" charge8="0.0603*elementary_charge" charge9="0.5973*elementary_charge" charge10="-0.5679*elementary_charge"/>
    </LibraryCharges>
</SMIRNOFF>
"""

xml_spec_docs_tip3p_library_charges_xml = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
       <LibraryCharge name="TIP3P" smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]" charge1="0.417*elementary_charge" charge2="-0.834*elementary_charge" charge3="0.417*elementary_charge"/>
    </LibraryCharges>
</SMIRNOFF>
"""

xml_spec_docs_charge_increment_model_xml = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ChargeIncrementModel version="0.4" number_of_conformers="1" partial_charge_method="AM1-Mulliken">
    <!-- A fractional charge can be moved along a single bond -->
    <ChargeIncrement smirks="[#6X4:1]-[#6X3a:2]" charge_increment1="-0.0073*elementary_charge" charge_increment2="0.0073*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1]-[#6X3a:2]-[#7]" charge_increment1="0.0943*elementary_charge" charge_increment2="-0.0943*elementary_charge"/>
    <!--- Alternatively, fractional charges can be redistributed among any number of bonded atoms -->
    <ChargeIncrement smirks="[N:1]([H:2])([H:3])" charge_increment1="0.02*elementary_charge" charge_increment2="-0.01*elementary_charge" charge_increment3="-0.01*elementary_charge"/>
    <!-- As of version 0.4 of the ChargeIncrementModel tag, it is possible to define one less charge_increment attribute than there are tagged atoms -->
    <!-- The final, undefined charge_increment will be calculated as to make the sum of the charge_increments equal 0 -->
    <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge_increment1="-0.0718*elementary_charge"/>
    <ChargeIncrement smirks="[N]-[C:1]-[C:2]-[Cl:3]" charge_increment1="-0.123*elementary_charge" charge_increment2="0.456*elementary_charge" />
  </ChargeIncrementModel>
</SMIRNOFF>
"""

xml_charge_increment_model_formal_charges = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ChargeIncrementModel version="0.3" number_of_conformers="0" partial_charge_method="formal_charge"/>
</SMIRNOFF>
"""

xml_ff_bo = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Bonds version="0.3" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
    <Bond smirks="[#6X4:1]~[#8X2:2]" id="bbo1" k_bondorder1="101.0 * kilocalories_per_mole/angstrom**2" k_bondorder2="123.0 * kilocalories_per_mole/angstrom**2" length_bondorder1="1.4 * angstrom" length_bondorder2="1.3 * angstrom"/>
  </Bonds>
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]~[#6X3:2]~[#6X3:3]~[*:4]" id="tbo1" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
    <Proper smirks="[*:1]~[#6X4:2]~[#8X2:3]~[*:4]" id="tbo2" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
</SMIRNOFF>
"""

xml_ff_torsion_bo_standard_supersede = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
    <Proper smirks="[*:1]~[#6X3:2]~[#6X3:3]~[*:4]" id="tbo1" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
    <Proper smirks="[*:1]~[#6X4:2]~[#8X2:3]~[*:4]" id="tbo2" periodicity1="2" phase1="0.0 * degree" k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
    <Proper smirks="[*:1]~[#6X4:2]-[#8X2:3]~[#1:4]" id="t1" periodicity1="2" phase1="0.0 * degree" k1="1.20*kilocalories_per_mole" idivf1="1.0"/>
  </ProperTorsions>
</SMIRNOFF>
"""

xml_tip5p = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
            <LibraryCharge name="tip5p" smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]" charge1="0.*elementary_charge" charge2="0.*elementary_charge" charge3="0.*elementary_charge"/>
    </LibraryCharges>
    <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" switch_width="1.0 * angstrom" cutoff="9.0 * angstrom" method="cutoff">
            <Atom smirks="[#1:1]-[#8X2H2+0]-[#1]" epsilon="0. * mole**-1 * kilojoule" id="n35" sigma="1 * nanometer"/>
            <Atom smirks="[#1]-[#8X2H2+0:1]-[#1]" epsilon="0.66944 * mole**-1 * kilojoule" id="n35" sigma="0.312 * nanometer"/>
    </vdW>
     <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
        <Bond smirks="[#1:1]-[#8X2H2+0:2]-[#1]" length="0.9572 * angstrom" k="462750.4 * nanometer**-2 * mole**-1 * kilojoule" id="b1" />
    </Bonds>
    <Angles version="0.3" potential="harmonic">
        <Angle smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]" angle="1.82421813418 * radian" k="836.8 * mole**-1 * radian**-2 * kilojoule" id="a1" />
    </Angles>
    <VirtualSites version="0.3" exclusion_policy="parents">
        <VirtualSite
            type="DivalentLonePair"
            name="EP"
            smirks="[#1:2]-[#8X2H2+0:1]-[#1:3]"
            distance="0.70 * angstrom"
            charge_increment2="0.1205*elementary_charge"
            charge_increment1="0.0*elementary_charge"
            charge_increment3="0.1205*elementary_charge"
            sigma="1.0*angstrom"
            epsilon="0.0*kilocalories_per_mole"
            outOfPlaneAngle="54.71384225*degree"
            match="all_permutations" >
        </VirtualSite>
    </VirtualSites>
    <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" scale15="1.0" switch_width="0.0 * angstrom" cutoff="9.0 * angstrom"/>
  <Constraints version="0.3">
    <Constraint smirks="[#1:1]-[#8X2H2+0:2]-[#1]" id="c1" distance="0.9572 * angstrom"/>
    <Constraint smirks="[#1:1]-[#8X2H2+0]-[#1:2]" id="c2" distance="1.5139006545247014 * angstrom"/>
  </Constraints>
</SMIRNOFF>
"""

xml_gbsa_ff = """<?xml version='1.0' encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <GBSA version="0.3" gb_model="HCT" solvent_dielectric="78.5" solute_dielectric="1" sa_model="None" surface_area_penalty="5.4*calories/mole/angstroms**2" solvent_radius="1.4*angstroms">
          <Atom smirks="[*:1]" radius="0.15*nanometer" scale="0.8"/>
    </GBSA>
</SMIRNOFF>
"""

xml_charge_increment_model_ff_no_missing_cis = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge_increment1="-0.06*elementary_charge" charge_increment2="0.06*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1]-[#1:2]" charge_increment1="-0.01*elementary_charge" charge_increment2="0.01*elementary_charge"/>
    <ChargeIncrement smirks="[C:1][C:2][O:3]" charge_increment1="0.2*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_one_less_ci = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge_increment1="-0.06*elementary_charge" charge_increment2="0.06*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1]-[#1:2]" charge_increment1="-0.01*elementary_charge"/>
    <ChargeIncrement smirks="[C:1][C:2][O:3]" charge_increment1="0.2*elementary_charge" charge_increment2="-0.1*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_ethanol = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]-[#8:2]" charge_increment1="-0.05*elementary_charge" charge_increment2="0.05*elementary_charge"/>
    <ChargeIncrement smirks="[C:1][C:2][O:3]" charge_increment1="0.2*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_net_charge = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X3:1]-[#8X1-1:2]" charge_increment1="-0.05*elementary_charge" charge_increment2="0.05*elementary_charge"/>
    <ChargeIncrement smirks="[#6X3:1]=[#8X1:2]" charge_increment1="0.2*elementary_charge" charge_increment2="-0.2*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_override = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#1:1]-[#6:2]([#1:3])([#1:4])" charge_increment1="0.123*elementary_charge" charge_increment2="0.369*elementary_charge" charge_increment3="-0.123*elementary_charge" charge_increment4="0.123*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])([#1:4])" charge_increment1="0.3*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge" charge_increment4="-0.1*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_both_apply = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="0" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])([#1:4])" charge_increment1="0.3*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge" charge_increment4="-0.1*elementary_charge"/>
    <ChargeIncrement smirks="[#6X4:1][#6X4:2][#8]" charge_increment1="0.05*elementary_charge" charge_increment2="-0.05*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_match_once = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]([#1:2])[#6][#8]" charge_increment1="0.1*elementary_charge" charge_increment2="-0.1*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_match_two = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])[#6][#8]" charge_increment1="0.1*elementary_charge" charge_increment2="-0.05*elementary_charge" charge_increment3="-0.05*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_charge_increment_model_ff_match_all = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Electrostatics version="0.3" method="PME" scale12="0.0" scale13="0.0" scale14="0.833333" cutoff="9.0 * angstrom"/>
  <ChargeIncrementModel version="0.3" number_of_conformers="1" partial_charge_method="formal_charge">
    <ChargeIncrement smirks="[#6X4:1]([#1:2])([#1:3])([#1:4])" charge_increment1="0.3*elementary_charge" charge_increment2="-0.1*elementary_charge" charge_increment3="-0.1*elementary_charge" charge_increment4="-0.1*elementary_charge"/>
  </ChargeIncrementModel>
</SMIRNOFF>"""

xml_ff_virtual_sites_bondcharge_match_once = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[*:1]~[*:2]"
            distance="0.1*angstrom"
            charge_increment1="0.1*elementary_charge"
            charge_increment2="0.1*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.1*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[#7:1]~[#7:2]"
            distance="0.2*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.2*elementary_charge"
            sigma="0.2*angstrom"
            epsilon="0.2*kilocalories_per_mole"
            match="all_permutations" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[#7:1]~[#7:2]"
            distance="0.2*nanometers"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.2*elementary_charge"
            sigma="0.2*angstrom"
            epsilon="0.2*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

xml_ff_virtual_sites_bondcharge_match_all = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[*:1]~[*:2]"
            distance="0.1*angstrom"
            charge_increment1="0.1*elementary_charge"
            charge_increment2="0.1*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.1*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[#7:1]~[#7:2]"
            distance="0.2*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.2*elementary_charge"
            sigma="0.2*angstrom"
            epsilon="0.2*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP"
            smirks="[#7:1]~[#7:2]"
            distance="0.2*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.2*elementary_charge"
            sigma="0.2*angstrom"
            epsilon="0.2*kilocalories_per_mole"
            match="all_permutations" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

xml_ff_virtual_sites_bondcharge_match_once_two_names = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="BondCharge"
            name="EP1"
            smirks="[*:1]~[*:2]"
            distance="0.1*angstrom"
            charge_increment1="0.1*elementary_charge"
            charge_increment2="0.1*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.1*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
        <VirtualSite
            type="BondCharge"
            name="EP2"
            smirks="[*:1]~[*:2]"
            distance="0.2*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.2*elementary_charge"
            sigma="0.2*angstrom"
            epsilon="0.2*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

xml_ff_virtual_sites_monovalent = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="MonovalentLonePair"
            name="EP"
            smirks="[#8:1]~[#6:2]~[#6:3]"
            distance="0.1*angstrom"
            charge_increment1="0.1*elementary_charge"
            charge_increment2="0.1*elementary_charge"
            charge_increment3="0.1*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.1*kilocalories_per_mole"
            inPlaneAngle="110.*degree"
            outOfPlaneAngle="41*degree"
            match="all_permutations" >
        </VirtualSite>
        <VirtualSite
            type="MonovalentLonePair"
            name="EP"
            smirks="[#8:1]=[#6:2]-[#6:3]"
            distance="0.2*angstrom"
            charge_increment1="0.2*elementary_charge"
            charge_increment2="0.2*elementary_charge"
            charge_increment3="0.2*elementary_charge"
            sigma="0.2*angstrom"
            epsilon="0.2*kilocalories_per_mole"
            inPlaneAngle="120.*degree"
            outOfPlaneAngle="42*degree"
            match="all_permutations" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

xml_ff_virtual_sites_divalent_match_all = """<?xml version="1.0" encoding="ASCII"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="DivalentLonePair"
            name="EP"
            smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]"
            distance="0.70 * angstrom"
            charge_increment1="0.241*elementary_charge"
            charge_increment2="0.0*elementary_charge"
            charge_increment3="0.241*elementary_charge"
            sigma="3.12*angstrom"
            epsilon="0.16*kilocalories_per_mole"
            outOfPlaneAngle="54.71384225*degree"
            match="all_permutations" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

xml_ff_virtual_sites_trivalent_match_once = """<?xml version="1.0" encoding="ASCII"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="TrivalentLonePair"
            name="EP"
            smirks="[*:1]-[#7X3:2](-[*:3])-[*:4]"
            distance="0.50 * angstrom"
            charge_increment1="0.0*elementary_charge"
            charge_increment2="1.0*elementary_charge"
            charge_increment3="0.0*elementary_charge"
            charge_increment4="0.0*elementary_charge"
            sigma="0.0*angstrom"
            epsilon="0.0*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""

xml_ff_virtual_sites_trivalent_match_all = """
<?xml version="1.0" encoding="ASCII"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
      <Bond smirks="[*:1]~[*:2]" id="b999" k="500.0 * kilocalories_per_mole/angstrom**2" length="1.1 * angstrom"/>
    </Bonds>
    <VirtualSites version="0.3">
        <VirtualSite
            type="TrivalentLonePair"
            name="EP"
            smirks="[*:1]-[#7X3:2]-([*:3])-[*:4]"
            distance="0.70 * angstrom"
            charge_increment1="0.1*elementary_charge"
            charge_increment2="0.1*elementary_charge"
            charge_increment3="0.1*elementary_charge"
            charge_increment4="0.1*elementary_charge"
            sigma="0.1*angstrom"
            epsilon="0.1*kilocalories_per_mole"
            match="once" >
        </VirtualSite>
    </VirtualSites>
</SMIRNOFF>
"""


def round_charge(xml):
    """Round charge fields in a serialized OpenMM system to 2 decimal places"""
    # Example Particle line:                <Particle eps=".4577296" q="-.09709000587463379" sig=".1908"/>
    xmlsp = xml.split(' q="')
    for index, chunk in enumerate(xmlsp):
        # Skip file before first q=
        if index == 0:
            continue
        chunksp = chunk.split('" sig')
        chunksp[0] = f"{float(chunksp[0]):.2f}"
        chunk = '" sig'.join(chunksp)
        xmlsp[index] = chunk
    return ' q="'.join(xmlsp)


def generate_monatomic_ions():
    return (
        ("Li+", +1 * openmm_unit.elementary_charge),
        ("Na+", +1 * openmm_unit.elementary_charge),
        ("K+", +1 * openmm_unit.elementary_charge),
        ("Rb+", +1 * openmm_unit.elementary_charge),
        ("Cs+", +1 * openmm_unit.elementary_charge),
        ("F-", -1 * openmm_unit.elementary_charge),
        ("Cl-", -1 * openmm_unit.elementary_charge),
        ("Br-", -1 * openmm_unit.elementary_charge),
        ("I-", -1 * openmm_unit.elementary_charge),
    )


nonbonded_resolution_matrix = [
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic_potential": "Coulomb",
        "has_periodic_box": True,
        "omm_force": None,
        "exception": SMIRNOFFSpecUnimplementedError,
        "exception_match": "",
    },
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic_potential": "Coulomb",
        "has_periodic_box": False,
        "omm_force": openmm.NonbondedForce.NoCutoff,
        "exception": None,
        "exception_match": "",
    },
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic_potential": "reaction-field",
        "has_periodic_box": True,
        "omm_force": None,
        "exception": SMIRNOFFSpecUnimplementedError,
        "exception_match": "reaction-field",
    },
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic_potential": "reaction-field",
        "has_periodic_box": False,
        "omm_force": openmm.NonbondedForce.NoCutoff,
        "exception": None,
        "exception_match": "reaction-field",
    },
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic_potential": "PME",
        "has_periodic_box": True,
        "omm_force": openmm.NonbondedForce.PME,
        "exception": None,
        "exception_match": "",
    },
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic_potential": "PME",
        "has_periodic_box": False,
        "omm_force": openmm.NonbondedForce.NoCutoff,
        "exception": None,
        "exception_match": "",
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic_potential": "Coulomb",
        "has_periodic_box": True,
        "omm_force": None,
        "exception": IncompatibleParameterError,
        "exception_match": "",
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic_potential": "Coulomb",
        "has_periodic_box": False,
        "omm_force": None,
        "exception": SMIRNOFFSpecError,
        "exception_match": "vdw method PME.* is only valid for periodic systems",
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic_potential": "reaction-field",
        "has_periodic_box": True,
        "omm_force": None,
        "exception": IncompatibleParameterError,
        "exception_match": "must also be PME",
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic_potential": "reaction-field",
        "has_periodic_box": False,
        "omm_force": None,
        "exception": SMIRNOFFSpecError,
        "exception_match": "vdw method PME.* is only valid for periodic systems",
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic_potential": "PME",
        "has_periodic_box": True,
        "omm_force": openmm.NonbondedForce.LJPME,
        "exception": None,
        "exception_match": "",
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic_potential": "PME",
        "has_periodic_box": False,
        "omm_force": openmm.NonbondedForce.NoCutoff,
        "exception": SMIRNOFFSpecError,
        "exception_match": "vdw method PME.* is only valid for periodic systems",
    },
]

partial_charge_method_resolution_matrix = [
    {
        "toolkit": AmberToolsToolkitWrapper,
        "partial_charge_method": "AM1-Mulliken",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": AmberToolsToolkitWrapper,
        "partial_charge_method": "Gasteiger",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": AmberToolsToolkitWrapper,
        "partial_charge_method": "Madeup-ChargeMethod",
        "exception": ChargeMethodUnavailableError,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "AM1-Mulliken",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "Gasteiger",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "MMFF94",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "am1bcc",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "am1bccnosymspt",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "am1bccelf10",
        "exception": None,
        "exception_match": "",
    },
    {
        "toolkit": OpenEyeToolkitWrapper,
        "partial_charge_method": "Madeup-ChargeMethod",
        "exception": ChargeMethodUnavailableError,
        "exception_match": "",
    },
]


toolkit_registries = []
if OpenEyeToolkitWrapper.is_available():
    toolkit_registries.append(
        ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper]),
    )
if RDKitToolkitWrapper.is_available() and AmberToolsToolkitWrapper.is_available():
    toolkit_registries.append(
        ToolkitRegistry(
            toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        )
    )


class _ForceFieldFixtures:
    @pytest.fixture()
    def alkethoh_forcefield(self):
        return ForceField(
            get_data_file_path("test_forcefields/Frosst_AlkEthOH_parmAtFrosst.offxml")
        )

    @pytest.fixture()
    def force_field(self):
        return ForceField(get_data_file_path("test_forcefields/test_forcefield.offxml"))


class TestForceField(_ForceFieldFixtures):
    """Test the ForceField class"""

    def test_get_available_force_fields(self):
        """Ensure get_available_force_fields returns some expected data"""
        available_force_fields = get_available_force_fields(full_paths=False)

        # Incomplete list of some expected force fields
        expected_force_fields = [
            "smirnoff99Frosst-1.0.0.offxml",
            "smirnoff99Frosst-1.1.0.offxml",
            "openff-1.0.0.offxml",
            "openff_unconstrained-1.0.0.offxml",
            "openff-1.3.0.offxml",
            "openff-2.0.0.offxml",
            "ff14sb_off_impropers_0.0.3.offxml",
        ]

        for ff in expected_force_fields:
            assert ff in available_force_fields

    @pytest.mark.parametrize("full_path", [(True, False)])
    @pytest.mark.parametrize("force_field_file", [*get_available_force_fields()])
    def test_get_available_force_fields_loadable(self, full_path, force_field_file):
        """Ensure get_available_force_fields returns load-able files"""
        if "ff14sb" in force_field_file and "off_imp" not in force_field_file:
            pytest.skip(
                "Only the variants of the ff14SB port with SMIRNOFF-style impropers "
                "can be loaded by the toolkit by default. Those with Amber-style impropers"
                "require AmberImproperTorsionHandler."
            )
        ForceField(force_field_file)

    def test_force_field_case(self):
        """Ensure forcefield paths are loaded in a case-insensitive manner"""
        default_case = ForceField("smirnoff99Frosst-1.1.0.offxml")
        lower_case = ForceField("smirnoff99frosst-1.1.0.offxml")

        assert hash(default_case) == hash(lower_case)

    def test_do_not_load_in_child_dir(self, tmp_path):
        """Ensure force field XML files in nested subdirectories are not loaded
        when not explicitly pointed to."""
        nested_directory = tmp_path / os.path.join("a", "b", "c")
        os.makedirs(nested_directory, exist_ok=True)

        # Create a FF in a nested directory
        ForceField("openff-1.0.0.offxml").to_file(
            os.path.join(nested_directory, "force-field.offxml")
        )

        # Check that the file does not exist in the current working directory.
        assert not os.path.isfile("force-field.offxml")

        with pytest.raises(
            OSError, match="Source 'force-field.offxml' could not be read."
        ):
            ForceField("force-field.offxml")

    def test_load_bad_version(self):
        with pytest.raises(SMIRNOFFVersionError, match="99.3"):
            ForceField(
                get_data_file_path(
                    "test_forcefields/unsupported_smirnoff_version.offxml"
                ),
                disable_version_check=False,
            )

        ForceField(
            get_data_file_path("test_forcefields/unsupported_smirnoff_version.offxml"),
            disable_version_check=True,
        )

    def test_create_forcefield_no_args(self):
        """Test empty constructor"""
        forcefield = ForceField()

        # Should find BondHandler and AngleHandler, since they're default classes
        forcefield.get_parameter_handler("Bonds")
        forcefield.get_parameter_handler("Angles")

        # Shouldn't find InvalidKey handler, since it doesn't exist
        with pytest.raises(KeyError):
            forcefield.get_parameter_handler("InvalidKey")

        # Verify the aromatocitiy model is not None
        assert forcefield.aromaticity_model == "OEAroModel_MDL"

    def test_create_forcefield_aromaticity_model(self):
        """Test the aromaticiy_model argument of the constructor"""
        mdl = "OEAroModel_MDL"

        assert ForceField(aromaticity_model=mdl).aromaticity_model == mdl

        with pytest.raises(SMIRNOFFAromaticityError):
            ForceField(aromaticity_model="foobar")

    def test_create_forcefield_custom_handler_classes(self):
        """Test constructor given specific classes to register"""
        from openff.toolkit.typing.engines.smirnoff import BondHandler

        forcefield = ForceField(parameter_handler_classes=[BondHandler])

        # Should find BondHandler, since we registered it
        forcefield.get_parameter_handler("Bonds")

        # Shouldn't find AngleHandler, since we didn't allow that to be registered
        with pytest.raises(KeyError):
            forcefield.get_parameter_handler("Angles")

    def test_create_forcefield_from_file(self, force_field):
        """Test basic file loading in constructor"""
        assert len(force_field._parameter_handlers["Bonds"]._parameters) == 87
        assert len(force_field._parameter_handlers["Angles"]._parameters) == 38
        assert len(force_field._parameter_handlers["ProperTorsions"]._parameters) == 158
        assert len(force_field._parameter_handlers["ImproperTorsions"]._parameters) == 4
        assert len(force_field._parameter_handlers["vdW"]._parameters) == 35
        assert force_field.aromaticity_model == "OEAroModel_MDL"

    def test_load_bad_string(self):
        with pytest.raises(IOError) as exception_info:
            ForceField("1234")

        # This may need to be updated if the `openforcefields` package changes name;
        # searching through `site-packages/` and other paths is probably unrelaiable
        for match in [
            "Source '1234' could not be read",
            "SMIRNOFFParseError",
            "syntax error",
            "openforcefields",
        ]:
            assert match in str(exception_info.value)

    def test_load_bad_bytes(self):
        with pytest.raises(IOError) as exception_info:
            ForceField(b"the holy grail of computational chemistry")

        for match in [
            "Source 'b'the holy grail",
            "could not be read",
            "SMIRNOFFParseError",
            "syntax error",
        ]:
            assert match in str(exception_info.value)

        # See note in test_load_bad_string
        assert "openforcefields" not in str(exception_info.value)

    def test_load_filelike_object(self):
        with open(get_data_file_path("test_forcefields/test_forcefield.offxml")) as f:
            ForceField(f)

    def test_load_bad_filelike_object(self):
        with open(get_data_file_path("test_forcefields/mangled.offxml")) as f:
            with pytest.raises(IOError) as exception_info:
                ForceField(f)

            assert "while trying to parse source as an object" in str(
                exception_info.value
            )

    def test_create_forcefield_from_xml_string(self):
        forcefield = ForceField(xml_simple_ff)
        assert len(forcefield._parameter_handlers["Bonds"]._parameters) == 3
        assert len(forcefield._parameter_handlers["Angles"]._parameters) == 2
        assert len(forcefield._parameter_handlers["ProperTorsions"]._parameters) == 3
        assert len(forcefield._parameter_handlers["ImproperTorsions"]._parameters) == 2
        assert len(forcefield._parameter_handlers["vdW"]._parameters) == 2

    def test_load_do_not_convert_non_quantity_strings(self):
        """Reproduce part of #1493"""
        sage = ForceField("openff-2.0.0.offxml")

        for parameter_handler_name in sage.registered_parameter_handlers:
            parameter_handler = sage.get_parameter_handler(parameter_handler_name)

            for parameter in parameter_handler.parameters:
                assert isinstance(parameter.smirks, str)
                assert not isinstance(parameter.smirks, unit.Quantity)

                # Ensure that, for example, F isn't converted to Farad
                if (
                    parameter_handler_name == "LibraryCharges"
                    and parameter.name is not None
                ):
                    assert isinstance(parameter.name, str)
                    assert not isinstance(parameter.name, unit.Quantity)

    def test_pickle(self):
        """
        Test pickling and unpickling a force field
        """
        import pickle

        forcefield_1 = ForceField(xml_simple_ff)
        pickled = pickle.dumps(forcefield_1)
        forcefield_2 = pickle.loads(pickled)
        assert forcefield_1.to_string() == forcefield_2.to_string()

    def test_pickle_with_cosmetic_attributes(self):
        """
        Test pickling and unpickling a force field with cosmetic attributes
        """
        import pickle

        forcefield_1 = ForceField(
            xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True
        )
        pickled = pickle.dumps(forcefield_1)
        forcefield_2 = pickle.loads(pickled)
        assert forcefield_1.to_string() == forcefield_2.to_string()
        # Ensure that the cosmetic attributes stuck around
        assert "blah=blah2" in forcefield_2.to_string()

    def test_xml_string_roundtrip(self):
        """
        Test writing a ForceField to an XML string
        """
        forcefield_1 = ForceField(xml_simple_ff)
        string_1 = forcefield_1.to_string("XML")
        # Ensure that we have spaces instead of tabs
        assert "    " in string_1
        assert "\t" not in string_1
        forcefield_2 = ForceField(string_1)
        string_2 = forcefield_2.to_string("XML")
        assert string_1 == string_2

    def test_xml_string_roundtrip_keep_cosmetic(self):
        """
        Test roundtripping a force field to an XML string with and without retaining cosmetic elements
        """
        # Ensure an exception is raised if we try to read the XML string with cosmetic attributes
        with pytest.raises(
            SMIRNOFFSpecError,
            match="Unexpected kwarg [(]parameters: k, length[)]  passed",
        ):
            ForceField(xml_ff_w_cosmetic_elements)

        # Create a force field from XML successfully, by explicitly permitting cosmetic attributes
        forcefield_1 = ForceField(
            xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True
        )

        # Convert the force field back to XML
        string_1 = forcefield_1.to_string("XML", discard_cosmetic_attributes=False)

        # Ensure that the new XML string has cosmetic attributes in it
        assert 'cosmetic_element="why not?"' in string_1
        assert 'parameterize_eval="blah=blah2"' in string_1
        with pytest.raises(
            SMIRNOFFSpecError,
            match="Unexpected kwarg [(]parameters: k, length[)]  passed",
        ):
            ForceField(string_1, allow_cosmetic_attributes=False)

        # Complete the forcefield_1 --> string --> forcefield_2 roundtrip
        forcefield_2 = ForceField(string_1, allow_cosmetic_attributes=True)

        # Ensure that the force field remains the same after the roundtrip
        string_2 = forcefield_2.to_string("XML", discard_cosmetic_attributes=False)
        assert string_1 == string_2

        # Discard the cosmetic attributes and ensure that the string is different
        string_3 = forcefield_2.to_string("XML", discard_cosmetic_attributes=True)
        assert string_1 != string_3
        # Ensure that the new XML string does NOT have cosmetic attributes in it
        assert 'cosmetic_element="why not?"' not in string_3
        assert 'parameterize_eval="blah=blah2"' not in string_3

    def test_read_0_1_smirnoff(self):
        """Test reading an 0.1 spec OFFXML file"""
        ForceField(
            get_data_file_path(
                "test_forcefields/smirnoff99Frosst_reference_0_1_spec.offxml"
            )
        )

    def test_read_0_1_smirff(self):
        """Test reading an 0.1 spec OFFXML file, enclosed by the legacy "SMIRFF" tag"""
        ForceField(
            get_data_file_path(
                "test_forcefields/smirff99Frosst_reference_0_1_spec.offxml"
            )
        )

    def test_read_0_2_smirnoff(self):
        """Test reading an 0.2 spec OFFXML file"""
        ForceField(
            get_data_file_path(
                "test_forcefields/smirnoff99Frosst_reference_0_2_spec.offxml"
            )
        )

    @pytest.mark.parametrize("file_path_extension", ["xml", "XML", "offxml", "OFFXML"])
    @pytest.mark.parametrize(
        "specified_format",
        [
            None,
            "xml",
            "XML",
            ".xml",
            ".XML",
            "offxml",
            "OFFXML",
            ".offxml",
            ".OFFXML",
            XMLParameterIOHandler(),
        ],
    )
    def test_xml_file_roundtrip(self, file_path_extension, specified_format):
        """
        Test roundtripping a ForceField to and from an XML file
        """
        # These files will be deleted once garbage collection runs (end of this function)
        iofile1 = NamedTemporaryFile(suffix="." + file_path_extension)
        iofile2 = NamedTemporaryFile(suffix="." + file_path_extension)
        forcefield_1 = ForceField(xml_simple_ff)
        forcefield_1.to_file(iofile1.name, io_format=specified_format)
        forcefield_2 = ForceField(iofile1.name)
        forcefield_2.to_file(iofile2.name, io_format=specified_format)
        assert open(iofile1.name).read() == open(iofile2.name).read()

    @pytest.mark.parametrize("file_path_extension", ["xml", "XML", "offxml", "OFFXML"])
    @pytest.mark.parametrize(
        "specified_format",
        [
            None,
            "xml",
            "XML",
            ".xml",
            ".XML",
            "offxml",
            "OFFXML",
            ".offxml",
            ".OFFXML",
            XMLParameterIOHandler(),
        ],
    )
    def test_xml_file_roundtrip_keep_cosmetic(
        self, file_path_extension, specified_format
    ):
        """
        Test roundtripping a force field to an XML file with and without retaining cosmetic elements
        """
        # These files will be deleted once garbage collection runs (end of this function)
        iofile1 = NamedTemporaryFile(suffix="." + file_path_extension)
        iofile2 = NamedTemporaryFile(suffix="." + file_path_extension)
        iofile3 = NamedTemporaryFile(suffix="." + file_path_extension)

        # Ensure an exception is raised if we try to read the XML string with cosmetic attributes
        with pytest.raises(
            SMIRNOFFSpecError,
            match="Unexpected kwarg [(]parameters: k, length[)]  passed",
        ):
            ForceField(xml_ff_w_cosmetic_elements)

        # Create a force field from XML successfully
        forcefield_1 = ForceField(
            xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True
        )

        # Convert the force field back to XML, keeping cosmetic attributes
        forcefield_1.to_file(
            iofile1.name, discard_cosmetic_attributes=False, io_format=specified_format
        )

        # Ensure that the new XML string has cosmetic attributes in it
        assert 'cosmetic_element="why not?"' in open(iofile1.name).read()
        assert 'parameterize_eval="blah=blah2"' in open(iofile1.name).read()
        with pytest.raises(
            SMIRNOFFSpecError,
            match="Unexpected kwarg [(]parameters: k, length[)]  passed",
        ):
            ForceField(iofile1.name, allow_cosmetic_attributes=False)

        # Complete the forcefield_1 --> file --> forcefield_2 roundtrip
        forcefield_2 = ForceField(iofile1.name, allow_cosmetic_attributes=True)

        # Ensure that the force field remains the same after the roundtrip
        forcefield_2.to_file(
            iofile2.name, discard_cosmetic_attributes=False, io_format=specified_format
        )
        assert open(iofile1.name).read() == open(iofile2.name).read()

        # Discard the cosmetic attributes and ensure that the string is different
        forcefield_2.to_file(
            iofile3.name, discard_cosmetic_attributes=True, io_format=specified_format
        )
        assert open(iofile1.name).read() != open(iofile3.name).read()

        # Ensure that the new XML string does NOT have cosmetic attributes in it
        assert 'cosmetic_element="why not?"' not in open(iofile3.name).read()
        assert 'parameterize_eval="blah=blah2"' not in open(iofile3.name).read()

    def test_load_section_without_section_version(self):
        """Ensure that a SMIRNOFFSpecError is raised if we try to load a SMIRNOFF section without a version.
        Section versions are a requirement added in the 0.3 spec."""
        with pytest.raises(
            SMIRNOFFSpecError,
            match="Missing version while trying to construct "
            "<class 'openff.toolkit.typing.engines."
            "smirnoff.parameters.ToolkitAM1BCCHandler'>.",
        ):
            ForceField(xml_without_section_version)

    def test_load_two_sources(self):
        """Test loading data from two SMIRNOFF data sources"""
        ff = ForceField(
            xml_simple_ff, xml_ff_w_cosmetic_elements, allow_cosmetic_attributes=True
        )
        assert len(ff.get_parameter_handler("Bonds").parameters) == 5

    def test_load_two_sources_authors_dates(self):
        """Test that authors and dates are handled properly"""
        ff = ForceField(
            xml_ff_w_cosmetic_elements,
            xml_ff_w_comments,
            allow_cosmetic_attributes=True,
        )
        xml_str = ff.to_string("XML")
        assert (
            "<Author>Alice and Bob AND C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, "
            "UC Irvine; D. L. Mobley, UC Irvine</Author>" in xml_str
        )
        assert "<Date>MMXVIII-VII-XIV AND 2018-07-14</Date>" in xml_str

        # Test property getters
        assert (
            "Alice and Bob AND C. I. Bayly, OpenEye/UC Irvine; C. C. Bannan, "
            "UC Irvine; D. L. Mobley, UC Irvine" == ff.author
        )
        assert "MMXVIII-VII-XIV AND 2018-07-14" == ff.date

        # Test property setters
        ff.author = "Me"
        ff.date = "yesteryear"
        xml_str = ff.to_string("XML")
        assert "<Author>Me</Author>" in xml_str
        assert "<Date>yesteryear</Date>" in xml_str

        # Unset both author and date and ensure they don't get written out.
        ff.author = None
        ff.date = None
        xml_str = ff.to_string("XML")
        assert "<Author>" not in xml_str
        assert "<Date>" not in xml_str

    def test_load_two_sources_incompatible_tags(self):
        """Test loading data from two SMIRNOFF data sources which have incompatible physics"""
        # Make an XML force field with a modified vdW 1-4 scaling factor
        nonstandard_xml_ff = xml_ff_w_comments.replace('scale14="0.5"', 'scale14="1.0"')
        with pytest.raises(
            IncompatibleParameterError,
            match="handler value: 0.5, incompatible value: 1.0",
        ):
            ForceField(xml_simple_ff, nonstandard_xml_ff)

    def test_gbsahandler_sa_model_none(self):
        """
        Ensure that string values of "None" are correctly interpreted in the GBSAHandler's sa_model field
        """
        ForceField(xml_gbsa_ff)

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_parameterize_ethanol(self, toolkit_registry, force_field):
        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        molecules = [create_ethanol()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        force_field.create_openmm_system(
            topology,
            toolkit_registry=toolkit_registry,
        )

    @pytest.fixture()
    def create_circular_handler_dependencies(self):
        from openff.toolkit.typing.engines.smirnoff.parameters import (
            AngleHandler,
            BondHandler,
            ConstraintHandler,
        )

        # Modify the BondHandler and AngleHandler classes to depend on the other one running first during
        # system parameterization. Unfortunately, I can't figure out how to do this just to these _instances_
        # of the Handlers, so I modify them at the class level, and then un-modify them at the end of the test.
        orig_bh_depends = copy.deepcopy(BondHandler._DEPENDENCIES)
        orig_ah_depends = copy.deepcopy(AngleHandler._DEPENDENCIES)
        BondHandler._DEPENDENCIES = [ConstraintHandler, AngleHandler]
        AngleHandler._DEPENDENCIES = [ConstraintHandler, BondHandler]

        # The tests run here. Regardless of outcome, the code after `yield` runs after the test completes
        yield

        # Return handler dependencies to their original states
        BondHandler._DEPENDENCIES = orig_bh_depends
        AngleHandler._DEPENDENCIES = orig_ah_depends

    def test_parameterize_ethanol_missing_torsion(self):
        from openff.toolkit.typing.engines.smirnoff.parameters import (
            UnassignedProperTorsionParameterException,
        )

        forcefield = ForceField(xml_missing_torsion)
        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        molecules = [create_ethanol()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        with pytest.raises(
            UnassignedProperTorsionParameterException,
            match="- Topology indices [(]5, 0, 1, 6[)]: "
            r"names and elements [(](H\d+)? H[)], [(](C\d+)? C[)], [(](C\d+)? C[)], [(](H\d+)? H[)],",
        ):
            forcefield.create_openmm_system(topology)

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_parameterize_1_cyclohexane_1_ethanol(
        self,
        toolkit_registry,
        force_field,
    ):
        """Test parameterizing a periodic system of two distinct molecules"""
        pdbfile = app.PDBFile(
            get_data_file_path("systems/test_systems/1_cyclohexane_1_ethanol.pdb")
        )
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [create_ethanol(), create_cyclohexane()]
        # molecules = [Molecule.from_file(get_data_file_path(name)) for name in ('molecules/ethanol.mol2',
        #                                                                      'molecules/cyclohexane.mol2')]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        force_field.create_openmm_system(topology)

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_parameterize_1_cyclohexane_1_ethanol_vacuum(
        self,
        toolkit_registry,
        force_field,
    ):
        """Test parametrizing a nonperiodic system of two distinct molecules"""
        pdbfile = app.PDBFile(
            get_data_file_path("systems/test_systems/1_cyclohexane_1_ethanol.pdb")
        )
        molecules = [create_ethanol(), create_cyclohexane()]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        topology.box_vectors = None

        force_field.create_openmm_system(topology)

    @pytest.mark.slow
    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    @pytest.mark.parametrize(
        "box",
        [
            "ethanol_water.pdb",
            "cyclohexane_water.pdb",
            "cyclohexane_ethanol_0.4_0.6.pdb",
            "propane_methane_butanol_0.2_0.3_0.5.pdb",
        ],
    )
    def test_parameterize_large_system(
        self,
        toolkit_registry,
        box,
        force_field,
    ):
        """Test parameterizing a large system of several distinct molecules.
        This test is very slow, so it is only run if the --runslow option is provided to pytest.
        """
        box_file_path = get_data_file_path(
            os.path.join("systems", "packmol_boxes", box)
        )
        pdbfile = app.PDBFile(box_file_path)
        mol_names = ["water", "cyclohexane", "ethanol", "propane", "methane", "butanol"]
        sdf_files = [
            get_data_file_path(os.path.join("systems", "monomers", name + ".sdf"))
            for name in mol_names
        ]
        molecules = [Molecule.from_file(sdf_file) for sdf_file in sdf_files]
        topology = Topology.from_openmm(
            pdbfile.topology,
            unique_molecules=molecules,
        )

        force_field.create_openmm_system(
            topology,
            toolkit_registry=toolkit_registry,
        )
        # TODO: Add check to ensure system energy is finite

    @requires_openeye
    def test_parameterize_ethanol_different_reference_ordering_openeye(
        self, force_field
    ):
        """
        Test parameterizing the same PDB, using reference mol2s that have different atom orderings.
        The results of both should be identical.
        """
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        # Load the unique molecules with one atom ordering
        molecules1 = [Molecule.from_file(get_data_file_path("molecules/ethanol.sdf"))]
        topology1 = Topology.from_openmm(
            pdbfile.topology,
            unique_molecules=molecules1,
        )

        omm_system1 = force_field.create_openmm_system(
            topology1,
            toolkit_registry=toolkit_registry,
        )
        # Load the unique molecules with a different atom ordering
        molecules2 = [
            Molecule.from_file(get_data_file_path("molecules/ethanol_reordered.sdf"))
        ]
        topology2 = Topology.from_openmm(
            pdbfile.topology,
            unique_molecules=molecules2,
        )
        omm_system2 = force_field.create_openmm_system(
            topology2,
            toolkit_registry=toolkit_registry,
        )

        serialized_1 = XmlSerializer.serialize(omm_system1)
        serialized_2 = XmlSerializer.serialize(omm_system2)

        serialized_1 = round_charge(serialized_1)
        serialized_2 = round_charge(serialized_2)

        assert serialized_1 == serialized_2

    @requires_rdkit
    def test_parameterize_ethanol_different_reference_ordering_rdkit(self, force_field):
        """
        Test parameterizing the same PDB, using reference sdfs that have different atom orderings.
        The results of both should be identical.
        """
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        )
        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))

        # Load the unique molecules with one atom ordering
        molecules1 = [Molecule.from_file(get_data_file_path("molecules/ethanol.sdf"))]
        topology1 = Topology.from_openmm(
            pdbfile.topology,
            unique_molecules=molecules1,
        )

        omm_system1 = force_field.create_openmm_system(
            topology1,
            toolkit_registry=toolkit_registry,
        )

        # Load the unique molecules with a different atom ordering
        molecules2 = [
            Molecule.from_file(get_data_file_path("molecules/ethanol_reordered.sdf"))
        ]
        topology2 = Topology.from_openmm(
            pdbfile.topology,
            unique_molecules=molecules2,
        )
        omm_system2 = force_field.create_openmm_system(
            topology2,
            toolkit_registry=toolkit_registry,
        )

        serialized_1 = XmlSerializer.serialize(omm_system1)
        serialized_2 = XmlSerializer.serialize(omm_system2)

        serialized_1 = round_charge(serialized_1)
        serialized_2 = round_charge(serialized_2)

        assert serialized_1 == serialized_2

    @requires_rdkit
    def test_parameterize_mol_missing_stereo_rdkit(self, force_field):
        """
        Test parameterizing a molecule with undefined stereochemsity using the RDKit/AmberTools backend.
        """
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        )

        molecule = Molecule.from_smiles("CC1CCC(=O)O1", allow_undefined_stereo=True)
        topology = Topology.from_molecules([molecule])

        force_field.create_openmm_system(
            topology,
            toolkit_registry=toolkit_registry,
        )

    @requires_openeye
    def test_parameterize_mol_missing_stereo_openeye(self, force_field):
        """
        Test parameterizing a molecule with undefined stereochemsity using the OpenEye backend.
        """
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])

        molecule = Molecule.from_smiles("CC1CCC(=O)O1", allow_undefined_stereo=True)
        topology = Topology.from_molecules([molecule])

        force_field.create_openmm_system(
            topology,
            toolkit_registry=toolkit_registry,
        )

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_pass_invalid_kwarg_to_create_openmm_system(
        self, toolkit_registry, force_field
    ):
        """Test to ensure an exception is raised when an unrecognized kwarg is passed"""
        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        molecules = []
        molecules.append(Molecule.from_smiles("CCO"))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        with pytest.raises(
            TypeError,
            match="got an unexpected keyword argument .*invalid_kwarg.*",
        ):
            # TODO: specify desired toolkit_registry behavior in Interchange
            force_field.create_openmm_system(
                topology,
                invalid_kwarg="aaa",
                toolkit_registry=toolkit_registry,
            )

    def test_electrostatics_switch_width_unsupported(self):
        handler = ElectrostaticsHandler(version=0.4)
        with pytest.raises(
            SMIRNOFFSpecUnimplementedError,
            match="not support an electrostatic switch width",
        ):
            handler.switch_width = 1.234 * unit.nanometer

    @pytest.mark.parametrize("mod_cuoff", [True, False])
    def test_nonbonded_cutoff_no_box_vectors(self, mod_cuoff, force_field):
        """Ensure that the NonbondedForce objects use the cutoff specified in the
        ParameterHandler, not the OpenMM defaults"""
        top = Topology.from_molecules(create_ethanol())
        assert top.box_vectors is None

        if mod_cuoff:
            # Ensure a modified, non-default cutoff will be propogated through
            force_field["vdW"].cutoff = 0.777 * unit.nanometer
            force_field["Electrostatics"].cutoff = 0.777 * unit.nanometer

        omm_sys = force_field.create_openmm_system(top)

        for f in omm_sys.getForces():
            if isinstance(f, openmm.NonbondedForce):
                nonbonded_force = f

        found_cutoff = from_openmm(nonbonded_force.getCutoffDistance())
        vdw_cutoff = force_field["vdW"].cutoff
        e_cutoff = force_field["Electrostatics"].cutoff
        assert (found_cutoff - vdw_cutoff).m_as(unit.nanometer) < 1e-6
        assert (found_cutoff - e_cutoff).m_as(unit.nanometer) < 1e-6

    @pytest.mark.skip(
        reason="periodic_potential='Coulomb' not supported by Interchange"
    )
    def test_vdw_cutoff_overrides_electrostatics(self):
        topology = Molecule.from_smiles("[#18]").to_topology()
        topology.box_vectors = [3, 3, 3] * unit.nanometer

        force_field = ForceField()

        vdw_handler = vdWHandler(version=0.3)
        vdw_handler.method = "cutoff"
        vdw_handler.cutoff = 6.0 * unit.angstrom
        vdw_handler.scale14 = 1.0

        vdw_handler.add_parameter(
            {
                "smirks": "[#18:1]",
                "epsilon": 1.0 * unit.kilojoules_per_mole,
                "sigma": 1.0 * unit.angstrom,
            }
        )
        force_field.register_parameter_handler(vdw_handler)

        electrostatics_handler = ElectrostaticsHandler(version=0.3)
        electrostatics_handler.cutoff = 7.0 * unit.angstrom
        electrostatics_handler.periodic_potential = "PME"
        force_field.register_parameter_handler(electrostatics_handler)

        library_charges = LibraryChargeHandler(version=0.3)
        library_charges.add_parameter(
            {
                "smirks": "[#18:1]",
                "charge1": 0.0 * unit.elementary_charge,
            }
        )
        force_field.register_parameter_handler(library_charges)

        system = force_field.create_openmm_system(topology)

        assert np.isclose(
            system.getForce(0).getCutoffDistance().value_in_unit(openmm_unit.angstrom),
            6.0,
        )

        # Ensure an exception is raised when the electrostatics cutoff is meaningful
        # and mismatched
        force_field.deregister_parameter_handler(force_field["Electrostatics"])

        # TODO: Don't change the box vectors once this case is supported
        topology.box_vectors = None
        electrostatics_handler.periodic_potential = "Coulomb"
        force_field.register_parameter_handler(electrostatics_handler)

        with pytest.raises(IncompatibleParameterError, match="cutoff must equal"):
            force_field.create_openmm_system(topology)

    def test_nondefault_nonbonded_cutoff(self):
        """Test that the cutoff of the NonbondedForce is set properly when vdW and Electrostatics cutoffs
        are identical but not the psuedo-default value of 9.0 A."""
        topology = Molecule.from_smiles("[#18]").to_topology()
        topology.box_vectors = [3, 3, 3] * unit.nanometer

        force_field = ForceField()

        vdw_handler = vdWHandler(version=0.3)
        vdw_handler.method = "cutoff"
        vdw_handler.cutoff = 7.89 * unit.angstrom
        vdw_handler.scale14 = 1.0

        vdw_handler.add_parameter(
            {
                "smirks": "[#18:1]",
                "epsilon": 1.0 * unit.kilojoules_per_mole,
                "sigma": 1.0 * unit.angstrom,
            }
        )
        force_field.register_parameter_handler(vdw_handler)

        electrostatics_handler = ElectrostaticsHandler(version=0.3)
        electrostatics_handler.cutoff = 7.89 * unit.angstrom
        electrostatics_handler.periodic_potential = "PME"
        force_field.register_parameter_handler(electrostatics_handler)

        library_charges = LibraryChargeHandler(version=0.3)
        library_charges.add_parameter(
            {
                "smirks": "[#18:1]",
                "charge1": 0.0 * unit.elementary_charge,
            }
        )
        force_field.register_parameter_handler(library_charges)

        system = force_field.create_openmm_system(topology)

        found_cutoff = (
            system.getForce(0).getCutoffDistance().value_in_unit(openmm_unit.angstrom)
        )

        assert abs(found_cutoff - 7.89) < 1e-6

    def test_registered_parameter_handlers(self, force_field):
        """Test registered_parameter_handlers property"""
        registered_handlers = force_field.registered_parameter_handlers

        expected_handlers = [
            "Bonds",
            "Angles",
            "ProperTorsions",
            "ImproperTorsions",
            "vdW",
            "Electrostatics",
            "ToolkitAM1BCC",
        ]

        for expected_handler in expected_handlers:
            assert expected_handler in registered_handlers

        assert "LibraryChrages" not in registered_handlers

    def test_parameter_handler_lookup(self, force_field):
        """Ensure __getitem__ lookups work"""

        handlers_before = sorted(force_field._parameter_handlers)

        for val in handlers_before:
            looked_up_handler = force_field[val]
            assert isinstance(looked_up_handler, ParameterHandler)

        handlers_after = sorted(force_field._parameter_handlers)

        assert handlers_before == handlers_after

    @pytest.mark.parametrize("unregistered_handler", ["LibraryCharges", "foobar"])
    def test_unregistered_parameter_handler_lookup(
        self, unregistered_handler, force_field
    ):
        """Ensure ForceField.__getitem__ lookups do not register new handlers"""
        assert unregistered_handler not in force_field._parameter_handlers
        with pytest.raises(KeyError, match=unregistered_handler):
            force_field[unregistered_handler]
        assert unregistered_handler not in force_field._parameter_handlers

    def test_lookup_parameter_handler_object(self, force_field):
        """Ensure ForceField.__getitem__ raises NotImplemented when passed a ParameterHandler object"""
        bonds = force_field["Bonds"]
        with pytest.raises(NotImplementedError):
            force_field[bonds]
        with pytest.raises(NotImplementedError):
            force_field[type(bonds)]

    def test_lookup_parameter_type(self, force_field):
        """Test both ForceField and ParameterHandler __getitem__ methods"""
        smirks = "[#6X4:1]-[#6X3:2]=[#8X1+0]"

        param = force_field["Bonds"][smirks]
        assert param.smirks == smirks

        # Look up the same param by its index in the ParameterList
        param_idx = 2
        assert param == force_field["Bonds"][param_idx]

        with pytest.raises(
            ParameterLookupError,
            match="Lookup by instance is not supported",
        ):
            force_field["Bonds"][param]
        with pytest.raises(
            ParameterLookupError,
            match=r" not found in ParameterList",
        ):
            force_field["vdW"][smirks]

    @pytest.mark.parametrize(
        "to_deregister",
        [
            "ToolkitAM1BCC",
            ToolkitAM1BCCHandler,
            ToolkitAM1BCCHandler(skip_version_check=True),
        ],
    )
    def test_deregister_parameter_handler(self, to_deregister):
        """Ensure that ForceField.deregister_parameter_handler behaves correctly"""
        ff = ForceField(xml_simple_ff)
        # Make sure the handler is present in the test force field
        handler_found = False
        for handler in ff._parameter_handlers.values():
            if isinstance(handler, ToolkitAM1BCCHandler):
                handler_found = True
        assert handler_found
        ff.deregister_parameter_handler(to_deregister)
        # Make sure the handler is absent after deletion
        handler_found = False
        for handler in ff._parameter_handlers.values():
            if isinstance(handler, ToolkitAM1BCCHandler):
                handler_found = True
        assert not (handler_found)

        with pytest.raises(KeyError):
            ff.deregister_parameter_handler(to_deregister)

    def test_hash(self):
        """Test hashes on all available force fields"""
        ffs = get_available_force_fields()

        for ff1, ff2 in itertools.combinations(ffs, 2):
            assert hash(ff1) != hash(ff2)

    def test_hash_cosmetic(self, force_field):
        """Test that adding a cosmetic attribute does not change the hash"""
        hash_without_cosmetic = hash(force_field)

        force_field.get_parameter_handler("Bonds").add_cosmetic_attribute("foo", 4)
        assert "foo" in force_field["Bonds"]._cosmetic_attribs

        hash_with_cosmetic = hash(force_field)

        assert hash_with_cosmetic == hash_without_cosmetic

    def test_hash_strip_author_date(self):
        """Ensure that author and date are ignored in hashing"""
        ff_no_data = ForceField()
        ff_with_author_date = ForceField()
        ff_with_author_date.author = "John Doe"
        ff_with_author_date.date = "2020-01-01"

        assert hash(ff_no_data) == hash(ff_with_author_date)

    def test_hash_strip_ids(self):
        """Test the behavior of strip_ids arg to __hash__()"""
        from openff.toolkit.typing.engines.smirnoff import BondHandler

        length = 1 * unit.angstrom
        k = 10 * unit.kilocalorie_per_mole / unit.angstrom**2

        param_with_id = {
            "smirks": "[*:1]-[*:2]",
            "length": length,
            "k": k,
            "id": "b1",
        }

        param_without_id = {
            "smirks": "[*:1]-[*:2]",
            "length": length,
            "k": k,
        }

        ff_with_id = ForceField()
        ff_without_id = ForceField()

        ff_with_id.register_parameter_handler(BondHandler(version=0.3))
        ff_without_id.register_parameter_handler(BondHandler(version=0.3))

        ff_with_id.get_parameter_handler("Bonds").add_parameter(param_with_id)
        ff_without_id.get_parameter_handler("Bonds").add_parameter(param_without_id)

        assert hash(ff_with_id) == hash(ff_without_id)

    def test_issue_1216(self):
        """Test that the conflict between vsitehandler and librarychargehandler
        identified in issue #1216 remains resolved.

        https://github.com/openforcefield/openff-toolkit/issues/1216
        """
        from openff.toolkit.typing.engines.smirnoff.parameters import VirtualSiteHandler

        force_field = ForceField()
        force_field.get_parameter_handler("Electrostatics")

        vsite_handler: VirtualSiteHandler = force_field.get_parameter_handler(
            "VirtualSites"
        )
        vsite_handler.add_parameter(
            {
                "smirks": "[#6:1][#9:2]",
                "name": "EP",
                "type": "BondCharge",
                "distance": 1.0 * unit.angstrom,
                "match": "all_permutations",
                "charge_increment1": 0.2 * unit.elementary_charge,
                "charge_increment2": 0.1 * unit.elementary_charge,
                "sigma": 1.0 * unit.angstrom,
                "epsilon": 0.0 * unit.kilocalorie_per_mole,
            }
        )

        library_handler: LibraryChargeHandler = force_field.get_parameter_handler(
            "LibraryCharges"
        )
        library_handler.add_parameter(
            {
                "smirks": "[F:2][C:1]([H:3])([H:4])([H:5])",
                "charge": [
                    0.3 * unit.elementary_charge,
                    -0.15 * unit.elementary_charge,
                    -0.05 * unit.elementary_charge,
                    -0.05 * unit.elementary_charge,
                    -0.05 * unit.elementary_charge,
                ],
            }
        )
        force_field.label_molecules(Molecule.from_smiles("CF").to_topology())

    def test_issue_1475(self):
        """Reproduce issue #1475."""
        force_field = ForceField("openff-2.0.0.offxml")

        class BogusHandler(ParameterHandler):
            _TAGNAME = "bogus"

        bogus_handler = BogusHandler(version=0.3)

        force_field.register_parameter_handler(bogus_handler)

        assert "bogus" in force_field._parameter_handlers.keys()
        assert bogus_handler in force_field._parameter_handlers.values()

        assert "bogus" in force_field._parameter_handler_classes.keys()
        assert BogusHandler in force_field._parameter_handler_classes.values()

        assert force_field["bogus"] is not None


class TestForceFieldPluginLoading:
    def test_handlers_tracked_if_already_loaded(self):
        """Reproduce issue #1542."""
        from openff.toolkit.typing.engines.smirnoff.plugins import load_handler_plugins

        plugins = load_handler_plugins()

        assert (
            len(plugins) > 0
        ), "Test assumes that some ParameterHandler plugins are available"

        assert ForceField(load_plugins=False)._plugin_parameter_handler_classes == []
        assert ForceField(load_plugins=True)._plugin_parameter_handler_classes == [
            *plugins
        ]


class TestForceFieldSerializaiton(_ForceFieldFixtures):
    def test_json_dump(self, force_field):
        """Test that at a ForceField can at least be dumped to JSON without error."""
        import json

        json.dumps(force_field._to_smirnoff_data())


class TestForceFieldChargeAssignment(_ForceFieldFixtures):
    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_charges_from_molecule(self, toolkit_registry, force_field):
        """Test skipping charge generation and instead getting charges from the original Molecule"""
        # Create an ethanol molecule without using a toolkit
        molecules = [create_ethanol()]

        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        omm_system = force_field.create_openmm_system(
            topology, charge_from_molecules=molecules, toolkit_registry=toolkit_registry
        )
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = (
            (0, -0.4 * openmm_unit.elementary_charge),
            (1, -0.3 * openmm_unit.elementary_charge),
            (2, -0.2 * openmm_unit.elementary_charge),
        )
        for particle_index, expected_charge in expected_charges:
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_charges_from_molecule_reordered(self, toolkit_registry, force_field):
        """A copy of test_charges_from_molecule but with the same molecule in a different atom order."""
        molecules = [create_ethanol()]

        # In 1_ethanol_reordered.pdb, the first three atoms go O-C-C instead of C-C-O. This part of the test ensures
        # that the charges are correctly mapped according to this PDB in the resulting system.
        pdbfile = app.PDBFile(
            get_data_file_path("systems/test_systems/1_ethanol_reordered.pdb")
        )
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        omm_system = force_field.create_openmm_system(
            topology,
            charge_from_molecules=molecules,
            toolkit_registry=toolkit_registry,
        )
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = (
            (0, -0.2 * openmm_unit.elementary_charge),
            (1, -0.4 * openmm_unit.elementary_charge),
            (2, -0.3 * openmm_unit.elementary_charge),
        )
        for particle_index, expected_charge in expected_charges:
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_nonintegral_charge_exception(self, toolkit_registry, force_field):
        # Create an ethanol molecule without using a toolkit
        from openff.interchange.exceptions import NonIntegralMoleculeChargeException

        ethanol = create_ethanol()
        ethanol.partial_charges[0] = 1.0 * unit.elementary_charge

        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=[ethanol])

        with pytest.raises(
            NonIntegralMoleculeChargeException, match="Molecule .* has a net charge"
        ):
            force_field.create_openmm_system(
                topology,
                charge_from_molecules=[ethanol],
                toolkit_registry=toolkit_registry,
            )

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_nonintegral_charge_override(self, toolkit_registry, force_field):
        ethanol = create_ethanol()
        ethanol.partial_charges[0] = 1.0 * unit.elementary_charge

        pdbfile = app.PDBFile(get_data_file_path("systems/test_systems/1_ethanol.pdb"))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=[ethanol])

        # Pass when the `allow_nonintegral_charges` keyword is included
        force_field.create_openmm_system(
            topology,
            charge_from_molecules=[ethanol],
            toolkit_registry=toolkit_registry,
            allow_nonintegral_charges=True,
        )

    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_some_charges_from_molecule(self, toolkit_registry, force_field):
        """
        Test creating an OpenMM system where some charges come from a Molecule, but others come from toolkit
        calculation
        """
        ethanol = create_ethanol()
        cyclohexane = create_cyclohexane()
        molecules = [ethanol, cyclohexane]

        pdbfile = app.PDBFile(
            get_data_file_path("systems/test_systems/1_cyclohexane_1_ethanol.pdb")
        )
        topology = Topology.from_openmm(
            pdbfile.topology,
            unique_molecules=molecules,
        )

        omm_system = force_field.create_openmm_system(
            topology, charge_from_molecules=[ethanol], toolkit_registry=toolkit_registry
        )
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = (
            (18, -0.4 * openmm_unit.elementary_charge),
            (19, -0.3 * openmm_unit.elementary_charge),
            (20, -0.2 * openmm_unit.elementary_charge),
        )
        for particle_index, expected_charge in expected_charges:
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge
        for atom_index in range(topology.n_atoms):
            q, _, _ = nonbondedForce.getParticleParameters(atom_index)
            assert q != (0.0 * unit.elementary_charge)

    @pytest.mark.parametrize(
        "ff_inputs",
        [
            [
                get_data_file_path("test_forcefields/test_forcefield.offxml"),
                get_data_file_path("test_forcefields/tip3p.offxml"),
            ],
            ["openff-2.0.0.offxml"],
        ],
    )
    def test_library_charges_to_single_water(self, ff_inputs):
        """Test assigning charges to one water molecule using library charges"""
        ff = ForceField(*ff_inputs)
        mol = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "water.sdf"))
        )
        omm_system = ff.create_openmm_system(mol.to_topology())
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [-0.834, 0.417, 0.417] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_parse_library_charges_from_spec_docs(self):
        """Ensure that the examples for librarycharges in the SMIRNOFF spec page are still valid"""
        # TODO: This test is practically useless while the XML strings are hard-coded at the top of this file.
        #       We should implement something like doctests for the XML snippets on the SMIRNOFF spec page.
        ForceField(xml_spec_docs_ala_library_charges_xml)
        ForceField(xml_spec_docs_tip3p_library_charges_xml)

    def test_parse_charge_increment_model_from_spec_docs(self):
        """Ensure that the examples for librarycharges in the SMIRNOFF spec page are still valid"""
        # TODO: This test is practically useless while the XML strings are hard-coded at the top of this file.
        #       We should implement something like doctests for the XML snippets on the SMIRNOFF spec page.
        ForceField(xml_spec_docs_charge_increment_model_xml)

    def test_charge_increment_model_forward_and_reverse_ethanol(self):
        """Test application of ChargeIncrements to the same molecule with different orderings in the topology"""
        file_path = get_data_file_path("test_forcefields/test_forcefield.offxml")
        ff = ForceField(file_path, xml_charge_increment_model_ff_ethanol)
        del ff._parameter_handlers["ToolkitAM1BCC"]
        top = Topology.from_molecules([create_ethanol(), create_reversed_ethanol()])
        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0.2,
            -0.15,
            -0.05,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.05,
            -0.15,
            0.2,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

    def test_charge_increment_model_one_less_ci_than_tagged_atom(self):
        """
        Ensure that we support the behavior where a ChargeIncrement is initialized with one less chargeincrement value
        than tagged atom. We test this by making two equivalent (one with fully explicit CIs, the other with some
        implicit CIa) FFs and ensuring that both perform the same parameterization.
        """
        # Make a FF from each OFFXML string
        file_path = get_data_file_path("test_forcefields/test_forcefield.offxml")
        ff1 = ForceField(file_path, xml_charge_increment_model_ff_one_less_ci)
        del ff1._parameter_handlers["ToolkitAM1BCC"]
        ff2 = ForceField(file_path, xml_charge_increment_model_ff_no_missing_cis)
        del ff2._parameter_handlers["ToolkitAM1BCC"]
        top = Topology.from_molecules([create_ethanol()])
        # Make a system from each FF
        sys1 = ff2.create_openmm_system(top)
        sys2 = ff2.create_openmm_system(top)
        # Extract the nonbonded force from each system
        nonbonded_force1 = [
            force
            for force in sys1.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        nonbonded_force2 = [
            force
            for force in sys2.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]

        # Ensure that the systems both have the correct charges assigned
        expected_charges = [
            0.17,
            -0.18,
            -0.04,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge1, _, _ = nonbonded_force1.getParticleParameters(idx)
            charge2, _, _ = nonbonded_force2.getParticleParameters(idx)
            assert (
                abs(charge1 - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )
            assert charge1 == charge2

    def test_charge_increment_model_invalid_number_of_cis(self, force_field):
        """
        Ensure that we support the behavior where a ChargeIncrement with an incorrect number of tagged atoms
        and chargeincrementX values riases an error
        """
        del force_field._parameter_handlers["ToolkitAM1BCC"]
        cimh = force_field.get_parameter_handler(
            "ChargeIncrementModel",
            handler_kwargs={"version": "0.3", "partial_charge_method": "formal_charge"},
        )
        cimh.add_parameter(
            {
                "smirks": "[C:1][C:2][O:3]",
                "charge_increment1": 0.3 * unit.elementary_charge,
                "charge_increment2": -0.2 * unit.elementary_charge,
                "charge_increment3": -0.1 * unit.elementary_charge,
            }
        )

        # Add ONE MORE chargeincrement parameter than there are tagged atoms and ensure an exception is raised
        cimh.parameters[0].charge_increment.append(0.01 * unit.elementary_charge)
        top = Topology.from_molecules([create_ethanol()])
        with pytest.raises(
            SMIRNOFFSpecError,
            match="number of chargeincrements must be either the same",
        ):
            force_field.create_openmm_system(top)

        # Ensure that parameterization with the correct number of increments DOES NOT raise an exception
        cimh.parameters[0].charge_increment = cimh.parameters[0].charge_increment[:2]
        force_field.create_openmm_system(top)

        # Add TWO LESS chargeincrement parameters than there are tagged atoms and ensure an exception is raised
        cimh.parameters[0].charge_increment = cimh.parameters[0].charge_increment[:1]
        with pytest.raises(
            SMIRNOFFSpecError,
            match="number of chargeincrements must be either the same",
        ):
            force_field.create_openmm_system(top)

    def test_charge_increment_model_initialize_with_no_elements(self):
        """Ensure that we can initialize a ForceField object from an OFFXML with a ChargeIncrementModel header, but no
        ChargeIncrement elements"""
        ForceField(xml_charge_increment_model_formal_charges)

    def test_charge_increment_model_net_charge(self):
        """Test application of charge increments on a molecule with a net charge"""
        file_path = get_data_file_path("test_forcefields/test_forcefield.offxml")
        ff = ForceField(file_path, xml_charge_increment_model_ff_net_charge)
        del ff._parameter_handlers["ToolkitAM1BCC"]

        acetate = create_acetate()
        top = acetate.to_topology()

        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0,
            0.15,
            -0.2,
            -0.95,
            0,
            0,
            0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

    def test_charge_increment_model_deduplicate_symmetric_matches(self):
        """Test that chargeincrementmodelhandler deduplicates symmetric matches"""
        ethanol = create_ethanol()
        top = ethanol.to_topology()

        file_path = get_data_file_path("test_forcefields/test_forcefield.offxml")
        # Test a charge increment that matches all C-H bonds at once
        # (this should be applied once: C0-H3-H4-H5)
        ff = ForceField(file_path, xml_charge_increment_model_ff_match_all)
        del ff._parameter_handlers["ToolkitAM1BCC"]

        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0.3,
            0,
            0,
            -0.1,
            -0.1,
            -0.1,
            0.0,
            0.0,
            0.0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

        # Test a charge increment that matches two C-H bonds at a time
        # (this should be applied 3 times: C0-H3-H4, C0-H3-H5, C0-H4-H5)
        ff = ForceField(file_path, xml_charge_increment_model_ff_match_two)
        del ff._parameter_handlers["ToolkitAM1BCC"]

        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0.3,
            0,
            0,
            -0.1,
            -0.1,
            -0.1,
            0.0,
            0.0,
            0.0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

        # Test a charge increment that matches ONE C-H bond at a time
        # (this should be applied three times: C0-H3, C0-H4, C0-H5)
        ff = ForceField(file_path, xml_charge_increment_model_ff_match_once)
        del ff._parameter_handlers["ToolkitAM1BCC"]

        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0.3,
            0,
            0,
            -0.1,
            -0.1,
            -0.1,
            0.0,
            0.0,
            0.0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

    def test_charge_increment_model_completely_overlapping_matches_override(self):
        """Ensure that DIFFERENT chargeincrements override one another if they apply to the
        same atoms, regardless of order"""
        file_path = get_data_file_path("test_forcefields/test_forcefield.offxml")
        ff = ForceField(file_path, xml_charge_increment_model_ff_override)
        del ff._parameter_handlers["ToolkitAM1BCC"]

        ethanol = create_ethanol()
        top = ethanol.to_topology()
        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0.3,
            0,
            0,
            -0.1,
            -0.1,
            -0.1,
            0.0,
            0.0,
            0.0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

    def test_charge_increment_model_partially_overlapping_matches_both_apply(self):
        """Ensure that DIFFERENT chargeincrements BOTH get applied if they match
        a partially-overlapping set of atoms"""
        file_path = get_data_file_path("test_forcefields/test_forcefield.offxml")
        ff = ForceField(file_path, xml_charge_increment_model_ff_both_apply)
        del ff._parameter_handlers["ToolkitAM1BCC"]

        ethanol = create_ethanol()
        top = ethanol.to_topology()
        sys = ff.create_openmm_system(top)
        nonbonded_force = [
            force
            for force in sys.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]
        expected_charges = [
            0.35,
            -0.05,
            0,
            -0.1,
            -0.1,
            -0.1,
            0.0,
            0.0,
            0.0,
        ] * openmm_unit.elementary_charge
        for idx, expected_charge in enumerate(expected_charges):
            charge, _, _ = nonbonded_force.getParticleParameters(idx)
            assert (
                abs(charge - expected_charge) < 1.0e-6 * openmm_unit.elementary_charge
            )

    @pytest.mark.parametrize("inputs", partial_charge_method_resolution_matrix)
    def test_partial_charge_resolution(self, inputs):
        """Check that the proper partial charge methods are available, and that unavailable partial charge methods
        raise an exception.
        """
        toolkit_wrapper_class = inputs["toolkit"]
        if not (toolkit_wrapper_class.is_available()):
            pytest.skip(f"{toolkit_wrapper_class} is not available.")
        toolkit_wrapper = toolkit_wrapper_class()
        partial_charge_method = inputs["partial_charge_method"]
        expected_exception = inputs["exception"]
        expected_exception_match = inputs["exception_match"]
        ethanol = create_ethanol()
        ethanol.generate_conformers()
        if expected_exception is None:
            ethanol.assign_partial_charges(
                partial_charge_method=partial_charge_method,
                toolkit_registry=toolkit_wrapper,
            )
            abs_charge_sum = 0.0 * unit.elementary_charge

            # Ensure that nonzero charges were assigned
            for pc in ethanol.partial_charges:
                abs_charge_sum += abs(pc)
            assert abs_charge_sum > 0.5 * unit.elementary_charge

        else:
            with pytest.raises(expected_exception, match=expected_exception_match):
                ethanol.assign_partial_charges(
                    partial_charge_method=partial_charge_method,
                    toolkit_registry=toolkit_wrapper,
                )

    def test_library_charge_hierarchy(self):
        """Test assigning charges to one water molecule using library charges, where two LCs match and the
        assignment is determined by order they are added to the force field"""
        # Test with xml_OH_library_charges_xml loaded last, which should assign dummy partial charges
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            get_data_file_path("test_forcefields/tip3p.offxml"),
            xml_OH_library_charges_xml,
        )
        mol = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "water.sdf"))
        )
        omm_system = ff.create_openmm_system(mol.to_topology())
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [-2.0, 1.0, 1.0] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

        # Test again, but with tip3p.offxml loaded last (loading the correct partial charges)
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_OH_library_charges_xml,
            get_data_file_path("test_forcefields/tip3p.offxml"),
        )
        omm_system = ff.create_openmm_system(mol.to_topology())
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [-0.834, 0.417, 0.417] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_library_charges_to_two_waters(self):
        """Test assigning charges to two water molecules using library charges"""
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            get_data_file_path("test_forcefields/tip3p.offxml"),
        )
        mol = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "water.sdf"))
        )
        top = Topology.from_molecules([mol, mol])
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [
            -0.834,
            0.417,
            0.417,
            -0.834,
            0.417,
            0.417,
        ] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_library_charges_to_three_ethanols_different_atom_ordering(self):
        """Test assigning charges to three ethanols with different atom orderings"""
        # Define a library charge parameter for ethanol (C1-C2-O3) where C1 has charge -0.2, and its Hs have -0.02,
        # C2 has charge -0.1 and its Hs have -0.01, and O3 has charge 0.3, and its H has charge 0.08

        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ethanol_library_charges_ff,
        )

        # ethanol.sdf
        #      H5   H8
        #      |    |
        # H6 - C1 - C2 - O3 - H4
        #      |    |
        #      H7   H9
        #
        # ethanol_reordered.sdf (The middle C and O switch indices)
        #      H5   H8
        #      |    |
        # H6 - C1 - C3 - O2 - H4
        #      |    |
        #      H7   H9
        #
        # create_reversed_ethanol()
        #      H5   H2
        #      |    |
        # H4 - C8 - C7 - O6 - H0
        #      |    |
        #      H3   H1

        molecules = [
            Molecule.from_file(get_data_file_path("molecules/ethanol.sdf")),
            Molecule.from_file(get_data_file_path("molecules/ethanol_reordered.sdf")),
            create_reversed_ethanol(),
        ]
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [
            -0.2,
            -0.1,
            0.3,
            0.08,
            -0.02,
            -0.02,
            -0.02,
            -0.01,
            -0.01,
            -0.2,
            0.3,
            -0.1,
            0.08,
            -0.02,
            -0.02,
            -0.02,
            -0.01,
            -0.01,
            0.08,
            -0.01,
            -0.01,
            -0.02,
            -0.02,
            -0.02,
            0.3,
            -0.1,
            -0.2,
        ] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    @pytest.mark.parametrize("monatomic_ion,formal_charge", generate_monatomic_ions())
    def test_library_charges_monatomic_ions(self, monatomic_ion, formal_charge):
        """Test assigning library charges to each of the monatomic ions in openff-1.1.0.xml"""
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            get_data_file_path("test_forcefields/ion_charges.offxml"),
        )
        mol = Molecule.from_smiles("[{}]".format(monatomic_ion))
        omm_system = ff.create_openmm_system(mol.to_topology())

        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        q, _, _ = nonbondedForce.getParticleParameters(0)
        assert q == formal_charge

    def test_charge_method_hierarchy(self):
        """Ensure that molecules are parameterized by charge_from_molecules first, then library charges
        if not applicable, then AM1BCC otherwise"""
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_CH_zeroes_library_charges_xml,
            get_data_file_path("test_forcefields/tip3p.offxml"),
            xml_charge_increment_model_formal_charges,
        )
        # Cyclohexane will be assigned nonzero charges based on `charge_from_molecules` kwarg
        cyclohexane = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "cyclohexane.sdf"))
        )
        # Butanol will be assigned zero-valued charges based on `charge_from_molecules` kwarg
        butanol = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "butanol.sdf"))
        )
        # Propane will be assigned charges from AM1BCC
        propane = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "propane.sdf"))
        )
        # Water will be assigned TIP3P librarycharges
        water = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "water.sdf"))
        )
        # Ethanol should be caught by ToolkitAM1BCC
        ethanol = Molecule.from_file(
            get_data_file_path(os.path.join("systems", "monomers", "ethanol.sdf"))
        )
        # iodide should fail to be assigned charges by AM1-BCC, and instead be caught by ChargeIncrementHandler
        # (and assigned formal charge by dummy charge method)
        iodide = Molecule.from_smiles("[I-1]")

        # Assign dummy partial charges to cyclohexane, which we expect to find in the final system since it
        # is included in the charge_from_molecules kwarg to create_openmm_system
        cyclohexane.partial_charges = unit.Quantity(
            np.array(
                [
                    -0.2,
                    -0.2,
                    -0.2,
                    -0.2,
                    -0.2,
                    -0.2,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                    0.1,
                ]
            ),
            unit.elementary_charge,
        )

        # There were previously known issues when parameterizing molecules with all zero charges,
        # so test this explicitly with butanol. Since butanol will be in the charge_from_molecules kwarg,
        # we expect to find these charges in the final system.
        butanol.partial_charges = np.array([0.0] * 15) * unit.elementary_charge

        # Add dummy partial charges to propane, which should be IGNORED since it
        # isn't in the charge_from_molecules kwarg
        propane.partial_charges = np.array([99.0] * 11) * unit.elementary_charge

        # Add dummy partial charges to water, which should be IGNORED since it
        # isn't in the charge_from_molecules kwarg
        water.partial_charges = np.array([99.0] * 3) * unit.elementary_charge

        #            molecule       correct charge method
        molecules = [
            cyclohexane,  # charge_from_molecules kwarg
            butanol,  # charge_from_molecules kwarg
            propane,  # library charges
            water,  # library charges
            ethanol,  # AM1-BCC
            iodide,
        ]  # charge increment model (formal charge)
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(
            top, charge_from_molecules=[cyclohexane, butanol]
        )
        existing = [f for f in omm_system.getForces() if type(f) == NonbondedForce]

        # Ensure that the handlers do not make multiple NonbondedForce objects
        assert len(existing) == 1
        nonbondedForce = existing[0]
        expected_charges = [  # cyclohexane (18 atoms) should have the following values from charge_from_mols
            -0.2,
            -0.2,
            -0.2,
            -0.2,
            -0.2,
            -0.2,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            # butanol (15 atoms) should have the following values from charge_from_mols
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            # propane (11 atoms) should have the following values from xml_CH_zeroes_library_charges_xml
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            # water (3 atoms) should have the following charges from tip3p.offxml
            -0.834,
            0.417,
            0.417,
        ] * openmm_unit.elementary_charge

        # Ensure that the first four molecules have exactly the charges we intended
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

        # Ensure the last molecule (ethanol) had _some_ nonzero charge assigned by an AM1BCC implementation
        for particle_index in range(len(expected_charges), top.n_atoms - 1):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q != 0 * unit.elementary_charge

        # Ensure that iodine has a charge of -1, specified by charge increment model charge_method="formal charge"
        q, _, _ = nonbondedForce.getParticleParameters(top.n_atoms - 1)
        assert q == -1.0 * openmm_unit.elementary_charge

    def test_assign_charges_to_molecule_in_parts_using_multiple_library_charges(self):
        """Test assigning charges to parts of a molecule using two library charge lines. Note that these LibraryCharge
        SMIRKS have partial overlap, so this also tests that the hierarchy is correctly obeyed.
        """
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ethanol_library_charges_in_parts_ff,
        )

        molecules = [
            Molecule.from_file(get_data_file_path("molecules/ethanol.sdf")),
            Molecule.from_file(get_data_file_path("molecules/ethanol_reordered.sdf")),
        ]
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [
            -0.2,
            -0.1,
            0.3,
            0.08,
            -0.02,
            -0.02,
            -0.02,
            -0.01,
            -0.01,
            -0.2,
            0.3,
            -0.1,
            0.08,
            -0.02,
            -0.02,
            -0.02,
            -0.01,
            -0.01,
        ] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_assign_charges_using_library_charges_by_single_atoms(self):
        """Test assigning charges to parts of a molecule using per-atom library charges. Note that these LibraryCharge
        SMIRKS will match multiple atoms, so this is also a test of correct usage of the parameter hierarchy..
        """
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ethanol_library_charges_by_atom_ff,
        )

        molecules = [
            Molecule.from_file(get_data_file_path("molecules/ethanol.sdf")),
            Molecule.from_file(get_data_file_path("molecules/ethanol_reordered.sdf")),
        ]
        top = Topology.from_molecules(molecules)
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        expected_charges = [
            -0.2,
            -0.1,
            0.3,
            0.08,
            -0.02,
            -0.02,
            -0.02,
            -0.01,
            -0.01,
            -0.2,
            0.3,
            -0.1,
            0.08,
            -0.02,
            -0.02,
            -0.02,
            -0.01,
            -0.01,
        ] * openmm_unit.elementary_charge
        for particle_index, expected_charge in enumerate(expected_charges):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q == expected_charge

    def test_library_charges_dont_parameterize_molecule_because_of_incomplete_coverage(
        self,
    ):
        """Fail to assign charges to a molecule because not all atoms can be assigned"""
        molecules = [Molecule.from_file(get_data_file_path("molecules/toluene.sdf"))]
        top = Topology.from_molecules(molecules)

        # The library charges in the FF should not be able to fully cover toluene
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ethanol_library_charges_by_atom_ff,
        )
        # Delete the ToolkitAM1BCCHandler so the molecule won't get charges from anywhere
        del ff._parameter_handlers["ToolkitAM1BCC"]
        with pytest.raises(
            RuntimeError, match="Cc1ccccc1 could not be fully assigned charges"
        ):
            omm_system = ff.create_openmm_system(top)

        # If we do NOT delete the ToolkiAM1BCCHandler, then toluene should be assigned some nonzero partial charges.
        # The exact value will vary by toolkit, so we don't test that here.
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ethanol_library_charges_by_atom_ff,
        )
        omm_system = ff.create_openmm_system(top)
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        for particle_index in range(top.n_atoms):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert q != 0 * unit.elementary_charge

    @pytest.mark.parametrize(
        "charge_method,additional_offxmls",
        [
            ("ToolkitAM1BCC", []),
            ("LibraryCharges", [xml_ethanol_library_charges_ff]),
            ("ChargeIncrementHandler", [xml_charge_increment_model_formal_charges]),
            ("charge_from_molecules", []),
        ],
    )
    def test_charges_on_ref_mols_when_using_return_topology(
        self, charge_method, additional_offxmls
    ):
        """Ensure that charges are set on returned topology if the user specifies 'return_topology=True' in
        create_openmm_system"""
        # TODO: Should this test also cover multiple unique molecules?
        mol = create_acetate()
        ff = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            *additional_offxmls,
        )
        charge_mols = []
        if charge_method == "charge_from_molecules":
            mol.partial_charges = unit.Quantity(
                np.array([-1.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]), unit.elementary_charge
            )
            charge_mols = [mol]
        omm_system, ret_top = ff.create_openmm_system(
            mol.to_topology(), charge_from_molecules=charge_mols, return_topology=True
        )
        nonbondedForce = [
            f for f in omm_system.getForces() if type(f) == NonbondedForce
        ][0]
        ref_mol_from_ret_top = [i for i in ret_top.reference_molecules][0]

        # Make sure the charges on the molecule are all nonzero, and that the molecule's
        # partial_charges array matches the charges in the openmm system
        all_charges_zero = True
        for particle_index, ref_mol_charge in enumerate(
            ref_mol_from_ret_top.partial_charges
        ):
            q, _, _ = nonbondedForce.getParticleParameters(particle_index)
            assert from_openmm(q) == ref_mol_charge
            if q != 0.0 * unit.elementary_charge:
                all_charges_zero = False
        assert not (all_charges_zero)

    def test_library_charges_from_molecule_manual(self, force_field):
        """Test that constructing a LibraryChargeHandler from partial charges on Molecule objects
        produces the same result as using the `charge_from_molecules` kwarg. while manually
        setting the molecule's partial charges to arbitrary non-physical values"""
        # TODO: Remove this test if `charge_from_molecules` is depcreated (#806)
        mol = Molecule.from_mapped_smiles("[Cl:1][C:2]#[C:3][F:4]")
        mol.partial_charges = np.linspace(-0.3, 0.3, 4) * unit.elementary_charge

        using_kwarg = force_field.create_openmm_system(
            topology=mol.to_topology(), charge_from_molecules=[mol]
        )

        library_charges = LibraryChargeHandler.LibraryChargeType.from_molecule(mol)
        force_field.register_parameter_handler(LibraryChargeHandler(version=0.3))
        force_field["LibraryCharges"].add_parameter(parameter=library_charges)
        using_library_charges = force_field.create_openmm_system(
            topology=mol.to_topology()
        )

        compare_partial_charges(using_kwarg, using_library_charges)

    def test_library_charges_from_molecule_assigned(self, force_field):
        """Test that constructing a LibraryChargeHandler from partial charges on Molecule objects
        produces the same result as using the `charge_from_molecules` kwarg. while manually
        setting the molecule's partial charges to arbitrary non-physical values"""
        mol = Molecule.from_smiles("CCO")
        mol.assign_partial_charges(partial_charge_method="mmff94")

        using_kwarg = force_field.create_openmm_system(
            topology=mol.to_topology(), charge_from_molecules=[mol]
        )

        library_charges = LibraryChargeHandler.LibraryChargeType.from_molecule(mol)
        force_field.register_parameter_handler(LibraryChargeHandler(version=0.3))
        force_field["LibraryCharges"].add_parameter(parameter=library_charges)
        using_library_charges = force_field.create_openmm_system(
            topology=mol.to_topology()
        )

        compare_partial_charges(using_kwarg, using_library_charges)

    @requires_openeye
    def test_toolkit_am1bcc_uses_elf10_if_oe_is_available(self, force_field):
        """Ensure that the ToolkitAM1BCCHandler assigns ELF10 charges if OpenEye is available."""
        # Can't just use CCO - Molecule needs to be big enough to realistically
        # result in more than one conformer when ELF10 is requested
        mol = Molecule.from_smiles("OCCCCCCO")

        # The test forcefield "out of the box" just uses ToolkitAM1BCC, so this
        # should produce ELF10 charges
        sys = force_field.create_openmm_system(topology=mol.to_topology())

        # Make another system with a force field that has the ToolkitAM1BCCHandler
        # removed, and is explicitly told to use a ChargeIncrementModelHandler
        # with AM1BCCELF10 charges. Ensure that the resulting charges are identical.
        force_field.deregister_parameter_handler("ToolkitAM1BCC")
        cimh = force_field.get_parameter_handler(
            "ChargeIncrementModel",
            {"version": 0.3, "partial_charge_method": "am1bccelf10"},
        )

        sys_a1b_elf10 = force_field.create_openmm_system(topology=mol.to_topology())

        compare_partial_charges(sys, sys_a1b_elf10)

        # Now modify the ChargeIncrementModelHandler to specifically use single-conformer
        # AM1BCC and ensure that the resulting charges are different from the first system.
        cimh.partial_charge_method = "am1bcc"
        assert force_field["ChargeIncrementModel"].partial_charge_method == "am1bcc"

        sys_a1b_single_conf = force_field.create_openmm_system(
            topology=mol.to_topology()
        )
        # Ensure that we get different results with AM1BCC ELF10 and single conf AM1BCC
        with pytest.raises(AssertionError):
            compare_partial_charges(sys_a1b_elf10, sys_a1b_single_conf)


class TestForceFieldConstraints:
    """Tests that constraints are correctly applied and behave correctly."""

    @classmethod
    def check_molecule_constraints(cls, molecule, system, bond_elements, bond_length):
        """Check that the bonds in the molecule is correctly constrained."""
        for constraint_idx in range(system.getNumConstraints()):
            atom1_idx, atom2_idx, distance = system.getConstraintParameters(
                constraint_idx
            )
            atom_elements = {
                molecule.atoms[atom1_idx].symbol,
                molecule.atoms[atom2_idx].symbol,
            }
            assert atom_elements == bond_elements
            assert np.isclose(
                distance.value_in_unit(openmm_unit.angstrom),
                bond_length.m_as(unit.angstrom),
            )

    def test_constraints_hbonds(self):
        """Test that hydrogen bonds constraints are applied correctly to a ethane molecule."""
        # Parametrize an ethane molecule.
        ethane = Molecule.from_smiles("CC")
        topology = Topology.from_molecules([ethane])
        ff = ForceField(
            XML_FF_GENERICS, get_data_file_path("test_forcefields/old/hbonds.offxml")
        )
        system = ff.create_openmm_system(topology)

        # Check that all C-H bonds have been constrained to the FF bond length.
        self.check_molecule_constraints(
            ethane, system, bond_elements={"C", "H"}, bond_length=1.09 * unit.angstrom
        )


def generate_alkethoh_parameters_assignment_cases():
    """Create dynamically all test cases that should be ran for the AlkEthOH set."""
    # These AlkEthOH molecules are always run by test_alkethoh_parameters_assignment.
    fast_test_cases = ["r0", "r12", "r118", "c38", "c100", "c1161", "c1266"]

    def extract_id(file_path):
        """Extract the AlkEthOH molecule ID from the file path."""
        # An example of file path is AlkEthOH_tripos/AlkEthOH_chain_filt1/AlkEthOH_c555.crd
        return os.path.splitext(os.path.basename(file_path))[0][9:]

    # Get all the molecules ids from the tarfiles. The tarball is extracted
    # in conftest.py if slow tests are activated.
    import tarfile

    alkethoh_tar_file_path = get_data_file_path(
        os.path.join("molecules", "AlkEthOH_tripos.tar.gz")
    )
    with tarfile.open(alkethoh_tar_file_path, "r:gz") as tar:
        # Collect all the files discarding the duplicates in the test_filt1 folder.
        slow_test_cases = {
            extract_id(m.name)
            for m in tar.getmembers()
            if "crd" in m.name and "test_filt1" not in m.name
        }

    # Remove fast test cases from slow ones to avoid duplicate tests.
    # Remove also water (c1302), which was reparameterized in AlkEthOH
    # to be TIP3P (not covered by Frosst_AlkEthOH_parmAtFrosst.
    for fast_test_case in fast_test_cases + ["c1302"]:
        slow_test_cases.remove(fast_test_case)

    # Mark all slow cases as slow.
    slow_test_cases = [
        pytest.param(case, marks=pytest.mark.slow) for case in sorted(slow_test_cases)
    ]

    # Isolate the AlkEthOH ID.
    return fast_test_cases + slow_test_cases


def generate_freesolv_parameters_assignment_cases():
    """Create dynamically all test cases that should be ran for the FreeSolv set."""
    import tarfile

    # For these tests, UndefinedStereochemistryError is ignored.
    # The chirality was manually checked (see issue #175).
    ignore_undefined_stereo = {"2501588", "3266352", "7829570"}

    # These molecules are always tested by test_freesolv_parameters_assignment().
    # Each test case is (freesolv_id, force_field_version, allow_undefined_stereo).
    fast_test_cases = [
        ("1019269", "0_0_4_fixed", False),
        (
            "63712",
            "0_0_2",
            False,
        ),  # The XML was regenerated after fixing the issue described in #179.
        ("1723043", "0_0_2", False),
        ("2501588", "0_0_2", True),  # Test impropers and undefined stereochemistry.
        (
            "3323117",
            "0_0_2",
            False,
        ),  # The XML was regenerated after fixing the issue described in #179.
        ("1107178", "0_0_2", False),  # Molecule with iodine
        ("1036761", "0_0_2", False),  # Molecule with primary amine
    ]

    def extract_id(file_path):
        """Extract the FreeSolv ID and force field version from the file subpath."""
        # An example of file path is FreeSolv/xml_0_0_4_fixed/mobley_7913234_vacuum.xml
        freesolv_id = os.path.basename(file_path).split("_")[1]
        force_field_version = os.path.basename(os.path.dirname(file_path))[4:]
        allow_undefined_stereo = freesolv_id in ignore_undefined_stereo
        return (freesolv_id, force_field_version, allow_undefined_stereo)

    # Get all the tarball XML files available. The tarball is extracted
    # in conftest.py if slow tests are activated.
    freesolv_tar_file_path = get_data_file_path(
        os.path.join("molecules", "FreeSolv.tar.gz")
    )
    with tarfile.open(freesolv_tar_file_path, "r:gz") as tar:
        slow_test_cases = {
            extract_id(m.name) for m in tar.getmembers() if ".xml" in m.name
        }

    # Remove fast test cases from slow ones to avoid duplicate tests.
    for fast_test_case in fast_test_cases:
        slow_test_cases.remove(fast_test_case)

    # Mark all slow cases as slow.
    slow_test_cases = [
        pytest.param(*case, marks=pytest.mark.slow) for case in sorted(slow_test_cases)
    ]

    return fast_test_cases + slow_test_cases


class TestForceFieldParameterAssignment(_ForceFieldFixtures):
    """Regression tests checking that parameters are assigned correctly."""

    @requires_openeye_mol2
    @pytest.mark.parametrize(
        "alkethoh_id", generate_alkethoh_parameters_assignment_cases()
    )
    def test_alkethoh_parameters_assignment(self, alkethoh_id, alkethoh_forcefield):
        """Test that ForceField assign parameters correctly in the AlkEthOH set.
        The test compares the System parameters of a AlkEthOH molecule
        parameterized with AMBER and Frosst_AlkEthOH_parmAtFrosst.offxml.

        The AMBER files were prepared following the pipeline described here:
            https://github.com/openforcefield/open-forcefield-data/tree/master/Model-Systems/AlkEthOH_distrib/
        They were built for the SMIRNOFF parametrization to yield exact same
        parameters.

        The AlkEthOH set, however, does not have impropers, which should be
        tested separately. Currently, test_freesolv_parameters_assignment
        does the job.

        """
        from openff.toolkit.tests.utils import (
            compare_amber_smirnoff,
            get_alkethoh_file_path,
        )

        # Obtain the path to the input files.
        alkethoh_name = "AlkEthOH_" + alkethoh_id
        mol2_filepath, top_filepath, crd_filepath = get_alkethoh_file_path(
            alkethoh_name, get_amber=True
        )

        # Load molecule.
        molecule = Molecule.from_file(mol2_filepath)

        # Compare parameters. Skip the energy checks as the parameter check should be
        # sufficient. We test both energies and parameters in the slow test.
        # We ignore the charges for now as they are not included in the force field.
        # TODO: Reactivate the charge check when we'll be able to load charges from files.
        compare_amber_smirnoff(
            top_filepath,
            crd_filepath,
            alkethoh_forcefield,
            molecule,
            check_energies=False,
            ignore_charges=True,
        )

    @requires_openeye_mol2
    def test_multi_alkethoh_parameters_assignment(self, alkethoh_forcefield):
        """Test that systems with multiple reference molecules are parametrized correctly.

        The test relies on the fact that we have already verified we can
        parametrize correctly single AlkEthOH molecules in
        test_alkethoh_parameters_assignment(). We use ParmEd to merge
        the AMBER files to be used as reference parameters.

        """
        import parmed

        from openff.toolkit.tests.utils import (
            compare_system_energies,
            compare_system_parameters,
            get_alkethoh_file_path,
        )

        # The AlkEthOH molecule ids to mix in the systems.
        alketoh_ids = ["r0", "c38", "c1161"]

        # Load molecules and structures.
        molecules = []
        structures = []
        for alkethoh_id in alketoh_ids:
            mol2_filepath, top_filepath, crd_filepath = get_alkethoh_file_path(
                "AlkEthOH_" + alkethoh_id, get_amber=True
            )
            molecules.append(Molecule.from_file(mol2_filepath))
            amber_parm = parmed.load_file(top_filepath, crd_filepath)
            # Convert this into a real structure as mixing AmberParm objects is bugged (see ParmEd#1045).
            structures.append(amber_parm.copy(parmed.Structure))

        # Merge the structures into a single system with two copies of the last molecule.
        structure_mixture = (
            structures[0] + structures[1] + structures[2] + structures[-1]
        )
        amber_system = structure_mixture.createSystem(
            nonbondedMethod=openmm.app.NoCutoff
        )

        # Create the OpenFF System through ForceField.
        topology = Topology.from_openmm(
            structure_mixture.topology, unique_molecules=molecules
        )
        topology.box_vectors = None
        off_system = alkethoh_forcefield.create_openmm_system(topology)

        # Translate the molecules a little to avoid overlapping atoms.
        positions = copy.deepcopy(structure_mixture.positions)
        translate_vectors = [
            np.array([1.0, 0.0, 0.0]) * openmm_unit.nanometer,
            np.array([0.0, 1.0, 0.0]) * openmm_unit.nanometer,
            np.array([0.0, 0.0, 1.0]) * openmm_unit.nanometer,
            # Leave the fourth molecule where it is.
        ]
        current_atom_idx = 0
        for translate_vector, mol in zip(translate_vectors, molecules):
            n_mol_atoms = len(mol.atoms)
            positions[
                current_atom_idx : current_atom_idx + n_mol_atoms
            ] += translate_vector
            current_atom_idx += n_mol_atoms

        # Compare parameters and systems.
        # TODO: Reactivate charges comparison when we'll be able to read them from the file.
        compare_system_parameters(
            amber_system,
            off_system,
            systems_labels=("AMBER", "SMIRNOFF"),
            ignore_charges=True,
        )
        compare_system_energies(
            amber_system, off_system, positions, ignore_charges=True
        )

    @requires_openeye_mol2
    @pytest.mark.parametrize(
        ("freesolv_id", "forcefield_version", "allow_undefined_stereo"),
        generate_freesolv_parameters_assignment_cases(),
    )
    def test_freesolv_parameters_assignment(
        self, freesolv_id, forcefield_version, allow_undefined_stereo
    ):
        """Regression test on parameters assignment based on the FreeSolv set used in the 0.1 paper.

        This, contrarily to the similar AlkEthOH test, checks also constraints
        and improper torsions.

        """
        from openff.toolkit.tests.utils import (
            compare_system_parameters,
            get_freesolv_file_path,
        )

        mol2_file_path, xml_file_path = get_freesolv_file_path(
            freesolv_id, forcefield_version
        )

        # Load molecules.
        molecule = Molecule.from_file(
            mol2_file_path, allow_undefined_stereo=allow_undefined_stereo
        )

        # Create OpenFF System with the current toolkit.
        forcefield_file_path = get_data_file_path(
            "test_forcefields/old/test_ff_" + forcefield_version + "_spec_0_2.offxml"
        )
        ff = ForceField(
            forcefield_file_path,
            get_data_file_path("test_forcefields/old/hbonds.offxml"),
        )
        ff_system = ff.create_openmm_system(molecule.to_topology())

        # Load OpenMM System created with the 0.1 version of the toolkit.
        with open(xml_file_path, "r") as f:
            xml_system = openmm.XmlSerializer.deserialize(f.read())

        # Compare parameters. We ignore the improper folds as in 0.0.3 we
        # used a six-fold implementation while we now use a three-fold way.
        # TODO: Reactivate charge comparison once we'll be able to read them from file.
        compare_system_parameters(
            ff_system,
            xml_system,
            systems_labels=("current OpenFF", "SMIRNOFF 0.0.4"),
            ignore_charges=True,
            ignore_improper_folds=True,
        )

    @unimplemented_interchange
    @requires_openeye_mol2
    @pytest.mark.parametrize(("is_periodic"), (False, True))
    @pytest.mark.parametrize(("gbsa_model"), ["HCT", "OBC1", "OBC2"])
    @pytest.mark.parametrize(
        ("freesolv_id", "forcefield_version", "allow_undefined_stereo"),
        generate_freesolv_parameters_assignment_cases(),
    )
    def test_freesolv_gbsa_energies(
        self,
        gbsa_model,
        is_periodic,
        freesolv_id,
        forcefield_version,
        allow_undefined_stereo,
    ):
        """
        Regression test on HCT, OBC1, and OBC2 GBSA models. This test ensures that the
        SMIRNOFF-based GBSA models match the equivalent OpenMM implementations.
        """
        import parmed as pmd

        from openff.toolkit.tests.utils import (
            compare_system_energies,
            create_system_from_amber,
            get_context_potential_energy,
            get_freesolv_file_path,
        )

        mol2_file_path, _ = get_freesolv_file_path(freesolv_id, forcefield_version)

        # Load molecules.
        molecule = Molecule.from_file(
            mol2_file_path, allow_undefined_stereo=allow_undefined_stereo
        )

        # Give each atom a unique name, otherwise OpenMM will complain
        for idx, atom in enumerate(molecule.atoms):
            atom.name = f"{atom.symbol}{idx}"
        positions = to_openmm(molecule.conformers[0])

        off_gbsas = {
            "HCT": "test_forcefields/GBSA_HCT-1.0.offxml",
            "OBC1": "test_forcefields/GBSA_OBC1-1.0.offxml",
            "OBC2": "test_forcefields/GBSA_OBC2-1.0.offxml",
        }
        # Create OpenFF System with the current toolkit.
        ff = ForceField(
            "test_forcefields/test_forcefield.offxml", off_gbsas[gbsa_model]
        )

        # OpenMM 7.7 and older don't properly handle parsing prmtop files with a GBSA model and
        # a switching function; this bug may be fixed but for testing purposes, simply turn off
        # the switching function
        ff["vdW"].switch_width = 0.0 * unit.angstrom

        off_top = molecule.to_topology()
        if is_periodic:
            off_top.box_vectors = (
                (30.0, 0, 0),
                (0, 30.0, 0),
                (0, 0, 30.0),
            ) * unit.angstrom

        else:
            off_top.box_vectors = None

        off_omm_system = ff.create_openmm_system(
            off_top, charge_from_molecules=[molecule]
        )

        off_nonbonded_force = [
            force
            for force in off_omm_system.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]

        omm_top = off_top.to_openmm()
        pmd_struct = pmd.openmm.load_topology(omm_top, off_omm_system, positions)
        prmtop_file = NamedTemporaryFile(suffix=".prmtop")
        inpcrd_file = NamedTemporaryFile(suffix=".inpcrd")
        pmd_struct.save(prmtop_file.name, overwrite=True)
        pmd_struct.save(inpcrd_file.name, overwrite=True)

        openmm_gbsas = {
            "HCT": openmm.app.HCT,
            "OBC1": openmm.app.OBC1,
            "OBC2": openmm.app.OBC2,
        }

        # The functional form of the nonbonded force will change depending on whether the cutoff
        # is None during initialization. Therefore, we need to figure that out here.

        # WARNING: The NonbondedMethod enums at openmm.app.forcefield and openmm.CustomGBForce
        # aren't necessarily the same, and could be misinterpreted if the wrong one is used. For
        # create_system_from_amber, we must provide the app.forcefield version.

        if is_periodic:
            amber_nb_method = openmm.app.forcefield.CutoffPeriodic
            amber_cutoff = off_nonbonded_force.getCutoffDistance()
        else:
            amber_nb_method = openmm.app.forcefield.NoCutoff
            amber_cutoff = None

        # OpenMM will process False to mean "don't use one" and a quantity as "use with this distance"
        switching_distance = off_nonbonded_force.getUseSwitchingFunction()
        if switching_distance:
            switching_distance = off_nonbonded_force.getSwitchingDistance()

        (
            amber_omm_system,
            amber_omm_topology,
            amber_positions,
        ) = create_system_from_amber(
            prmtop_file.name,
            inpcrd_file.name,
            implicitSolvent=openmm_gbsas[gbsa_model],
            nonbondedMethod=amber_nb_method,
            nonbondedCutoff=amber_cutoff,
            gbsaModel="ACE",
            implicitSolventKappa=0.0,
            switchDistance=switching_distance,
        )

        # Retrieve the GBSAForce from both the AMBER and OpenForceField systems
        off_gbsa_forces = list()
        for force in off_omm_system.getForces():
            if isinstance(force, openmm.GBSAOBCForce) or isinstance(
                force, openmm.openmm.CustomGBForce
            ):
                off_gbsa_forces.append(force)

        amber_gbsa_forces = list()
        for force in amber_omm_system.getForces():
            if isinstance(force, openmm.GBSAOBCForce) or isinstance(
                force, openmm.openmm.CustomGBForce
            ):
                amber_gbsa_forces.append(force)

        assert len(amber_gbsa_forces) == len(off_gbsa_forces) == 1

        off_gbsa_force = off_gbsa_forces[0]
        amber_gbsa_force = amber_gbsa_forces[0]

        # We get radius and screen values from each model's getStandardParameters method
        if gbsa_model == "HCT":
            gb_params = (
                openmm.app.internal.customgbforces.GBSAHCTForce.getStandardParameters(
                    omm_top
                )
            )
        elif gbsa_model == "OBC1":
            gb_params = (
                openmm.app.internal.customgbforces.GBSAOBC1Force.getStandardParameters(
                    omm_top
                )
            )
        elif gbsa_model == "OBC2":
            gb_params = (
                openmm.app.internal.customgbforces.GBSAOBC2Force.getStandardParameters(
                    omm_top
                )
            )

        # Use GB params from OpenMM GBSA classes to populate parameters
        for idx, (radius, screen) in enumerate(gb_params):
            # Keep the charge, but throw out the old radius and screen values
            q, _, _ = amber_gbsa_force.getParticleParameters(idx)
            if isinstance(amber_gbsa_force, openmm.GBSAOBCForce):
                # Note that in GBSAOBCForce, the per-particle parameters are separate
                # arguments, while in CustomGBForce they're a single iterable
                amber_gbsa_force.setParticleParameters(idx, q, radius, screen)

            elif isinstance(amber_gbsa_force, openmm.CustomGBForce):
                # !!! WARNING: CustomAmberGBForceBase expects different per-particle parameters
                # depending on whether you use addParticle or setParticleParameters. In
                # setParticleParameters, we have to apply the offset and scale BEFORE setting
                # parameters, whereas in addParticle, it is applied afterwards, and the particle
                # parameters are not set until an auxillary finalize() method is called. !!!
                amber_gbsa_force.setParticleParameters(
                    idx, (q, radius - 0.009, screen * (radius - 0.009))
                )

        # Put the GBSA force into a separate group so we can specifically compare GBSA energies
        amber_gbsa_force.setForceGroup(1)
        off_gbsa_force.setForceGroup(1)

        # Some manual overrides to get the OFF system's NonbondedForce matched up with the AMBER system
        if is_periodic:
            off_nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            off_nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

        off_nonbonded_force.setReactionFieldDielectric(1.0)

        # Create Contexts
        integrator = openmm.VerletIntegrator(1.0 * openmm_unit.femtoseconds)
        platform = Platform.getPlatformByName("Reference")
        amber_context = openmm.Context(amber_omm_system, integrator, platform)
        off_context = openmm.Context(
            off_omm_system, copy.deepcopy(integrator), platform
        )

        # Get context energies
        amber_energy = get_context_potential_energy(amber_context, positions)
        off_energy = get_context_potential_energy(off_context, positions)
        # Very handy for debugging
        # print(openmm.XmlSerializer.serialize(off_gbsa_force))
        # print(openmm.XmlSerializer.serialize(amber_gbsa_force))
        # Ensure that the GBSA energies (which we put into ForceGroup 1) are identical
        # For Platform=OpenCL, we do get "=="-level identical numbers, but for "Reference", we don't.
        # assert amber_energy[1] == off_energy[1]
        assert (
            abs(amber_energy[1] - off_energy[1]) < 1e-5 * openmm_unit.kilojoule_per_mole
        )

        # Ensure that all system energies are the same
        compare_system_energies(
            off_omm_system,
            amber_omm_system,
            positions,
            by_force_type=False,
            atol=1.0e-4,
        )

    def test_tip5p_dimer_energy(self):
        pytest.skip("fails until migration to openff-interchange")

        from openff.toolkit.tests.utils import evaluate_molecules_off

        off_ff = ForceField("test_forcefields/test_forcefield.offxml", xml_tip5p)

        molecule1 = create_water()
        molecule1.atoms[1].name = "O"
        molecule1.atoms[0].name = "H1"
        molecule1.atoms[2].name = "H2"
        molecule1.generate_conformers(n_conformers=1)

        molecule1.conformers[0] = unit.Quantity(
            np.array(
                [
                    [-0.78900161, -0.19816432, -0.0],
                    [-0.00612716, 0.39173634, -0.0],
                    [0.79512877, -0.19357202, 0.0],
                ]
            ),
            unit.angstrom,
        )

        molecule2 = create_water()
        molecule2.atoms[1].name = "O"
        molecule2.atoms[0].name = "H1"
        molecule2.atoms[2].name = "H2"
        molecule2.generate_conformers(n_conformers=1)

        # This is mol1 + 10 angstrom
        molecule2.conformers[0] = unit.Quantity(
            np.array(
                [
                    [9.21099839, 9.80183568, 10.0],
                    [9.99387284, 10.39173634, 10.0],
                    [10.79512877, 9.80642798, 10.0],
                ]
            ),
            unit.angstrom,
        )

        off_crds, off_ene = evaluate_molecules_off(
            [molecule1, molecule2], off_ff, minimize=False
        )

        ref_crds_with_vsite = unit.Quantity(
            np.array(
                [
                    [-0.0789001605855, -0.0198164316973, -0.0],
                    [-0.0006127160367, 0.0391736336129, -0.0],
                    [0.0795128766222, -0.0193572019156, 0.0],
                    [0.9210998394144, 0.9801835683026, 1.0],
                    [0.9993872839632, 1.0391736336129, 1.0],
                    [1.0795128766222, 0.9806427980843, 1.0],
                    [-0.0012451030403, 0.0796049187395, -0.0571394020767],
                    [-0.0012451030403, 0.0796049187395, 0.0571394020767],
                    [0.9987548969596, 1.0796049187395, 0.9428605979232],
                    [0.9987548969596, 1.0796049187395, 1.0571394020767],
                ]
            ),
            unit.nanometer,
        )

        ref_ene = 0.0011797690240 * openmm_unit.kilojoule_per_mole

        assert np.allclose(
            off_crds.m_as(unit.angstrom),
            ref_crds_with_vsite.m_as(unit.angstrom),
        )
        # allow 1% error in energy difference (default is .001%)
        assert np.allclose(
            off_ene.value_in_unit(openmm_unit.kilocalorie_per_mole),
            ref_ene.value_in_unit(openmm_unit.kilocalorie_per_mole),
            rtol=0.05,
        )

        # skip direct comparison for now because of particle creation and atom
        # ordering complications

        # # Create the oMM paramterized system
        # omm_ff = app.ForceField("tip5p.xml")

        # omm_top = off_top.to_openmm()
        # omm_modeller = app.Modeller(omm_top, positions)
        # omm_modeller.addExtraParticles(omm_ff)
        # omm_top = omm_modeller.getTopology()

        # omm_omm_system = omm_ff.createSystem(omm_top, nonbondedMethod=app.NoCutoff)

        # compare_system_energies(
        #     off_omm_system,
        #     omm_omm_system,
        #     positions,
        #     by_force_type=False,
        #     modify_system=False,
        # )

    @unimplemented_interchange
    @requires_openeye_mol2
    @pytest.mark.parametrize("zero_charges", [True, False])
    @pytest.mark.parametrize(("gbsa_model"), ["HCT", "OBC1", "OBC2"])
    def test_molecule_energy_gb_no_sa(self, zero_charges, gbsa_model):
        """Test creating a GBSA system without a surface energy term, and validate its energy
        against the same system made using OpenMM's AMBER GBSA functionality"""
        import parmed as pmd

        from openff.toolkit.tests.utils import (
            compare_system_energies,
            create_system_from_amber,
            get_context_potential_energy,
        )

        # Load an arbitrary molecule from the freesolv set
        molecule = Molecule.from_file(
            get_data_file_path("molecules/FreeSolv/mol2files_sybyl/mobley_1036761.mol2")
        )

        molecule.name = (
            "mobley_1036761"  # Name the molecule, otherwise OpenMM will complain
        )
        if zero_charges:
            molecule.partial_charges = (
                np.zeros(molecule.n_atoms) * unit.elementary_charge
            )

        # Give each atom a unique name, otherwise OpenMM will complain
        for idx, atom in enumerate(molecule.atoms):
            atom.name = f"{atom.symbol}{idx}"

        positions = np.concatenate(
            (molecule.conformers[0], molecule.conformers[0] + (10 * unit.angstrom))
        )
        # Create OpenFF System with the current toolkit.

        off_gbsas = {
            "HCT": "test_forcefields/GBSA_HCT-1.0.offxml",
            "OBC1": "test_forcefields/GBSA_OBC1-1.0.offxml",
            "OBC2": "test_forcefields/GBSA_OBC2-1.0.offxml",
        }

        ff = ForceField(
            "test_forcefields/test_forcefield.offxml", off_gbsas[gbsa_model]
        )
        ff.get_parameter_handler("GBSA").sa_model = None
        off_top = Topology.from_molecules([molecule, molecule])
        off_omm_system = ff.create_openmm_system(
            off_top, charge_from_molecules=[molecule]
        )

        omm_top = off_top.to_openmm()
        pmd_struct = pmd.openmm.load_topology(omm_top, off_omm_system, positions)
        prmtop_file = NamedTemporaryFile(suffix=".prmtop")
        inpcrd_file = NamedTemporaryFile(suffix=".inpcrd")
        pmd_struct.save(prmtop_file.name, overwrite=True)
        pmd_struct.save(inpcrd_file.name, overwrite=True)

        openmm_gbsas = {
            "HCT": openmm.app.HCT,
            "OBC1": openmm.app.OBC1,
            "OBC2": openmm.app.OBC2,
        }

        (
            amber_omm_system,
            amber_omm_topology,
            amber_positions,
        ) = create_system_from_amber(
            prmtop_file.name,
            inpcrd_file.name,
            implicitSolvent=openmm_gbsas[gbsa_model],
            nonbondedMethod=openmm.app.forcefield.NoCutoff,
            nonbondedCutoff=None,
            gbsaModel=None,
            implicitSolventKappa=0.0,
        )

        # Retrieve the GBSAForce from both the AMBER and OpenForceField systems
        off_gbsa_forces = list()
        for force in off_omm_system.getForces():
            if isinstance(force, openmm.GBSAOBCForce) or isinstance(
                force, openmm.openmm.CustomGBForce
            ):
                off_gbsa_forces.append(force)

        amber_gbsa_forces = list()
        for force in amber_omm_system.getForces():
            if isinstance(force, openmm.GBSAOBCForce) or isinstance(
                force, openmm.openmm.CustomGBForce
            ):
                amber_gbsa_forces.append(force)

        assert len(amber_gbsa_forces) == len(off_gbsa_forces) == 1

        off_gbsa_force = off_gbsa_forces[0]
        amber_gbsa_force = amber_gbsa_forces[0]

        # We get radius and screen values from each model's getStandardParameters method
        if gbsa_model == "HCT":
            gb_params = (
                openmm.app.internal.customgbforces.GBSAHCTForce.getStandardParameters(
                    omm_top
                )
            )
        elif gbsa_model == "OBC1":
            gb_params = (
                openmm.app.internal.customgbforces.GBSAOBC1Force.getStandardParameters(
                    omm_top
                )
            )
        elif gbsa_model == "OBC2":
            gb_params = (
                openmm.app.internal.customgbforces.GBSAOBC2Force.getStandardParameters(
                    omm_top
                )
            )
            # This is only necessary until https://github.com/openmm/openmm/pull/2362 is bundled into a conda release
            amber_gbsa_force.setSurfaceAreaEnergy(0)

        # Use GB params from OpenMM GBSA classes to populate parameters
        for idx, (radius, screen) in enumerate(gb_params):
            # Keep the charge, but throw out the old radius and screen values
            q, _, _ = amber_gbsa_force.getParticleParameters(idx)

            if isinstance(amber_gbsa_force, openmm.GBSAOBCForce):
                # Note that in GBSAOBCForce, the per-particle parameters are separate
                # arguments, while in CustomGBForce they're a single iterable
                amber_gbsa_force.setParticleParameters(idx, q, radius, screen)

            elif isinstance(amber_gbsa_force, openmm.CustomGBForce):
                # !!! WARNING: CustomAmberGBForceBase expects different per-particle parameters
                # depending on whether you use addParticle or setParticleParameters. In
                # setParticleParameters, we have to apply the offset and scale BEFORE setting
                # parameters, whereas in addParticle, it is applied afterwards, and the particle
                # parameters are not set until an auxillary finalize() method is called. !!!
                amber_gbsa_force.setParticleParameters(
                    idx, (q, radius - 0.009, screen * (radius - 0.009))
                )

        # Put the GBSA force into a separate group so we can specifically compare GBSA energies
        amber_gbsa_force.setForceGroup(1)
        off_gbsa_force.setForceGroup(1)

        # Create Contexts
        integrator = openmm.VerletIntegrator(1.0 * openmm_unit.femtoseconds)
        platform = Platform.getPlatformByName("Reference")
        amber_context = openmm.Context(amber_omm_system, integrator, platform)
        off_context = openmm.Context(
            off_omm_system, copy.deepcopy(integrator), platform
        )

        # Get context energies
        amber_energy = get_context_potential_energy(amber_context, to_openmm(positions))
        off_energy = get_context_potential_energy(off_context, to_openmm(positions))

        # Very handy for debugging
        # print(openmm.XmlSerializer.serialize(off_gbsa_force))
        # print(openmm.XmlSerializer.serialize(amber_gbsa_force))

        # Ensure that the GBSA energies (which we put into ForceGroup 1) are identical
        # For Platform=OpenCL, we do get "=="-level identical numbers, but for "Reference", we don't.
        # assert amber_energy[1] == off_energy[1]
        assert (
            abs(amber_energy[1] - off_energy[1]) < 1e-5 * openmm_unit.kilojoule_per_mole
        )

        # If charges are zero, the GB energy component should be 0, so the total GBSA energy should be 0
        if zero_charges:
            assert amber_energy[1]._value == 0.0
        else:
            assert amber_energy[1]._value != 0.0

        # Ensure that all system energies are the same
        compare_system_energies(
            off_omm_system, amber_omm_system, to_openmm(positions), by_force_type=False
        )

    @pytest.mark.slow
    @requires_openeye_mol2
    @pytest.mark.parametrize("toolkit_registry", toolkit_registries)
    def test_parameterize_protein(
        self,
        toolkit_registry,
        force_field,
    ):
        """Test that ForceField assigns parameters correctly for a protein"""

        mol_path = get_data_file_path("proteins/T4-protein.mol2")
        molecule = Molecule.from_file(mol_path, allow_undefined_stereo=False)
        topology = Topology.from_molecules(molecule)

        labels = force_field.label_molecules(topology)[0]

        assert len(labels["Bonds"]) == 2654
        assert len(labels["Angles"]) == 4789
        assert len(labels["ProperTorsions"]) == 6973
        assert len(labels["ImproperTorsions"]) == 528

        force_field.create_openmm_system(
            topology,
            charge_from_molecules=[molecule],
            toolkit_registry=toolkit_registry,
        )

    def test_modified_14_factors(self, force_field):
        """Test that the 1-4 scaling factors for electrostatics and vdW handlers matche,
        to a tight precision, the values specified in the force field."""
        top = Molecule.from_smiles("CCCC").to_topology()
        default_14 = copy.deepcopy(force_field)
        e_mod_14 = copy.deepcopy(force_field)
        vdw_mod_14 = copy.deepcopy(force_field)

        e_mod_14["Electrostatics"].scale14 = 0.66
        assert e_mod_14["Electrostatics"].scale14 == 0.66

        vdw_mod_14["vdW"].scale14 = 0.777
        assert vdw_mod_14["vdW"].scale14 == 0.777

        default_omm_sys = default_14.create_openmm_system(top)
        e_mod_omm_sys = e_mod_14.create_openmm_system(top)
        vdw_mod_omm_sys = vdw_mod_14.create_openmm_system(top)

        for omm_sys, expected_vdw_14, expected_coul_14 in [
            [default_omm_sys, 0.5, 0.833333],
            [e_mod_omm_sys, 0.5, 0.66],
            [vdw_mod_omm_sys, 0.777, 0.833333],
        ]:
            found_coul_14, found_vdw_14 = get_14_scaling_factors(omm_sys)

            np.testing.assert_almost_equal(
                actual=found_vdw_14,
                desired=expected_vdw_14,
                decimal=10,
                err_msg="vdW 1-4 scaling factors do not match",
            )

            np.testing.assert_almost_equal(
                actual=found_coul_14,
                desired=expected_coul_14,
                decimal=10,
                err_msg="Electrostatics 1-4 scaling factors do not match",
            )

    def test_14_missing_nonbonded_handler(self, force_field):
        """Test that something sane happens with 1-4 scaling factors if a
        ForceField is missing a vdWHandler and/or ElectrostaticsHandler"""
        top = Molecule.from_smiles("CCCC").to_topology()

        ff_no_vdw = copy.deepcopy(force_field)
        ff_no_electrostatics = copy.deepcopy(force_field)

        ff_no_vdw.deregister_parameter_handler("vdW")

        ff_no_electrostatics.deregister_parameter_handler("Electrostatics")
        ff_no_electrostatics.deregister_parameter_handler("ToolkitAM1BCC")

        sys_no_electrostatics = ff_no_electrostatics.create_openmm_system(top)

        np.testing.assert_almost_equal(
            actual=get_14_scaling_factors(sys_no_electrostatics)[1],
            desired=ff_no_electrostatics["vdW"].scale14,
            decimal=8,
        )

        sys_no_vdw = ff_no_vdw.create_openmm_system(top)
        np.testing.assert_almost_equal(
            actual=get_14_scaling_factors(sys_no_vdw)[0],
            desired=ff_no_vdw["Electrostatics"].scale14,
            decimal=8,
        )

    def test_overwrite_bond_orders(self):
        """Test that previously-defined bond orders in the topology are overwritten"""
        mol = create_ethanol()
        mol.assign_fractional_bond_orders(bond_order_model="am1-wiberg")
        top = Topology.from_molecules(mol)

        mod_mol = create_ethanol()
        mod_mol.assign_fractional_bond_orders(bond_order_model="am1-wiberg")
        # populate the mol with garbage bond orders
        for bond in mod_mol.bonds:
            bond.fractional_bond_order = 3 - bond.fractional_bond_order
        mod_top = Topology.from_molecules(mod_mol)

        default_bo = [
            b.fractional_bond_order for b in [*top.reference_molecules][0].bonds
        ]
        mod_bo = [
            b.fractional_bond_order for b in [*mod_top.reference_molecules][0].bonds
        ]
        assert not default_bo == mod_bo

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"), xml_ff_bo
        )

        omm_system, omm_sys_top = forcefield.create_openmm_system(
            top, return_topology=True
        )
        mod_omm_system, mod_omm_sys_top = forcefield.create_openmm_system(
            mod_top, return_topology=True
        )

        default_bond_force = [
            f for f in omm_system.getForces() if isinstance(f, openmm.HarmonicBondForce)
        ][0]
        mod_bond_force = [
            f
            for f in mod_omm_system.getForces()
            if isinstance(f, openmm.HarmonicBondForce)
        ][0]

        kj_mol_nm2 = openmm_unit.kilojoule_per_mole / openmm_unit.nanometer**2

        for idx in range(default_bond_force.getNumBonds()):
            # atom1_index, atom2_index, length (nm), k (kj/mol/nm2)
            default = default_bond_force.getBondParameters(idx)
            modified = mod_bond_force.getBondParameters(idx)
            assert default[2].value_in_unit(openmm_unit.nanometer) == pytest.approx(
                modified[2].value_in_unit(openmm_unit.nanometer)
            )
            assert default[3].value_in_unit(kj_mol_nm2) == pytest.approx(
                modified[3].value_in_unit(kj_mol_nm2)
            )

        for bond1, bond2 in zip(omm_sys_top.bonds, mod_omm_sys_top.bonds):
            # 'approx()' because https://github.com/openforcefield/openff-toolkit/issues/994
            assert bond1.fractional_bond_order == pytest.approx(
                bond2.fractional_bond_order
            )

    def test_fractional_bond_order_ignore_existing_confs(self):
        """Test that previously-defined bond orders in the topology are overwritten"""
        mol = create_ethanol()
        top = Topology.from_molecules(mol)

        mod_mol = create_ethanol()
        mod_mol.generate_conformers()
        mod_mol._conformers[0][0][0] = (
            mod_mol._conformers[0][0][0] + 1.0 * unit.angstrom
        )
        mod_mol._conformers[0][1][0] = (
            mod_mol._conformers[0][1][0] - 1.0 * unit.angstrom
        )
        mod_mol._conformers[0][2][0] = (
            mod_mol._conformers[0][2][0] + 1.0 * unit.angstrom
        )

        mod_top = Topology.from_molecules(mod_mol)

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"), xml_ff_bo
        )

        omm_system, omm_sys_top = forcefield.create_openmm_system(
            top, return_topology=True
        )
        mod_omm_system, mod_omm_sys_top = forcefield.create_openmm_system(
            mod_top, return_topology=True
        )

        # Check that the assigned bond parameters are identical for both systems
        default_bond_force = [
            f for f in omm_system.getForces() if isinstance(f, openmm.HarmonicBondForce)
        ][0]
        mod_bond_force = [
            f
            for f in mod_omm_system.getForces()
            if isinstance(f, openmm.HarmonicBondForce)
        ][0]

        for idx in range(default_bond_force.getNumBonds()):
            assert all(
                np.isclose(
                    default_value.value_in_unit(default_value.unit),
                    mod_value.value_in_unit(default_value.unit),
                )
                for default_value, mod_value in zip(
                    default_bond_force.getBondParameters(idx)[2:],
                    mod_bond_force.getBondParameters(idx)[2:],
                )
            )

        # Check that the assigned torsion parameters are identical for both systems
        default_torsion_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]
        mod_torsion_force = [
            force
            for force in mod_omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]

        for idx in range(default_torsion_force.getNumTorsions()):
            default_k = default_torsion_force.getTorsionParameters(idx)[-1]
            mod_k = mod_torsion_force.getTorsionParameters(idx)[-1]
            assert np.isclose(
                default_k.value_in_unit(default_k.unit),
                mod_k.value_in_unit(default_k.unit),
            )

        # Check that the returned topology has the correct bond order
        for bond1, bond2 in zip(omm_sys_top.bonds, mod_omm_sys_top.bonds):
            assert np.isclose(bond1.fractional_bond_order, bond2.fractional_bond_order)

    @pytest.mark.parametrize(
        (
            "get_molecule",
            "k_torsion_interpolated",
            "k_bond_interpolated",
            "length_bond_interpolated",
            "central_atoms",
        ),
        [
            (create_ethanol, 4.953856, 44375.504, 0.13770, (1, 2)),
            (create_reversed_ethanol, 4.953856, 44375.504, 0.13770, (7, 6)),
        ],
    )
    def test_fractional_bondorder_from_molecule(
        self,
        get_molecule,
        k_torsion_interpolated,
        k_bond_interpolated,
        length_bond_interpolated,
        central_atoms,
    ):
        """Test the fractional bond orders are used to interpolate k and length values as we expect

        Values for bond and torsion parameters are inferred from the bond order specified in the input molecule (1.23)

        Parameter   | SMIRNOFF unit | OpenMM unit   | param values at bond orders 1, 2  | used bond order   | value in OpenMM force
        bond k        kcal/mol/A**2   kJ/mol/nm**2    101, 123                            1.23                44375.504
        bond length   A               nm              1.4, 1.3                            1.23                0.1377
        torsion k     kcal            kj              1, 1.8                              1.23                4.953856

        See test_fractional_bondorder_calculated_openeye test_fractional_bondorder_calculated_rdkit for
        how the bond orders are obtained.
        """
        mol = get_molecule()
        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"), xml_ff_bo
        )
        topology = Topology.from_molecules(mol)

        omm_system = forcefield.create_openmm_system(
            topology,
            charge_from_molecules=[mol],
            partial_bond_orders_from_molecules=[mol],
        )

        # Verify that the assigned bond parameters were correctly interpolated
        off_bond_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.HarmonicBondForce)
        ][0]

        for idx in range(off_bond_force.getNumBonds()):
            params = off_bond_force.getBondParameters(idx)

            atom1, atom2 = params[0], params[1]
            atom1_mol, atom2_mol = central_atoms

            if ((atom1 == atom1_mol) and (atom2 == atom2_mol)) or (
                (atom1 == atom2_mol) and (atom2 == atom1_mol)
            ):
                k = params[-1]
                length = params[-2]
                assert_almost_equal(k / k.unit, k_bond_interpolated)
                assert_almost_equal(length / length.unit, length_bond_interpolated)

        # Verify that the assigned torsion parameters were correctly interpolated
        off_torsion_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if ((atom2 == atom2_mol) and (atom3 == atom3_mol)) or (
                (atom2 == atom3_mol) and (atom3 == atom2_mol)
            ):
                k = params[-1]
                assert_almost_equal(k / k.unit, k_torsion_interpolated)

    def test_fractional_bondorder_multiple_same_mol(self):
        """Check that an error is thrown when essentially the same molecule is entered more than once
        for partial_bond_orders_from_molecules"""

        mol = create_ethanol()
        mol2 = create_reversed_ethanol()

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ff_bo,
        )
        topology = Topology.from_molecules([mol, mol2])

        with pytest.raises(ValueError):
            forcefield.create_openmm_system(
                topology,
                charge_from_molecules=[mol],
                partial_bond_orders_from_molecules=[mol, mol2],
            )

    @pytest.mark.parametrize(
        ("get_molecule", "central_atoms"),
        [(create_ethanol, (1, 2)), (create_reversed_ethanol, (7, 6))],
    )
    def test_fractional_bondorder_superseded_by_standard_torsion(
        self, get_molecule, central_atoms
    ):
        """Test that matching rules are still respected with fractional bondorder parameters.

        We check here that a standard torsion parameter can still get priority on a match
        if it is further down the list.
        """

        mol = get_molecule()

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ff_torsion_bo_standard_supersede,
        )
        topology = Topology.from_molecules([mol])

        omm_system = forcefield.create_openmm_system(
            topology,
            charge_from_molecules=[mol],
            partial_bond_orders_from_molecules=[mol],
        )

        off_torsion_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if ((atom2 == atom2_mol) and (atom3 == atom3_mol)) or (
                (atom2 == atom3_mol) and (atom3 == atom2_mol)
            ):
                k = params[-1]
                assert_almost_equal(k / k.unit, 5.0208)

    @requires_rdkit
    @pytest.mark.parametrize(
        (
            "get_molecule",
            "k_torsion_interpolated",
            "k_bond_interpolated",
            "bond_length_interpolated",
            "central_atoms",
        ),
        [
            (create_ethanol, 4.953856, 42266.9635, 1.39991, (1, 2)),
            (create_reversed_ethanol, 4.953856, 42266.9635, 1.39991, (7, 6)),
        ],
    )
    def test_fractional_bondorder_calculated_rdkit(
        self,
        get_molecule,
        k_torsion_interpolated,
        k_bond_interpolated,
        bond_length_interpolated,
        central_atoms,
    ):
        """Test that torsion barrier heights interpolated from fractional bond
        orders calculated with the Amber toolkit are assigned within our
        expectations.

        This test mirrors test_fractional_bondorder but uses RDKit to re-compute the bond orders
        on input molecules (and therefore disregards the bond orders on them).

        Note that RDKitToolkitWrapper does not provide a method for assigning fractional bond orders;
        instead, AmberTools is used in the case that OpenEye is not installed.

        To obtain the bond order:
        >>> mol = create_ethanol()
        >>> AmberToolsToolkitWrapper().assign_fractional_bond_orders(mol)
        >>> mol.get_bond_between(1, 2).fractional_bond_order
        1.00093033

        """

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        )
        mol = get_molecule()

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_ff_bo,
        )
        topology = Topology.from_molecules([mol])

        omm_system, ret_top = forcefield.create_openmm_system(
            topology,
            charge_from_molecules=[mol],
            return_topology=True,
            toolkit_registry=toolkit_registry,
        )

        off_torsion_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]

        off_bond_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.HarmonicBondForce)
        ][0]

        ret_mol = list(ret_top.reference_molecules)[0]
        bond = ret_mol.get_bond_between(*central_atoms)

        # ambertools appears to yield 1.00093033 for this bond
        assert abs(bond.fractional_bond_order - 1.0) > 1.0e-6
        assert 0.95 < bond.fractional_bond_order < 1.05

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if ((atom2 == atom2_mol) and (atom3 == atom3_mol)) or (
                (atom2 == atom3_mol) and (atom3 == atom2_mol)
            ):
                k = params[-1]

                # do a hand calculation as a sanity check
                slope = (7.5312 - 4.184) / (1.8 - 1.0)
                k_torsion_interpolated_ret = (
                    slope * (bond.fractional_bond_order - 1.0) + 4.184
                )
                assert_almost_equal(k_torsion_interpolated_ret, k / k.unit, 2)

                # check that we *are not* matching the values we'd get if we
                # had offered our molecules to `partial_bond_orders_from_molecules`
                with pytest.raises(AssertionError):
                    assert_almost_equal(k / k.unit, k_torsion_interpolated)

        for idx in range(off_bond_force.getNumBonds()):
            params = off_bond_force.getBondParameters(idx)

            atom1, atom2 = params[0], params[1]
            atom1_mol, atom2_mol = central_atoms

            if atom1 == atom1_mol and atom2 == atom2_mol:
                k = params[-1]

                assert_almost_equal(k / k.unit, k_bond_interpolated, 0)

                length = params[-2]
                assert_almost_equal(length / length.unit, bond_length_interpolated, 0)

    @requires_openeye
    @pytest.mark.parametrize(
        (
            "get_molecule",
            "k_torsion_interpolated",
            "k_bond_interpolated",
            "bond_length_interpolated",
            "central_atoms",
        ),
        [
            (create_ethanol, 4.953856, 42208.540, 0.140054, (1, 2)),
            (create_reversed_ethanol, 4.953856, 42208.540, 0.140054, (7, 6)),
        ],
    )
    def test_fractional_bondorder_calculated_openeye(
        self,
        get_molecule,
        k_torsion_interpolated,
        k_bond_interpolated,
        bond_length_interpolated,
        central_atoms,
    ):
        """Test that torsion barrier heights interpolated from fractional bond
        orders calculated with the OpenEye toolkit are assigned within our
        expectations.

        This test mirrors test_fractional_bondorder but uses OpenEye to re-compute the bond orders
        on input molecules (and therefore disregards the bond orders on them).

        To obtain the bond order:
        >>> mol = create_ethanol()
        >>> OpenEyeToolkitWrapper().assign_fractional_bond_orders(mol)
        >>> mol.get_bond_between(1, 2).fractional_bond_order
        0.9945832743995565
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        mol = get_molecule()

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"), xml_ff_bo
        )
        topology = Topology.from_molecules([mol])

        omm_system, ret_top = forcefield.create_openmm_system(
            topology,
            charge_from_molecules=[mol],
            return_topology=True,
            toolkit_registry=toolkit_registry,
        )

        off_torsion_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]

        off_bond_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.HarmonicBondForce)
        ][0]

        ret_mol = list(ret_top.reference_molecules)[0]
        bond = ret_mol.get_bond_between(*central_atoms)

        # openeye toolkit appears to yield around .9945 for this bond
        assert abs(bond.fractional_bond_order - 1.0) > 1.0e-6
        assert 0.95 < bond.fractional_bond_order < 1.05

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if ((atom2 == atom2_mol) and (atom3 == atom3_mol)) or (
                (atom2 == atom3_mol) and (atom3 == atom2_mol)
            ):
                k = params[-1]

                # do a hand calculation as a sanity check
                slope = (7.5312 - 4.184) / (1.8 - 1.0)
                k_torsion_interpolated_ret = (
                    slope * (bond.fractional_bond_order - 1.0) + 4.184
                )
                assert_almost_equal(k_torsion_interpolated_ret, k / k.unit, 2)

                # check that we *are not* matching the values we'd get if we
                # had offered our molecules to `partial_bond_orders_from_molecules`
                with pytest.raises(AssertionError):
                    assert_almost_equal(k / k.unit, k_torsion_interpolated)

        for idx in range(off_bond_force.getNumBonds()):
            params = off_bond_force.getBondParameters(idx)

            atom1, atom2 = params[0], params[1]
            atom1_mol, atom2_mol = central_atoms

            if ((atom1 == atom1_mol) and (atom2 == atom2_mol)) or (
                (atom1 == atom2_mol) and (atom2 == atom1_mol)
            ):
                k = params[-1]
                assert_almost_equal(k / k.unit, k_bond_interpolated, 0)

                length = params[-2]
                assert_almost_equal(length / length.unit, bond_length_interpolated, 0)

    def test_fractional_bondorder_invalid_interpolation_method(self):
        """
        Ensure that requesting an invalid interpolation method leads to a
        FractionalBondOrderInterpolationMethodUnsupportedError
        """
        mol = create_ethanol()

        forcefield = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"), xml_ff_bo
        )
        forcefield.get_parameter_handler(
            "ProperTorsions"
        )._fractional_bondorder_interpolation = "invalid method name"
        topology = Topology.from_molecules([mol])

        # If important, this can be a custom exception instead of a verbose ValidationError
        with pytest.raises(
            ValidationError,
            match="given=invalid method name",
        ):
            forcefield.create_openmm_system(
                topology,
                charge_from_molecules=[mol],
            )


class TestForceFieldWithToolkits(_ForceFieldFixtures):
    """Test interactions between ``ForceField`` methods and wrapped toolkits."""

    # TODO: If `_toolkit_registry_manager` is made public or used for other parts of the API,
    #       these tests should be moved/adapted into more unit tests that call it directly

    def test_toolkit_registry_bogus_argument(self, force_field):
        topology = create_ethanol().to_topology()
        with pytest.raises(
            NotImplementedError,
            match="Only .*ToolkitRegistry.*ToolkitWrapper.* are supported",
        ):
            force_field.create_openmm_system(
                topology,
                toolkit_registry="rdkit",
            )

    def test_toolkit_registry_no_charge_methods(self, force_field):
        topology = create_ethanol().to_topology()
        with pytest.raises(
            ValueError, match="No registered toolkits can provide .*find_smarts_matches"
        ):
            force_field.create_openmm_system(
                topology, toolkit_registry=BuiltInToolkitWrapper()
            )

    @pytest.mark.skip(reason="Broken until Interchange supports Electrostatics 0.4")
    @requires_rdkit
    def test_toolkit_registry_bad_charge_method(self):
        topology = create_ethanol().to_topology()
        force_field = ForceField(
            get_data_file_path("test_forcefields/test_forcefield.offxml"),
            xml_charge_increment_model_ff_ethanol,
        )
        force_field.deregister_parameter_handler("ToolkitAM1BCC")
        force_field["ChargeIncrementModel"].partial_charge_method = "am1bccelf10"

        with pytest.raises(
            ValueError, match="No registered toolkits can provide .*elf10"
        ):
            force_field.create_openmm_system(
                topology,
                toolkit_registry=ToolkitRegistry(
                    [RDKitToolkitWrapper(), AmberToolsToolkitWrapper()]
                ),
            )


class TestInterchangeReturnTopology(_ForceFieldFixtures):
    # TODO: Remove these tests when `return_topology` is deprecated in version 0.12.0
    # TODO: Turn this test back on once Interchange is tagged with a pre-release supporting
    #       only Electrostatics version 0.4
    @pytest.mark.skip(reason="Interchange needs to be updated")
    def test_deprecation_warning_raised(self, force_field):
        """
        Ensure that the deprecation warning is raised when `return_topology` is
        used.
        """
        topology = create_ethanol().to_topology()
        topology.box_vectors = unit.Quantity([4, 4, 4], unit.nanometer)

        with pytest.warns(
            DeprecationWarning, match="DEPRECATED and will be removed in version 0.12.0"
        ):
            force_field.create_openmm_system(
                topology,
                return_topology=True,
            )


class TestSmirnoffVersionConverter:
    @requires_openeye_mol2
    @pytest.mark.parametrize(
        ("freesolv_id", "forcefield_version", "allow_undefined_stereo"),
        generate_freesolv_parameters_assignment_cases(),
    )
    @pytest.mark.parametrize(("spec"), ["0_1", "0_2", "0_3"])
    def test_read_smirnoff_spec_freesolv(
        self, freesolv_id, forcefield_version, allow_undefined_stereo, spec
    ):
        """
        Test reading an 0.2 smirnoff spec file, by reading an 0.1 spec representation of a set of parameters,
        and ensuring that it parameterizes molecules identically to the same FF in the most recent spec

        Regression test on parameters assignment based on the FreeSolv set used in the 0.1 paper.

        This, contrarily to the similar AlkEthOH test, checks also constraints
        and improper torsions.

        The AlkEthOH set, however, does not have impropers, which should be
        tested separately. Currently, test_freesolv_parameters_assignment
        does the job.
        """
        from openff.toolkit.tests.utils import (
            compare_system_parameters,
            get_freesolv_file_path,
        )

        mol2_file_path, xml_file_path = get_freesolv_file_path(
            freesolv_id, forcefield_version
        )

        # Load molecules.
        molecule = Molecule.from_file(
            mol2_file_path, allow_undefined_stereo=allow_undefined_stereo
        )

        # Create OpenFF System with the current toolkit.
        # forcefield_file_path = 'test_forcefields/old/smirnoff99Frosst_1_0_8_spec_0_2.offxml'
        forcefield_file_path = get_data_file_path(
            f"test_forcefields/old/test_ff_{forcefield_version}_spec_{spec}.offxml"
        )
        ff = ForceField(
            forcefield_file_path,
            get_data_file_path("test_forcefields/old/hbonds.offxml"),
        )
        ff_system = ff.create_openmm_system(molecule.to_topology())

        # Load OpenMM System created with the 0.1 version of the toolkit.
        with open(xml_file_path, "r") as f:
            xml_system = openmm.XmlSerializer.deserialize(f.read())

        # Compare parameters. We ignore the improper folds as in 0.0.3 we
        # used a six-fold implementation while we now use a three-fold way.
        # TODO: Reactivate charge comparison once we'll be able to read them from file.
        compare_system_parameters(
            ff_system,
            xml_system,
            systems_labels=("current OpenFF", "SMIRNOFF 0.0.4"),
            ignore_charges=True,
            ignore_improper_folds=True,
        )


class TestForceFieldGetPartialCharges(_ForceFieldFixtures):
    """Tests for the ForceField.get_partial_charges method."""

    @staticmethod
    def get_partial_charges_from_create_openmm_system(mol, force_field):
        """Helper method to compute partial charges from a generated openmm System."""
        system = force_field.create_openmm_system(mol.to_topology())
        nbforce = [
            f for f in system.getForces() if isinstance(f, openmm.openmm.NonbondedForce)
        ][0]

        n_particles = nbforce.getNumParticles()
        charges = [nbforce.getParticleParameters(i)[0] for i in range(n_particles)]

        return openmm_unit.Quantity(charges)

    def test_get_partial_charges(self, force_field):
        """Test that ethanol charges are computed correctly."""
        ethanol: Molecule = create_ethanol()

        ethanol_partial_charges = from_openmm(
            self.get_partial_charges_from_create_openmm_system(ethanol, force_field)
        )

        partial_charges = force_field.get_partial_charges(ethanol)

        assert (
            ethanol_partial_charges - partial_charges < 1.0e-6 * unit.elementary_charge
        ).all()
        assert partial_charges.shape == (ethanol.n_atoms,)


# TODO: test_store_cosmetic_attribs
# TODO: test_write_cosmetic_attribs
# TODO: test_store_cosmetic_elements (eg. Author)
# TODO: test_write_cosmetic_elements (eg. Author)
# TODO: add_handler_with_incompatible_kwargs (for example different scale14 vals)
# TODO: test_invalid_file_version
# TODO: test_create_gbsa
