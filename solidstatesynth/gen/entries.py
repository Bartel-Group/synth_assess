from pydmclab.utils.handy import read_json, write_json
import os
import math
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry, ConstantEnergyAdjustment
from pymatgen.ext.matproj import MPRester
from pymatgen.entries import computed_entries
from pydmclab.core.comp import CompTools
from pydmclab.data.thermochem import gas_thermo_data
from solidstatesynth.extract.mp import get_gases_data, get_MP_data
from rxn_network.entries.experimental import ExperimentalReferenceEntry
from itertools import combinations
import numpy as np


DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"


class Gibbs:
    def __init__(
        self,
        formula,
        solids_data=get_MP_data()["data"],
        gases_data=get_gases_data(),
        temperature=300,
        use_carbonate_correction=True,
    ):
        """
        Args:
            formula (str): formula of the target compound
            solids_data (list): list of dictionaries containing data to use for solids having that formula
                each dict must have the following keys:
                    {'formula' : <clean formula (CompTools(formula).clean)> (str),
                     'volume' : <volume of calculated structure> (float),
                     'nsites' : <number of sites in calculated structure> (int),
                     'formation_energy_per_atom' : <formation enthalpy per atom at 0 K> (float)}
                ideally:
                    - only ground state entries
                    - all formulas are clean
                note:
                    - this might not be the best default dictionary
            gases_data (dict): dictionary of dictionaries containing data to use for gases having that formula (from experiment)
                {formula (str) : {temperature (int) : free energy (float, eV/f.u.)}}
            temperature (int): temperature of interest
            use_carbonate_correction (bool) : whether to use carbonate correction for carbonates

        """
        # clean the formula coming in
        formula = CompTools(formula).clean
        self.formula = formula
        self.temperature = temperature
        self.use_carbonate_correction = use_carbonate_correction

        # retrieve gases data for this formula if its a gas
        gases_data = {
            f: gases_data[f] for f in gases_data if CompTools(f).clean == formula
        }

        # if gases_data isn't empty, we have a gas
        if gases_data:
            self.is_gas = True
            # set the compound data to the gases data for this formula
            self.compound_data = gases_data
        else:
            # if not a gas, look in solids_data
            self.is_gas = False

            # find formula in solids_data (ideally this is just one entry)
            ## NOTE: make sure dictionary passed here is clean then remove the cleaning done w/in this list
            solids_data = [
                e for e in solids_data if CompTools(e["formula"]).clean == formula
            ]

            # if no data, then we don't have data for this formula..
            if not solids_data:
                raise ValueError("No data for the target compound: %s" % formula)

            # if multiple entries, then something went awry b/c solid_data should just be ground-states
            if len(solids_data) > 1:
                print(
                    "WARNING: %i entries for the target compound: %s; filter ground-states before passing here; using first entry"
                    % (len(solids_data), formula)
                )
            # compound data is the only/first entry
            self.compound_data = solids_data[0]

    @property
    def is_carbonate(self):
        """
        Returns:
            True if the formula is a carbonate (C:O = 1:3)
        """
        formula = self.formula
        n_C, n_O = CompTools(formula).stoich("C"), CompTools(formula).stoich("O")
        if n_C and (n_O / n_C == 3):
            return True
        return False

    @property
    def carbonate_correction(self):
        """
        Returns:
            carbonate correction for the compound 0.83 * n_{CO3} groups (eV/atom)
        """
        if not self.use_carbonate_correction or not self.is_carbonate:
            return 0
        formula = self.formula
        return 0.830 * CompTools(formula).stoich("C") / CompTools(formula).n_atoms

    @property
    def gas_entry(self):
        """
        Returns:
            ExperimentalReferenceEntry for a gas at a temperature of interest


        Returns entry set for a gas at a temperature of interest. Note that the structure
        of this entry does not allow for modification at different temperatures (because
        temperatures are experimentally determined and G(T) cannot be extrapolated). The data
        used here in references is NIST experimental data for H2O and CO2 at different temperatures.
        NOTE THAT the functionality of this code requires that H2O and CO2 be written as such rather
        than clean formulas H2O1 and C1O2. The get_gases_data function in extract/mp.py is written
        such that keys are in this format as well. Note that temperature keys are ints (not strings)
        """
        formula = self.formula
        refs = self.compound_data
        ExperimentalReferenceEntry.REFERENCES = refs
        exp_entry = ExperimentalReferenceEntry(
            composition=Composition(formula), temperature=self.temperature
        )
        # print(exp_entry)
        return exp_entry

    @property
    def solid_entry(self):
        """
        Returns:
            GibbsComputedEntry for solid at a temperature of interest
        """
        compound_data = self.compound_data
        formula = self.formula
        # print(formula)
        vol_per_at = compound_data["volume"] / compound_data["nsites"]
        Ef_per_at = (
            compound_data["formation_energy_per_atom"] + self.carbonate_correction
        )
        return GibbsComputedEntry(
            volume_per_atom=vol_per_at,
            formation_energy_per_atom=Ef_per_at,
            composition=Composition(formula),
            temperature=self.temperature,
        )

    @property
    def entry(self):
        """
        Returns:
            GibbsComputedEntry if the target is a solid, ExperimentalReferenceEntry if the target is a gas
        """
        if self.is_gas:
            return self.gas_entry
        return self.solid_entry

    @property
    def dGf(self):
        """
        Returns:
            Gibbs free energy of formation per atom at a temperature of interest
        """
        if self.is_gas:
            return self.entry.energy_per_atom
        return self.entry.formation_energy_per_atom


class BuildGibbsEntrySet:
    """
    Builds an EntrySet for formulas associated with a target of interest. This entry set is comprised
    of GibbsComputedEntries for ground state polymorphs and ExperimentalReferenceEntries for gases.
    Note that initialization requires a determination of whether hypothetical compounds may be accounted for
    (with theoretical) and a stability filter (for energy above the hull)Depending on the with_theoretical initialization,
    MP with or without theoretical compounds may be employed

    """

    def __init__(
        self,
        els,
        temperature,
        with_theoretical=True,
        stability_filter=0.05,
        formulas=None,
        energies_at_300=True,
    ):

        self.with_theoretical = with_theoretical
        self.temperature = temperature
        self.stability_filter = stability_filter
        self.energies_at_300 = energies_at_300
        if formulas:
            self.formulas = formulas
        else:
            if not with_theoretical:
                self.data = read_json(
                    os.path.join(DATADIR, "241002_mp_experimental.json")
                )["data"]
            else:
                self.data = read_json(os.path.join(DATADIR, "241002_mp_gd.json"))[
                    "data"
                ]
            self.formulas = [entry["formula"] for entry in self.data]
        self.target_els = els

    def is_competing_formula(self, mp_formula, target_els):
        """
        Returns:
            True if the formula is in the chemical system (or sub chemical system) of the target
        """
        formula_els = CompTools(mp_formula).els
        if list(set(formula_els)) == list(set(["C", "H"])):
            return False
        if list(set(formula_els)) == list(set(["C", "H", "O"])):
            return False
        if all([el in target_els for el in formula_els]):
            return True
        allowed_els = []
        new_allowed_els = []
        for n in range(1, len(target_els)):
            allowed_els.extend(list(combinations(target_els, n)))
        if "O" in target_els:
            flexible_els = ["C", "H"]
        for el in flexible_els:
            for el_combo in allowed_els:
                el_combo = list(el_combo)
                # print(el_combo)
                if ("O" in el_combo) and (el not in el_combo) and (len(el_combo) > 1):
                    el_combo.append(el)
                    el_combo = tuple(sorted(el_combo))
                    new_allowed_els.append(el_combo)
        # print(new_allowed_els)
        allowed_els = set(allowed_els + new_allowed_els)
        # print(allowed_els)
        # filter our big list of precursors down to those that we deemed "possible"
        # precursors = [p for p in precursors if tuple(CompTools(p).els) in allowed_els]

        if tuple(CompTools(mp_formula).els) in allowed_els:
            return True
        return False

    def chemsys_competing_formulas(self):
        """
        Returns a list of all relevant competing formulas from Materials Project for the target (as defined above)
        Note that depending on which version of MP is used, this may account only for ground state experimental formulas or all ground state formulas.
        """
        target_els = self.target_els
        mp_data = self.data
        formulas_new = []
        formula_strings = []
        for entry in mp_data:
            formula_str = entry["formula"]
            if formula_str not in formula_strings:
                if self.is_competing_formula(formula_str, target_els):
                    formulas_new.append(entry)
                    formula_strings.append(formula_str)
        return formulas_new

    def build_entry_set(self):
        """
        Builds an entry set for the chemical space of interest. This entry set is comprised of GibbsComputedEntries for
        ground state polymorphs and ExperimentalReferenceEntries for gases. Temperature can be specified.
        """
        at_300 = self.energies_at_300
        temperature = self.temperature
        competing_formulas = self.chemsys_competing_formulas()
        print("competing formulas found")
        data = self.data
        entries = []
        for formula_data in competing_formulas:
            formula_str = formula_data["formula"]
            print(formula_data["material_id"])
            print("building entry")
            GibbsEntry = BuildGibbsEntry(formula_str, data)
            if GibbsEntry.is_NIST_gas:
                entry = GibbsEntry.gas_ExperimentalReferenceEntry_at_temp(temperature)
            else:
                entry = GibbsEntry.ground_GibbsComputedEntry_at_temp(
                    formula_data, temperature, energies_at_300=at_300
                )
            # print(entry)
            entries.append(entry)
            # print(entries)
            print("entry added")
        print("entries", entries)
        return GibbsEntrySet(entries)


def main():
    g = Gibbs("H2O", temperature=317)
    e = g.entry
    return g, e


if __name__ == "__main__":
    g, e = main()
