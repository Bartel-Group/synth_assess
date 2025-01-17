from pymatgen.core.composition import Composition
from pydmclab.core.comp import CompTools
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.experimental import ExperimentalReferenceEntry
from solidstatesynth.extract.mp import get_gases_data



DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"


class Gibbs:
    def __init__(
        self,
        formula,
        solids_data,
        gases_data=get_gases_data(),
        temperature=300,
        use_carbonate_correction=True,
        entry_id = None
    ):
        """
        Args:
            formula (str): formula of the target compound
            solids_data (dict): {formula (str, clean formula) : {< dict with mp data >}}
            gases_data (dict): dictionary of dictionaries containing data to use for gases having that formula (from experiment)
                {formula (str) : {temperature (int) : free energy (float, eV/f.u.)}}
                Refer to McDermott's rxn_network/entries/experimental.py for the structure of any desired
                experimental data
            temperature (int): temperature of interest
            use_carbonate_correction (bool) : whether to use carbonate correction for carbonates

        """
        # clean the formula coming in
        # self.entry = entry
        self.formula = CompTools(formula).clean
        self.temperature = temperature
        self.use_carbonate_correction = use_carbonate_correction
        self.entry_id = entry_id

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

            # # find formula in solids_data
            if formula not in solids_data:
                raise ValueError("No data for the target compound: %s" % formula)

            self.compound_data = solids_data[formula]
            # print(self.compound_data)

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
            carbonate correction for the compound 0.83 * n_{CO3} groups per total number of atoms (eV/atom)
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
            Note that energy is in eV, not eV/atom
        """
        entry = self.compound_data
        vol_per_at = entry["volume"] /entry["nsites"]
        Ef_per_at = (
            entry["formation_energy_per_atom"] + self.carbonate_correction
        )
        return GibbsComputedEntry(
            volume_per_atom=vol_per_at,
            formation_energy_per_atom=Ef_per_at,
            composition=Composition(self.formula),
            temperature=self.temperature,
            entry_id= self.entry_id
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


class FormulaChecker:
    def __init__(
        self, formula,chemsys_els,extend_with_hydroxides=True, extend_with_carbonates=True
    ):
        """
        Args:
            formula (str): formula to check
            chemsys [(str)]: [el1,el2,el3...]
            extend_with_hydroxides (bool): whether to include hydroxides as additional compounds in *chemsys*
            extend_with_carbonates (bool): whether to include carbonates as additional compounds in *chemsys*
        NOTE that extensions will occur only for oxide chemistries
        """
        self.formula = formula
        # self.chemsys = chemsys
        self.extend_with_hydroxides = extend_with_hydroxides
        self.extend_with_carbonates = extend_with_carbonates
        self.chemsys_els = chemsys_els

    @property
    def els_to_check(self):
        """
        We will check the formula in question for its relevance to the chemical system of interest
        """
        return CompTools(self.formula).els

    @property
    def same_chemsys(self):
        return self.els_to_check == self.chemsys_els

    @property
    def sub_chemsys(self):
        return set(self.els_to_check).issubset(set(self.chemsys_els))

    @property
    def is_hydroxide(self):
        formula = self.formula
        n_H, n_O = CompTools(formula).stoich("H"), CompTools(formula).stoich("O")
        if n_H and (n_O / n_H == 1):
            return True
        return False

    @property
    def is_carbonate(self):
        formula = self.formula
        n_C, n_O = CompTools(formula).stoich("C"), CompTools(formula).stoich("O")
        if n_C and (n_O / n_C == 3):
            return True
        return False

    @property
    def is_relevant(self):
        """
        Returns:
            True if formula of interest is in the chemical system or subsystem
            This allows for carbonates and hydroxides if the system has oxygen (and they are deemed allowed)
        """
        if self.same_chemsys:
            return True
        if self.sub_chemsys:
            return True
        els_to_check = self.els_to_check
        els = self.chemsys_els
        if "O" in els:
            if self.extend_with_carbonates and self.is_carbonate:
                non_CO_els = [el for el in els_to_check if el not in ["C", "O"]]
                if set(non_CO_els).issubset(set(els)):
                    return True
            if self.extend_with_hydroxides and self.is_hydroxide:
                non_OH_els = [el for el in els_to_check if el not in ["O", "H"]]
                if set(non_OH_els).issubset(set(els)):
                    return True
        return False


class GibbsSet:
    def __init__(
        self,
        chemsys_els,
        solids_data,
        gases_data=get_gases_data(),
        temperature=300,
        use_carbonate_correction=True,
        extend_with_hydroxides=True,
        extend_with_carbonates=True,
        with_theoretical=True,
        stability_threshold=0.5,
        entry_id_dict = {},
        include_only_these_formulas=[],
        exclude_these_formulas=[],
        add_these_formulas=[],
    ):
        """
         Args:
            chemsys_els [str]: ['el1', 'el2', 'el3'...]
            solids_data (dict): {formula (str, clean formula) : {< dict with data >}}
                < dict with data > must have these keys:
                    {'formula' : <clean formula (CompTools(formula).clean)> (str),
                     'volume' : <volume of calculated structure> (float),
                     'nsites' : <number of sites in calculated structure> (int),
                     'formation_energy_per_atom' : <formation enthalpy per atom at 0 K> (float),
                     'energy_above_hull' : <energy above hull> (float, eV/at)}
                notes:
                    - this should be all the ground-states you want to consider
                    - if solids data is not given, it will be extracted from the MP database (using only GS)
                    - the stability threshold is used only for default solids data. Users should impose their own
                    stability thresholds if they provide their own solids data
            gases_data (dict): dictionary of dictionaries containing data to use for gases having that formula (from experiment)
                {formula (str) : {temperature (int) : free energy (float, eV/f.u.)}}
            temperature (int): temperature of interest
            use_carbonate_correction (bool) : whether to use carbonate correction for carbonates
            extend_with_hydroxides (bool): whether to include hydroxides as additional compounds in *chemsys*
            extend_with_carbonates (bool): whether to include carbonates as additional compounds in *chemsys*
            stability_threshold (float): maximum energy above hull for a compound to be considered as "relevant"
            include_only_these_formulas (list): list of formulas to include (if you want to specify them)
            exclude_these_formulas (list): list of formulas to exclude (if you don't want them included)
            add_these_formulas (list): list of formulas to add (if you want to augment list with some)

        """
        self.solids_data = solids_data
        self.gases_data = gases_data
        self.temperature = temperature
        self.use_carbonate_correction = use_carbonate_correction
        self.extend_with_hydroxides = extend_with_hydroxides
        self.extend_with_carbonates = extend_with_carbonates
        self.stability_threshold = stability_threshold
        self.include_only_these_formulas = include_only_these_formulas
        self.exclude_these_formulas = exclude_these_formulas
        self.add_these_formulas = add_these_formulas
        self.els = chemsys_els
        self.entry_id_dict = entry_id_dict


    @property
    def formulas(self):
        """
        Returns:
            list of formulas to include in the entry set
        """
        if self.include_only_these_formulas:
            return self.include_only_these_formulas

        solids_data = self.solids_data
        all_formulas = list(solids_data.keys())
        chemsys = self.els
        extend_with_hydroxides, extend_with_carbonates = (
            self.extend_with_hydroxides,
            self.extend_with_carbonates,
        )

        formulas = [
            f
            for f in all_formulas
            if FormulaChecker(
                formula=f,
                chemsys_els=chemsys,
                extend_with_carbonates=extend_with_carbonates,
                extend_with_hydroxides=extend_with_hydroxides,
            ).is_relevant
        ]
        if self.extend_with_carbonates:
            formulas.append('C1O2')
        if self.extend_with_hydroxides:
            formulas.append('H2O1')

        if self.exclude_these_formulas:
            formulas = [f for f in formulas if f not in self.exclude_these_formulas]
        if self.add_these_formulas:
            formulas.extend(self.add_these_formulas)
        return formulas

    @property
    def entries(self):
        """
        Returns:
            List of GibbsEntry objects for all formulas of interest
        """
        formulas = self.formulas
        id_dict = self.entry_id_dict
        for formula in formulas:
            f = CompTools(formula).clean
            if f not in id_dict:
                id_dict[f] = None
        solids_data = self.solids_data
        gases_data = self.gases_data
        temperature = self.temperature
        use_carbonate_correction = self.use_carbonate_correction


        return [
            Gibbs(
                formula=f,
                solids_data=solids_data,
                gases_data=gases_data,
                temperature=temperature,
                use_carbonate_correction=use_carbonate_correction,
                entry_id= id_dict[CompTools(f).clean]
            ).entry
            for f in formulas
        ]

    @property
    def entry_set(self):
        """
        Returns:
            GibbsEntrySet for the chemical space of interest
        """
        return GibbsEntrySet(self.entries)


# def main():
#     # g = Gibbs("Fe2O3", temperature=1317)
#     # e = g.entry
#     # gs = GibbsSet(chemsys_els=['Fe','O'], temperature=1317)
#     # es = gs.entries
#     # gse = gs.entry_set
#     # return gs, es, gse


# if __name__ == "__main__":
#     return
#     # gs, es, gse = main()
