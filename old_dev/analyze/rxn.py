""" 
The purpose of this module is to quickly analyze a reaction
    AnalyzeRxn: analyze 1 reaction
        get information that may inform corrections
            gas concentrations
            presence of carbonate

Goal: make it easy to filter as needed

"""

# from solidstatesynth.dev.extract.tm import get_tm_precursors, get_tm_targets
import os
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json,write_json
from pymatgen.core.composition import Composition
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from solidstatesynth.analyze.compound import AnalyzeCompound, AnalyzeChemsys
from solidstatesynth.analyze.thermo import AnalyzeThermo
from solidstatesynth.core.utils import get_balanced_reaction_coefficients
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.reactions.basic import BasicReaction
DATADIR = "../data"
DATADIR_cemsbartel = "/Volumes/cems_bartel/projects/negative-examples/data"

class AnalyzeSynthesisRecipe(object):
    def __init__(self, precursors, target, temperature, environment):
        """
        Args:
            precursors (list):
                list of precursors in the reaction
                for the purposes of reaction generation, precusors may include all possible participants

            target (str):
                target compound in the reaction

            temperature (float):
                temperature of the reaction

            environment (str):
                atmosphere of the reaction

        Returns:
            (additionally)
            tm_precursors (list) : list of precursors in the text-mined dataset that are also in MP
            experimental materials project formulas (dict) : dict with targets as keys for all targets also found in ICSD


        """
        self.precursors = [CompTools(p).clean for p in precursors]
        self.target = CompTools(target).clean
        self.temperature = temperature
        self.environment = environment
        self.tm_precursors = read_json(os.path.join(DATADIR, 'tm_precursors.json'))['data']
        self.mp_experimental = read_json(os.path.join(DATADIR_cemsbartel, '241002_mp_experimental.json'))['data']
        # self.tm_targets = get_tm_targets(None, None)

    @property
    def precursors_in_mp_and_tm(self):
        """
        Returns:
            True if all precursors are in the text-mined dataset and MP
        """
        curr_precursors = self.precursors
        tm_precursors = self.tm_precursors
        return True if all([p in tm_precursors for p in curr_precursors]) else False

    @property
    def system_treatment(self):
        """
        Returns:
            the treatment of the system (str)

        Logic:
            we consider systems open to oxygen (return: 'open') if
                the target is an oxide and the atmosphere is air, oxygen, or not specified

            we consider systems closed (return: 'closed') if:
                the taget is not an oxide or the atmosphere is something else (e.g., inert, H2, etc)
        """
        environment = self.environment
        if AnalyzeCompound(self.target).is_oxide and (
            environment in ["air", "oxygen", None]
        ):
            return "open"
        else:
            return "closed"

    @property
    def gases_that_get_activity_correction(self):
        """
        Returns:
            list of gases where we want to correct their free energy
                to account for their increased activity at low concentrations

        Logic:
            O2, CO2, and H2O are the most common gases involved in reactions
                and we generally know their concentrations based on the atmosphere reasonably well
        """
        return ["O2", "C1O2", "H2O1"]

    @property
    def default_low_gas_concentration(self):
        """
        Returns:
            the default low gas concentration (float)

        Logic:
            for gases that are involved in the reaction but not inherent to the atmosphere
                we assume they are present in a low concentration (p_gas = 1e-5)
        """
        return 1e-5

    @property
    def gas_concentrations(self):
        """
        Returns:
            {gas (str) : concentration (float, partial pressure in atmospheres)}

        Logic:
            for air or unspecified atmospheres:
                pO2 = 21%, pCO2 = 400 ppm, pH2O = 1% (rough estimate)
        """
        environment = self.environment
        gases = self.gases_that_get_activity_correction

        if environment in ["air", None]:
            conc = {"O2": 0.21, "C1O2": 0.0004, "H2O1": 0.01}
        elif environment == "oxygen":
            conc = {"O2": 1.0}
        else:
            conc = {}

        for g in gases:
            if g not in conc:
            # if g not in gases:
                conc[g] = self.default_low_gas_concentration
        return conc
    


    @property
    def has_carbonate_precursor(self):
        """
        Returns:
            True if any of the precursors are carbonates

        Logic:
            may require special correction
        """
        precursors = self.precursors
        for p in precursors:
            n_carbon = CompTools(p).stoich("C")
            n_oxygen = CompTools(p).stoich("O")
            if n_carbon and n_oxygen and (n_oxygen / n_carbon == 3):
                return True
        return False
    
    @property
    def balanceable(self):
        """
        Returns:
            True if the reaction is balanceable
        """
        precursors = [Composition(p) for p in self.precursors]
        targets = [self.target,'C1O2','H2O1','O2']
        targets = [Composition(t) for t in targets]
        try:
            rxn = Reaction(precursors, targets)
            return True
        except ReactionError as e:
            if str(e) == "Reaction cannot be balanced.":
                return False
    

# def check():
#     precursors = ["BaCO3", "TiO2"]
#     target = "BaTiO3"
#     rxn = "BaCO3 + TiO2 -> BaTiO3 + CO2"
#     environment = "air"
#     temperature = 1000

#     ar = AnalyzeRxn(
#         precursors=precursors,
#         target=target,
#         rxn=rxn,
#         temperature=temperature,
#         environment=environment,
#     )

#     print("has carbonate precursor:", ar.has_carbonate_precursor)
#     print("gas concentrations:", ar.gas_concentrations)
#     print("system treatment:", ar.system_treatment)
#     print("precursors in mp and tm:", ar.precursors_in_mp_and_tm)
#     return ar


def main():

    # ar = check()
    return 


if __name__ == "__main__":
   main()
