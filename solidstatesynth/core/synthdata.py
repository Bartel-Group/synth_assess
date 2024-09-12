from pydmclab.core.comp import CompTools
from pydmclab.core.energies import ReactionEnergy, ChemPots
from pydmclab.data.thermochem import gas_thermo_data
from pydmclab.core.query import MPQuery
from pymatgen.analysis import interface_reactions
from pymatgen.analysis.phase_diagram import PhaseDiagram
# from rxn_network.costs import calculators
from solidstatesynth.utils import *
from solidstatesynth.competitions import *

class PositiveEntry(object):
    def __init__(self, entry, mp_data, gases_data, new_db, interpolation_data = None,):
        """
        This class is built to navigate one entry from the simplified text-mined database (new_db).
        Functions 'has_data', 'precursors','target','temperature', 'positive_reaction','chemsys', 'dG_rxn', 
        'dG_rxn_0', and 'interpolated' extract existing information from the entry. Meanwhile, functions 'has_gas', 
        'max_dG_target', 'limiting reagent', are used to extrapolate additional information
        """
        self.entry = entry
        self.mp_data = mp_data
        self.gases_data = gases_data
        self.new_db = new_db
        self.interpolation_data = interpolation_data
        # self.restructured_db = get_restructured_db()

    @property
    def has_data(self):
        """
        identifies whether a given entry is empty, or has reaction information.
        """
        if "common" in self.entry:
            return True
        else:
            return False

    @property
    def precursors(self):
        entry = self.entry
        return entry["common"]["precursors"]

    @property
    def target(self):
        entry = self.entry
        return entry["positive"]["target"]

    @property
    def temperature(self):
        entry = self.entry
        return entry["common"]["temperature"]

    @property
    def positive_reaction(self):
        entry = self.entry
        return entry["positive"]["reaction"]


    @property
    def chemsys(self):
        entry = self.entry
        return entry["common"]["chemsys"]

    @property
    def n_els_target(self):
        entry = self.entry
        return CompTools(entry["positive"]["target"]).n_els

    @property
    def interpolated(self):
        entry = self.entry
        interpolated = entry["common"]["interpolated"]
        if len(interpolated) > 0:
            return True
        else:
            return False

    @property
    def dG_rxn(self):
        entry = self.entry
        return entry["positive"]["dG_rxn"]
    
    @property
    def dG_rxn_0(self):
        """Note that dG_rxn_0 represents the reaction energy at 0K
        """
        entry = self.entry
        return entry["positive"]["dG_rxn_0"]
    
    @property
    def has_gas(self):
        """This function identifies common gasses present in the reaction string of a function. The energetics
        of reactions with gas evolution are more complex and therefore may require different treatment, or incur greater
        error.
        """
        entry = self.entry
        reaction = entry["positive"]["reaction"]
        gases = list(gas_thermo_data().keys())
        for gas in gases:
            if CompTools(gas).clean in reaction:
                return True
        alt_gases = ['CO2','H2O','O2']
        for gas in alt_gases:
            if gas in reaction:
                return True
        return False
    
    @property
    def amts_vars(self):
        entry = self.entry
        return entry['positive']['amts_vars'], entry['positive']['amts_vars_target']
    
    def dG_rxn_calc(self, precursor = None, target=None, corrected = False, amts_vars = None, amts_vars_target = None):
        """
        This function is used to calculated the reaction energy associated with a given entry. Note that the target and 
        precursor variables default to None-- when they are undefined, the target and precursors listed in the entry are used. 
        However, alternative precursors/targets may also be defined and the reaction energy associated with accessing these targets 
        from the entry precursors and temperature may be determined. The function can therefore be employed to identify the max-dG 
        target.
        """
        mp_data = self.mp_data
        interpolation_data = self.interpolation_data
        gases_data = self.gases_data
        temperature = self.temperature
        entry = self.entry
        environment = entry['common']['environment']
        if not target:
            target = self.target
        if not precursor:
            precursor = self.precursors
        if len(precursor) == 1:
            rxn_dict = {"reactants":[precursor], "products":[str(target)]}
        if corrected:
            amts_vars, amts_vars_target = self.amts_vars
        rxn_dict = {"reactants":[prec for prec in precursor if prec], "products":[str(target)]}
        dG_rxn,rxn = get_dGrxn_at_T(rxn_dict, temperature,mp_data, gases_data, 
                                    interpolation_data, corrected, environment,
                                    amts_vars, amts_vars_target, target)
        return dG_rxn,rxn
    
    # @property
    def max_dG_target(self, corrected = False):
        """
        This function is used to identify the 'max_dG_target' associated with a given entry. The max-dG target is
        that which forms with the lowest (most negative) reaction energy from the listed precursors at the listed temperature.
        The max-dG target may be the observed target. This function returns (the max-dG target, the associated reaction energy).
        """
        mp_data = self.mp_data
        gases_data = self.gases_data
        new_db = self.new_db
        chemsys = self.chemsys
        ce = ChemsysEntry(chemsys,mp_data, gases_data,new_db)
        chemsys_targets = ce.chemsys_targets
        max_dG = 0
        max_dG_target = None
        for t in chemsys_targets:
            dG_rxn,rxn = self.dG_rxn_calc(target=t, corrected=corrected)
            if not dG_rxn:
                continue
            elif dG_rxn:
                if dG_rxn < max_dG:
                    max_dG = dG_rxn
                    max_dG_target = t
        if max_dG_target == None:
            return None, None
        return max_dG_target, max_dG
    
    # @property
    def limiting_reagent(self, corrected = False):
        """
        This function is employed in identifying the limiting reagent in the preliminary reaction
        from precursors to the max-dG target. This reagent will be consumed completely in the first
        reaction and therefore the observed target is formed only from other precursors and the max-dG target.
        This function returns the limiting precursor (a string).
        """
        max_prec = None
        precursors = [CompTools(prec).clean for prec in self.precursors]
        max_dG_target,max_dG = self.max_dG_target(corrected = corrected)
        target = self.target
        max_dG_coefs = get_balanced_reaction_coefficients(precursors,str(max_dG_target))
        rxn_coefs = get_balanced_reaction_coefficients(precursors,target)
        if not max_dG_coefs:
            print(self.entry['common']['doi'])
            return max_prec
        if not rxn_coefs:
            return None
        prec_ratio_dict = {prec:max_dG_coefs[prec]/rxn_coefs[prec] for prec in precursors}
        max_prec_amt = max(prec_ratio_dict.values())
        max_prec_list = [prec for prec,amt in prec_ratio_dict.items() if amt == max_prec_amt]
        if max_prec_list:
            return max_prec_list[0]
        else:
            return None
        

    def driving_force(self, corrected=False):
        """
        This function determines the reaction energy associated with forming the observed target of 
        the entry from the remaining precursor(s) and the max-dG target (also known as the "driving force").
        Note that the function returns the reaction energy and the reaction string associate with the final reaction (?)
        """
        entry = self.entry
        max_dG_target, max_dG = self.max_dG_target(corrected = corrected)
        if not max_dG_target:
            return None, None
        if entry['positive']['target'] == str(max_dG_target):
            return 0, None
        else:
            limiting_reagent = self.limiting_reagent
            prec_new = [prec for prec in self.precursors if prec != limiting_reagent]
            prec_new.append(max_dG_target)
        dG_rxn,rxn = self.dG_rxn_calc(precursor=prec_new, corrected=corrected)
        return dG_rxn, rxn
    
    # @property
    def competitions(self, use_prec = True, open = True):
        entry = self.entry
        target = str(CompTools(entry['positive']['target']).clean)
        temperature = entry['common']['temperature']
        precursors = entry['common']['precursors']
        environment = entry['common']['environment']
        if use_prec:
            prec = precursors
        else:
            prec = None
        if not temperature:
            temperature = 1073
        # print('target:',target)
        # print('precursors:',prec)
        # print('temperature:',temperature)
        comp_data = get_competition_data(target, prec, temperature, environment, open)
        # if len(comp_data)>1:
        #     comp_data = get_target_reaction(precursors,comp_data)
        return (comp_data[0]['c1'], comp_data[0]['c2'])

class ChemsysEntry(object):
    def __init__(self, chemsys, mp_data, gases_data, new_db, interpolation_data=None):
        """
        This class takes in a given chemical system of the form 'A-B-C', where elements are organized in
        alphabetical order (by pydmclab convention). This class enables the identification of entries in the
        simplified database associated with this chemical sysem. This allows for the comparison of different targets
        in the same chemical space. Functions for determining the limiting reagent, driving force associated with a reaction
        """
        self.chemsys = chemsys
        self.mp_data = mp_data
        self.interpolation_data = interpolation_data
        self.gases_data = gases_data
        self.new_db = new_db

    @property
    def chemsys_entries(self):
        new_db = self.new_db
        chemsys = self.chemsys
        entries = []
        for entry in new_db:
            if entry["common"]["chemsys"] == chemsys:
                entries.append(entry)
        return entries
    
    @property
    def chemsys_targets(self):
        chemsys_entries = self.chemsys_entries
        target_list = []
        for entry in chemsys_entries:
            target = entry['positive']['target']
            if target not in target_list:
                target_list.append(target)
        return target_list

    @property
    def chemsys_precursors(self):
        mp_data = self.mp_data
        gases_data = self.gases_data
        new_db = self.new_db
        chemsys_entries = self.chemsys_entries
        precursors = []
        for entry in chemsys_entries:
            precursors.append(str(PositiveEntry(entry,mp_data,gases_data,new_db).precursors))
        return list(set(precursors))
    
    @property
    def precursor_targets(self, precursor):
        chemsys_entries = self.chemsys_entries
        mp_data = self.mp_data
        interpolation_data = self.interpolation_data 
        gases_data = self.gases_data
        textmined_data = self.textmined_data
        new_db = self.new_db
        targets = []
        for entry in chemsys_entries:
            pos_entry = PositiveEntry(entry,mp_data,gases_data,new_db)
            if precursor == str(pos_entry.precursors):
                targets.append(pos_entry.target)
        return list(set(targets))
    
    



def main():
    return


if __name__ == "__main__":
    main()