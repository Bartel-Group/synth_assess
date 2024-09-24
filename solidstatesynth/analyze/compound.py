""" 
The purpose of this module is to quickly analyze data
    AnalyzeTarget: analyze 1 target
        get possible precursors

Goal: make it easy to analyze a given target
    could come from text-mined data, MP, or elsewhere
"""
import os
# from solidstatesynth.data.tm import (
#     get_tm_precursors,
#     get_tm_targets,
#     get_mp_cmpds,
# )
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json
from itertools import combinations
from pymatgen.ext.matproj import MPRester

DATADIR = "../data"


class AnalyzeCompound(object):
    def __init__(self, formula):
        """
        Args:
            formula (str) : formula

        Returns:
            (additionally)
            tm_precursors (list) : list of precursors in the text-mined dataset that are also in MP
            tm_targets (list) : list of targets in the text-mined dataset that are also in MP
            mp_cmpds (list) : list of compounds in MP
        """
        self.formula = CompTools(formula).clean
        self.tm_precursors = read_json(os.path.join(DATADIR, 'tm_precursors.json'))['data']
        self.tm_targets = read_json(os.path.join(DATADIR, 'tm_targets.json'))['data']
        self.mp_icsd_cmpds = read_json(os.path.join(DATADIR, 'mp_icsd_cmpds.json'))

    @property
    def in_mp(self):
        """
        Returns:
            True if the target is in MP else False
        """
        return True if self.formula in self.mp_cmpds else False
    
    @property
    def in_icsd(self):
        """
        Returns:
            True if the target is in the ICSD database
        """
        return True if self.formula in self.mp_icsd_cmpds else False
    
    @property
    def in_tm(self):
        """
        Returns:
            True if the target is in the text-mined dataset (and MP) else Fasle
        """
        return True if self.formula in self.tm_targets else False
    
    @property
    def is_NIST_gas(self):
        """
        Returns:
            True if the target is a NIST gas
        """
        return True if self.formula in ['C1O2','H2O1','O1','H1'] else False
    
    @property
    def is_oxide(self):
        """
        Returns:
            True if the target has oxygen
        """
        return True if "O" in CompTools(self.target).els else False

class AnalyzeTarget(AnalyzeCompound):
    def __init__(self, target):
        """
        Args:
            target (str) : target compound

        Returns:
            (additionally)
            tm_precursors (list) : list of precursors in the text-mined dataset that are also in MP
            tm_targets (list) : list of targets in the text-mined dataset that are also in MP
            mp_cmpds (list) : list of compounds in MP
        """
        self.target = CompTools(target).clean
        self.tm_precursors = read_json(os.path.join(DATADIR, 'tm_precursors.json'))['data']
        self.tm_targets = read_json(os.path.join(DATADIR, 'tm_targets.json'))['data']
        self.mp_icsd_cmpds = read_json(os.path.join(DATADIR, 'mp_icsd_cmpds.json'))['formulas']


    @property
    def in_tm(self):
        """
        Returns:
            True if the target is in the text-mined dataset (and MP) else Fasle
        """
        return True if self.target in self.tm_targets else False

    @property
    def chemsys(self):
        """
        Returns:
            the chemical system (str, el1-el2-...)
        """
        return CompTools(self.target).chemsys
    
    @property
    def chemsys_targets(self):
        """
        Returns:
            list of targets in the same chemical system
        """
        chemsys_targets = [t for t in self.tm_targets if CompTools(t).chemsys == self.chemsys]
        chemsys_targets.extend([t for t in self.mp_icsd_cmpds if CompTools(t).chemsys == self.chemsys])
        return list(set(chemsys_targets))

    @property
    def n_els_in_target(self):
        """
        Returns:
            the number of elements in the target (int)
        """
        return CompTools(self.target).n_els

    @property
    def flexible_els(self):
        """
        Returns a list of elements
            these elements may not be part of the target chemical system but they could be part of a precursor's chemical system

        Logic:
            if the target is an oxide, we want to consider carbonates and hydroxides as possible precursors
        """
        if self.is_oxide:
            return ["H", "C"]
        else:
            return []

    def possible_precursors(self, restrict_to_tm=True):
        """
        Args:
            restrict_to_tm (bool) : restrict to text-mined precursors if True

        Returns:
            list of possible precursors in the text-mined dataset

        Logic:
            1) precursors should be a subset (not inclusive) of the target chemical system
                e.g., if target chemsys = La-Co-O, precursor should not contain all three of these elements
            2) if the target is an oxide, we want hydroxides and carbonates as precursors

        """
        # decide whether we want to consider all precursors or just text-mined precursors
        if restrict_to_tm:
            precursors = self.tm_precursors
        else:
            precursors = self.mp_icsd_cmpds

        # what elements are in the target
        target_els = CompTools(self.target).els

        # do we have additional elements to consider
        flexible_els = self.flexible_els

        # how many elements in the target (binary, ternary, etc.)
        nary = len(target_els)

        # determine allowed chemical systems for precursors
        ## first, just consider (n-1)ary systems that are subsets of the target chemical system
        allowed_els = []
        for n in range(1, nary):
            allowed_els.extend(list(combinations(target_els, n)))

        ## now incorporating our "flexible elements" (basically adding carbonates and hydroxides)
        flexible_els = self.flexible_els
        new_allowed_els = []
        for el in flexible_els:
            for el_combo in allowed_els:
                el_combo = list(el_combo)
                if ("O" in el_combo) and (el not in el_combo) and (len(el_combo) > 1):
                    el_combo.append(el)
                    el_combo = tuple(sorted(el_combo))
                    new_allowed_els.append(el_combo)
        allowed_els = set(allowed_els + new_allowed_els)
        allowed_el_list = target_els
        if 'O' in target_els:
            allowed_el_list.extend(['C','H'])
        # filter our big list of precursors down to those that we deemed "possible"
        precursors = [p for p in precursors if tuple(CompTools(p).els) in allowed_els]
        precursors.extend(allowed_el_list)
        return list(set(precursors))


def check():
    target = "BaTiO3"
    at = AnalyzeTarget(target)
    possible_precursors = at.possible_precursors(restrict_to_tm=True)
    print("Target = %s" % target)
    print("Possible precursors = %s" % possible_precursors)
    return at


def main():
    check()
    return


if __name__ == "__main__":
    at = main()
