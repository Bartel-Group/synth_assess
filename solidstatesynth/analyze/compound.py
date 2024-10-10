""" 
The purpose of this module is to quickly analyze data
    AnalyzeTarget: analyze 1 target
        get possible precursors

Goal: make it easy to analyze a given target
    could come from text-mined data, MP, or elsewhere
"""
import os
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json
from itertools import combinations
from pymatgen.ext.matproj import MPRester
from solidstatesynth.gen.build_entry import BuildGibbsEntrySet

DATADIR = "../data"
DATADIR_cemsbartel = "/Volumes/cems_bartel/projects/negative-examples/data"


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
        # any checks on the formula should be compared to the "clean" formula as the input formula is always cleaned
        self.tm_precursors = read_json(os.path.join(DATADIR, 'tm_precursors.json'))['data']
        self.tm_targets = read_json(os.path.join(DATADIR, 'tm_targets.json'))['data']
        self.gd_MP = read_json(os.path.join(DATADIR_cemsbartel, '241002_mp_gd.json'))['data']
        self.mp_experimental = read_json(os.path.join(DATADIR_cemsbartel, '241002_mp_experimental.json'))['data']


    @property
    def in_mp(self):
        """
        Returns:
            True if the target is in MP else False
        """
        return True if self.formula in self.gd_MP else False
    
    @property
    def in_icsd(self):
        """
        Returns:
            True if the target is in the ICSD database
        """
        return True if self.formula in self.mp_experimental else False
    
    @property
    def is_oxide(self):
        """
        Returns:
            True if the target has oxygen
        """
        return True if "O" in CompTools(self.formula).els else False
    
    @property
    def is_gas(self):
        """
        Returns:
            True if the target is a gas
        """
        return True if self.formula in ['C1O2','H2O1','O1','H1'] else False

    @property
    def in_tm(self):
        """
        Returns:
            True if the formula is in the text-mined dataset (prec or target) (and MP) else False
        """
        if self.formula in self.tm_targets:
            return True
        if self.formula in self.tm_precursors:
            return True
        return False

    @property
    def els(self):
        """
        Returns:
            list of elements in the formula
        """
        return CompTools(self.formula).els

class AnalyzeChemsys():
    def __init__(self, els):
        """

        Returns:
            (additionally)
            tm_precursors (list) : list of precursors in the text-mined dataset that are also in MP
            tm_targets (list) : list of targets in the text-mined dataset that are also in MP
            mp_cmpds (list) : list of compounds in MP
        """
        self.els = els
        # els can easily be extracted from a targeet of interest
        self.tm_precursors = read_json(os.path.join(DATADIR, 'tm_precursors.json'))['data']
        self.mp_data = read_json(os.path.join(DATADIR_cemsbartel, '241002_mp_experimental.json'))['data']
    
    @property
    def flexible_els(self):
        """
        Returns a list of elements
            these elements may not be part of the target chemical system but they could be part of a precursor's chemical system

        Logic:
            if the target is an oxide, we want to consider carbonates and hydroxides as possible precursors
        """
        if 'O' in self.els:
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
            data = self.mp_data
            precursors = list(data.keys())

        # what elements are in the target
        target_els = self.els

        # do we have additional elements to consider
        flexible_els = self.flexible_els

        # how many elements in the target (binary, ternary, etc.)
        nary = len(target_els)
        print(nary)

        # determine allowed chemical systems for precursors
        ## first, just consider (n-1)ary systems that are subsets of the target chemical system
        allowed_els = []
        for n in range(1, nary):
            allowed_els.extend(list(combinations(target_els, n)))
        print(allowed_els)

        # now incorporating our "flexible elements" (basically adding carbonates and hydroxides)
        new_allowed_els = []
        for el in flexible_els:
            for el_combo in allowed_els:
                el_combo = list(el_combo)
                print(el_combo)
                if ("O" in el_combo) and (el not in el_combo) and (len(el_combo) > 1):
                    el_combo.append(el)
                    el_combo = tuple(sorted(el_combo))
                    new_allowed_els.append(el_combo)
        print(new_allowed_els)
        allowed_els = set(allowed_els + new_allowed_els)
        print(allowed_els)
        # filter our big list of precursors down to those that we deemed "possible"
        # precursors = [p for p in precursors if tuple(CompTools(p).els) in allowed_els]
        precs_new = []
        # target_els.extend(flexible_els)
        for p in precursors:
            if tuple(CompTools(p).els) in allowed_els:
                precs_new.append(p)
            elif all([el in target_els for el in CompTools(p).els]):
                precs_new.append(p)
        return list(set(precs_new))


# def check():
    # target = "BaTiO3"
    # at = AnalyzeTarget(target)
    # possible_precursors = at.possible_precursors(restrict_to_tm=True)
    # print("Target = %s" % target)
    # print("Possible precursors = %s" % possible_precursors)
    # return at


# def main():
#     return


# if __name__ == "__main__":
#     main()
