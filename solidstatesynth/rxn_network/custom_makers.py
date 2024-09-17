from solidstatesynth.rxn_network.jobs.core import CalculateCompetitionMaker
from solidstatesynth.rxn_network.jobs.core import GetEntrySetMaker
from solidstatesynth.rxn_network.jobs.schema import CompetitionTaskDocument
from solidstatesynth.rxn_network.jobs.schema import EntrySetDocument
import ray
from solidstatesynth.rxn_network.jobs.utils import get_added_elem_data
from solidstatesynth.rxn_network.reactions.reaction_set import ReactionSet
from solidstatesynth.rxn_network.utils.ray import initialize_ray
from solidstatesynth.rxn_network.entries.entry_set import GibbsEntrySet
from solidstatesynth.rxn_network.core import Composition
from jobflow import SETTINGS, Maker, job
from solidstatesynth.rxn_network.utils.funcs import get_logger

from pymatgen.core.composition import Element

from collections.abc import Iterable
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from copy import deepcopy
from solidstatesynth.rxn_network.entries.utils import initialize_entry

logger = get_logger(__name__)



def process_entries(
    entries: Iterable[ComputedStructureEntry],
    temperature: float,
    e_above_hull: float,
    filter_at_temperature: int | None = None,
    include_nist_data: bool = True,
    include_freed_data: bool = False,
    include_polymorphs: bool = False,
    formulas_to_include: Iterable[str] | None = None,
    calculate_e_above_hulls: bool = False,
    ignore_nist_solids: bool = True,
) -> GibbsEntrySet:
    """Convenience function for processing a set of ComputedStructureEntry objects into a
    GibbsEntrySet with specified parameters. This is used when building entries in most
    of the jobs/flows.

    Args:
        entries: Iterable of ComputedStructureEntry objects. These can be downloaded
            from The Materials Project API or created manually with pymatgen.
        temperature: Temperature [K] for determining Gibbs Free Energy of
            formation, dGf(T).
        e_above_hull: Energy above hull (eV/atom) for thermodynamic stability threshold;
            i.e., include all entries with energies below this value.
        filter_at_temperature: Temperature (in Kelvin) at which entries are filtered for
            thermodynamic stability (e.g., room temperature). Generally, this often
            differs from the synthesis temperature.
        include_nist_data: Whether to include NIST-JANAF data in the entry set.
            Defaults to True.
        include_freed_data: Whether to include FREED data in the entry set. Defaults
            to False. WARNING: This dataset has not been thoroughly tested. Use at
            your own risk!
        include_polymorphs: Whether to include non-ground state polymorphs in the entry
            set. Defaults to False.
        formulas_to_include: An iterable of compositional formulas to ensure are
            included in the processed dataset. Sometimes, entries are filtered out that
            one would like to include, or entries don't exist for those compositions.
        calculate_e_above_hulls: Whether to calculate e_above_hull and store as an
            attribute in the data dictionary for each entry.
        ignore_nist_solids: Whether to ignore NIST data for solids with high melting
            points (Tm >= 1500 ºC). Defaults to True.

    Returns:
        A GibbsEntrySet object containing entry objects with the user-specified
        constraints.
    """
    temp = temperature
    if filter_at_temperature:
        temp = filter_at_temperature

    entry_set = GibbsEntrySet.from_computed_entries(
        entries=entries,
        temperature=temp,
        include_nist_data=include_nist_data,
        include_freed_data=include_freed_data,
        ignore_nist_solids=ignore_nist_solids,
        #############
        apply_atmospheric_co2_correction=False
        #############
    )
    included_entries = [initialize_entry(f, entry_set) for f in formulas_to_include] if formulas_to_include else []

    entry_set = entry_set.filter_by_stability(e_above_hull=e_above_hull, include_polymorphs=include_polymorphs)
    entry_set.update(included_entries)  # make sure these aren't filtered out

    if filter_at_temperature and (filter_at_temperature != temperature):
        entry_set = entry_set.get_entries_with_new_temperature(temperature)

    if calculate_e_above_hulls:
        entry_set = GibbsEntrySet(deepcopy(entry_set), calculate_e_above_hulls=True)

    return entry_set



class CustomESM(GetEntrySetMaker):
    
    @job(entries="entries", output_schema=EntrySetDocument)
    def make(self, chemsys: str):
        """Returns a job that acquires a GibbsEntrySet for the desired chemical system.

        NOTE: This job stores the entry set in an additional store called
        "entries". This needs to be configured through a user's jobflow.yaml file. See
        "additional_stores".

        Args:
            chemsys: The chemical system of the entry set to be acquired.
        """
        entry_db = SETTINGS.JOB_STORE.additional_stores.get(self.entry_db_name)  # pylint: disable=no-member

        if entry_db:
            logger.info(f"Using user-specified Entry DB: {self.entry_db_name}")
            property_data = self.property_data
            if property_data is None:
                property_data = ["theoretical"]
            elif "theoretical" not in property_data:
                property_data.append("theoretical")

            #Hope I don't need this
            entries = None
            # entries = get_all_entries_in_chemsys_from_entry_db(
            #     entry_db,
            #     chemsys,
            # )
        else:
            try:
                from mp_api.client import MPRester
            except ImportError as err:
                raise ImportError("You may need to install the Materials Project API: pip install -U mp-api") from err

            kwargs = {}
            if self.MP_API_KEY:
                kwargs["api_key"] = self.MP_API_KEY

            elems = {Element(i) for i in chemsys.split("-")}

            if len(elems) <= 25:
                with MPRester(**kwargs) as mpr:
                    entries = mpr.get_entries_in_chemsys(
                        elements=chemsys,
                        additional_criteria={"thermo_types": ["GGA_GGA+U"]},
                    )
            else:  # this approach is faster for big systems
                other_elems = self._get_exclude_elems(elems)

                with MPRester(**kwargs) as mpr:
                    docs = mpr.summary.search(exclude_elements=other_elems, all_fields=False, deprecated=False)
                mpids = [d.material_id for d in docs]

                with MPRester(**kwargs) as mpr:
                    entries = mpr.get_entries(mpids, additional_criteria={"thermo_types": ["GGA_GGA+U"]})

        if self.custom_entries:
            entries.extend(self.custom_entries)

        entries = process_entries(
            entries,
            temperature=self.temperature,
            e_above_hull=self.e_above_hull,
            include_nist_data=self.include_nist_data,
            include_freed_data=self.include_freed_data,
            filter_at_temperature=self.filter_at_temperature,
            include_polymorphs=self.include_polymorphs,
            formulas_to_include=self.formulas_to_include,
            calculate_e_above_hulls=self.calculate_e_above_hulls,
            ignore_nist_solids=self.ignore_nist_solids,
        )

        doc = EntrySetDocument(
            entries=entries,
            e_above_hull=self.e_above_hull,
            include_polymorphs=self.include_polymorphs,
            formulas_to_include=self.formulas_to_include,
        )
        doc.task_label = self.name

        return doc




class CustomCCM(CalculateCompetitionMaker):

    @job(rxns="rxns", output_schema=CompetitionTaskDocument)
    def make(self, rxn_sets: list[ReactionSet], entries: GibbsEntrySet, target_formula: str):
        """Returns a job that calculates competition scores and/or chemical potential
        distances for all synthesis reactions to a target phase given a provided list of
        reaction sets.

        NOTE: This job stores the reaction set in an additional store called
        "rxns". This needs to be configured through a user's jobflow.yaml file. See
        "additional_stores".

        Args:
            rxn_sets: a list of reaction sets making up all enumerated reactions in the
                chemical reaction network of interest. These will automatically be
                combined and reprocessed to match the specified conditions (open_elem +
                chempot).
            entries: The entry set used to enumerate all provided reactions. This will
                be used to facilitate selectivity calculations and ensure all reaction
                sets can be easily combined.
            target_formula: The formula of the desired target phase. This will be used
                to identify all synthesis reactions (i.e., those that produce the
                target).

        Returns:
            A job that returns synthesis reactions to the target phase, decorated with
            the relevant selectivity metrics.
        """
        if not ray.is_initialized():
            initialize_ray()

        target_formula = Composition(target_formula).reduced_formula
        added_elements, added_chemsys = get_added_elem_data(entries, [target_formula])

        logger.info("Loading reactions..")
        all_rxns = rxn_sets[0]
        for rxn_set in rxn_sets[1:]:
            all_rxns = all_rxns.add_rxn_set(rxn_set)

        if self.open_elem:  # reinitialize with open element
            all_rxns = ReactionSet(
                all_rxns.entries,
                all_rxns.indices,
                all_rxns.coeffs,
                self.open_elem,
                self.chempot,
                all_rxns.all_data,
            )

        size = len(all_rxns)  # need to get size before storing in ray

        logger.info("Identifying target reactions...")

        target_rxns = all_rxns.get_rxns_by_product(target_formula, return_set=True)

        logger.info(f"Identified {len(target_rxns)} target reactions out of {size} total reactions.")

        logger.info("Removing unnecessary reactions from total reactions to save memory...")

        all_target_reactants = {reactant.reduced_formula for r in target_rxns for reactant in r.reactants}

        ########### My changes
        if self.open_elem:
            all_target_reactants.add("O2")
        ###########


        all_rxns = all_rxns.get_rxns_by_reactants(
            all_target_reactants,
            return_set=True,  # type: ignore
        )

        logger.info(f"Keeping {len(all_rxns)} out of {size} total reactions...")

        size = len(all_rxns)

        logger.info("Placing reactions in ray object store...")

        all_rxns = ray.put(all_rxns.as_dict())
        logger.info("Beginning competition calculations...")

        decorated_rxns = target_rxns

        if self.calculate_competition:
            decorated_rxns = self._get_competition_decorated_rxns(target_rxns, all_rxns, size)

        logger.info("Calculating chemical potential distances...")

        if self.calculate_chempot_distances:
            decorated_rxns = self._get_chempot_decorated_rxns(decorated_rxns, entries)

        logger.info("Saving decorated reactions.")
        results = ReactionSet.from_rxns(decorated_rxns, entries=entries)

        data = {
            "rxns": results,
            "target_formula": target_formula,
            "open_elem": self.open_elem,
            "chempot": self.chempot,
            "added_elements": added_elements,
            "added_chemsys": added_chemsys,
            "calculate_competition": self.calculate_competition,
            "calculate_chempot_distances": self.calculate_chempot_distances,
            "cpd_kwargs": self.cpd_kwargs,
        }

        doc = CompetitionTaskDocument(**data)
        doc.task_label = self.name
        return doc