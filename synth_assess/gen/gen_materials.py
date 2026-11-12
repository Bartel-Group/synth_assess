from chemeleon.modules.chemeleon import Chemeleon
import os
from pymatgen.core import Structure
from mp_api.client import MPRester
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pydmclab.mlp.chgnet.dynamics import CHGNetRelaxer
from pydmclab.core.comp import CompTools
from pydmclab.core.struc import StrucTools

from pydmclab.utils.handy import read_json, write_json
import itertools
from smact.screening import smact_validity
from tqdm import tqdm
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from collections import defaultdict
from torch import Tensor
import random


class GenerateMaterials:
    def __init__(self, systems, exclude_formulas=None, n_samples=100, max_stoich=8, max_natoms=40, max_factor=8, max_comps=None, save_dir="./generated_materials", chgnet_filter=None):
        """
        Initializes the GenerateMaterials class.
        Args:
            systems (list): List of chemical systems as tuples of element strings to generate materials for. Ex: [("Li", "Mn", "O"), ("Y", "Ba", "Cu", "O")].
            exclude_formulas (list, optional): List of formulas to exclude from generation, i.e. known formulas in Materias Project. If None, will query them instead.
            n_samples (int, optional): Number of samples to generate for each composition (only will use lowest energy sample).
            max_stoich (int, optional): Maximum stoichiometry for the compositions.
            max_natoms (int, optional): Maximum number of atoms in the generated structures.
            max_factor (int, optional): Maximum factor for the composition generation.
            max_comps (int, optional): After all the possible compositions for the system are enumerated, randomly chooses this many to generate structures for. If None, all compositions are used.
            save_dir (str, optional): Directory to save the generated materials.
            chgnet_filter (float, optional): If float (eV/atom), does not save any structures above the filter. If None, no filter is applied.
        """

        self.systems = systems
        self.n_samples = n_samples
        self.max_stoich = max_stoich
        self.max_natoms = max_natoms
        self.max_factor = max_factor
        self.chgnet_filter = chgnet_filter
        self.max_comps = max_comps

        # If the excluded formulas are not specified, query the Materials Project for them
        if exclude_formulas is None:
            system_strings = ["-".join(sorted(system)) for system in systems]
            print("Querying Materials Project for excluded formulas...")
            with MPRester() as mpr:
                entries = mpr.materials.summary.search(chemsys=system_strings, all_fields=False, fields=["formula_pretty"])
            self.exclude_formulas = list(set([CompTools(i.formula_pretty).pretty for i in entries]))
        else:
            self.exclude_formulas = exclude_formulas

        self.save_dir = save_dir
        self.data = {}

    # Adapted from Chemeleon's navigate_chemical_system to filter out known material compositions
    def generate_materials(self):
        """
        Generates materials for the specified chemical systems using the Chemeleon model.
        The generated materials are saved in CIF format in the specified directory.
        """
        # Initialize Chemeleon model
        chemmy = Chemeleon.load_composition_model()
        exclude = self.exclude_formulas

        random.seed(611)

        # Loop througb each chemical system
        for sys_index, chemsys in enumerate(self.systems):
            # Create directory for the chemical system
            chemsys_str = "-".join(sorted(chemsys))
            path = os.path.join(self.save_dir, chemsys_str)
            os.makedirs(path, exist_ok=True)

            print(f"Generating materials for {chemsys} ({sys_index+1}/{len(self.systems)})")

            # Enumerates all the possible compositions for materials in the chemical system, within the specifications
            all_compositions = [
                Composition(
                    {el: amt for el, amt in zip(chemsys, amt_list)}
                ).reduced_composition
                for amt_list in itertools.product(range(1, self.max_stoich + 1), repeat=len(chemsys))
                if max(amt_list) > 0
            ]

            # Filter out compositions that are invalid to SMACT or are excluded 
            valid_compositions = [comp for comp in all_compositions if smact_validity(comp)]
            unique_valid_compositions = [i for i in set(valid_compositions) if CompTools(i).pretty not in exclude]

            # Filter out compositions that are too large for the max_natoms
            unique_valid_compositions = [i for i in unique_valid_compositions if i.num_atoms <= self.max_natoms]

            # Randomly choose a subset of the valid compositions to generate structures for
            if len(unique_valid_compositions) > self.max_comps:
                print(f"Selecting {self.max_comps} compositions from {len(unique_valid_compositions)} valid options")
                unique_valid_compositions = random.sample(unique_valid_compositions, self.max_comps)


            # Samples all the structures for the valid compositions
            sm = StructureMatcher()
            collections_gen_st = []


            for comp_index, comp in enumerate(tqdm(unique_valid_compositions)):
                
                print(f"\tSampling for {comp} ({comp_index+1}/{len(unique_valid_compositions)})")
                reduced_natoms = int(comp.num_atoms)
                comp = comp.reduced_composition.alphabetical_formula

                valid_gen_st_list = []
                for f in range(1, self.max_factor + 1):
                    if reduced_natoms * f > self.max_natoms:
                        break
                    n_atoms = reduced_natoms * f
                    text_input = comp
                    print(
                        f"\t\tSampling {self.n_samples} structures for {text_input} with {n_atoms} atoms"
                    )

                    # sampling
                    gen_atoms_list = chemmy.sample(
                        text_input=text_input,
                        n_atoms=n_atoms,
                        n_samples=self.n_samples,
                    )
                    if gen_atoms_list is None:
                        continue
                    gen_st_list = [
                        AseAtomsAdaptor.get_structure(atoms) for atoms in gen_atoms_list
                    ]

                    # validity check
                    for st in gen_st_list:
                        if max(st.lattice.abc) > 60:
                            continue
                        if st.composition.reduced_composition not in unique_valid_compositions:
                            continue
                        valid_gen_st_list.append(st)

                # unique structures
                unique_gen_st_list = [out[0] for out in sm.group_structures(valid_gen_st_list)]
                print(f"\tNumber of unique structures: {len(unique_gen_st_list)}")
                collections_gen_st.extend(unique_gen_st_list)

            # final unique structures
            collections_gen_st = [out[0] for out in sm.group_structures(collections_gen_st)]
            print(f"Number of final unique structures: {len(collections_gen_st)}")

            # Fix this
            # if self.chgnet_filter != None:
            #     energies = self.calculate_hull_energies(chemsys, collections_gen_st)
            #     collections_gen_st = [st for st in collections_gen_st if energies[st] < self.chgnet_filter]



            # save cif files
            idx_list = defaultdict(int)
            for st in collections_gen_st:
                # get composition
                comp = st.composition.reduced_composition.alphabetical_formula.replace(" ", "")
                idx_list[comp] += 1
                # save atoms to cif
                atoms = AseAtomsAdaptor.get_atoms(st)
                filename = f"gen_{comp}_{len(atoms)}_{idx_list[comp]}.cif"
                atoms.write(os.path.join(path, filename))
            print(f"Results saved in {path}")
                

    def calculate_hull_energies(self, chemsys, structs, verbose=False, save_relaxed_cifs=False):
        """
        Calculates the decompositions energies of all listed cifs in the given chemical system.
        First creates a pymatgen phase diagram for the chemical system, then applies the CHGNet model to each structure for the energies.
        Args:
            chemsys (tuple): Chemical system as a tuple of element strings. Ex: ("Li", "Mn", "O").
            structs (dict[str, Structure]): Dictionary of structures to calculate the decomposition energies for. The keys are the names of the structures and the values are the pymatgen Structure objects.
        Returns:
            energies (dict[str, float]): Dictionary of decomposition energies for each structure. The keys are the names of the structures and the values are the decomposition energies.
        """
        energies = {}

        if save_relaxed_cifs:
            dir = os.path.join(self.save_dir, "-".join(sorted(chemsys))+"_chgnet_relaxed")
            os.makedirs(dir, exist_ok=True)

        # Load the CHGNet model
        chgnet = CHGNetRelaxer()

        # Generate the phase diagram for the chemical system
        with MPRester() as mpr:
            entries = mpr.get_entries_in_chemsys(chemsys, additional_criteria={"thermo_types": ["GGA_GGA+U"]})

        pd = PhaseDiagram(entries=entries, elements=[Element(el) for el in chemsys])

        # Calculate the energy, and decomposition energy for each generated material
        for index, name in enumerate(structs):

            
            energy = None
            if os.path.exists(os.path.join(dir, name+"_chgnet_relaxed.cif")):
                st_relaxed = StrucTools(os.path.join(dir, name+"_chgnet_relaxed.cif")).structure
                energy = float(chgnet.predict_structure(st_relaxed)["e"])
            
            else:

                st = structs[name]
                
                chgnet_results = chgnet.relax(st, verbose=False)
                st_relaxed = chgnet_results["final_structure"]
                energy = chgnet_results["final_energy"] / st_relaxed.num_sites

            print(f"\t({index+1}/{len(structs)}) {st_relaxed.formula}") if verbose else None

            new_entry = PDEntry(st_relaxed.composition, st_relaxed.composition.num_atoms*energy, name=st_relaxed.formula)
            
            hull_energy = pd.get_phase_separation_energy(new_entry)

            energies[name] = hull_energy

            if save_relaxed_cifs and not os.path.exists(os.path.join(dir, name+"_chgnet_relaxed.cif")):
                StrucTools(st_relaxed).structure_to_cif(name+"_chgnet_relaxed", dir)

        return energies
    
    def get_struct_stabilities(self, verbose=False, save_relaxed_cifs=False):
        """
        Returns the stabilities of the generated structures for the specified systems.
        The stabilities (decomposition energies) are calculated using the CHGNet model (for energy) and the phase diagram.
        Returns:
            stability_dict (dict): Dictionary containing the stabilities (decomposition energies) of the generated structures.
        """
        stability_dict = {}

        # Loop through each chemical system and calculate stabilities
        for system_list in self.systems:
            system = "-".join(sorted(system_list))

            system_path = os.path.join(self.save_dir, system)

            print("Calculating stabilities for system:", system) if verbose else None
            
            structs = {}
            # Calculate the stability for each generated material
            for cif in os.listdir(system_path):

                name = cif.split(".")[0]
                structs[name] = StrucTools(os.path.join(system_path, cif)).structure
                

            energies = self.calculate_hull_energies(system.split("-"), structs, verbose=verbose, save_relaxed_cifs=save_relaxed_cifs)
            stability_dict[system] = energies

        return stability_dict
