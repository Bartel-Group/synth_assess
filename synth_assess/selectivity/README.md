# synth_assess.selectivity

This module enables users to generate reaction networks in order to identify target-forming reactions and to determine their associated selectivity. 
In the synth_assess.selectivity.rxn_metrics module, the GammaFromTarget class aggregates the entire workflow for reaction network construction and gamma computation and may be used to identify the most selective reaction to form a target of interest at a specified temperature. For further customization, users may refer to synth_assess.selectivity.entries, which is called in the rxn_metrics module.

For quick generation of target-forming reactions and selectivity determination, please use the following code.
```
# import class to compute Γ 
from synth-assess.selectivity.rxn_metrics import GammaFromTarget
```
To get all reactions to form a given target at a given temperature, run the following
```
all_rxns = GammaFromTarget(target, temperature).get_metrics(gen_data = None, is_gen = None)
```
where the target argument is a formula string and the temperature argument is a float (in Kelvin). If temperature is unspecified, 1073 K is used.
This line returns a list of dictionaries, each detailing a target-forming reaction and the associated Γ for the specified temperature.

To get the optimum reaction to form a given target at a given temperature:
```
opt_rxn = GammaFromTarget(target, temperature).opt_rxn(gen_data = None, is_gen = None)
```
where the target argument is a formula string and the temperature argument is a float (in Kelvin). If temperature is unspecified, 1073 K is used.
This line returns a dictionary, detailing the target-forming reaction with the lowest (most favorable) associated Γ for the specified temperature.

In both cases, if material is not in MP, additional data must be given-- for this purpose use is_gen = True and gen_data as input data. gen_dat must be of the same structure (and contain the same information) as mp_data (refer to synth_assess/data/README for details)
