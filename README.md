# synth-assess

A package for assessing a candidate target for solid-state synthesis. This tool allows users to generate and assess the selectivity of reactions to form a target of interest, and generate new targets, to reproduce our results, or for new targets. 

# Installation

Navigate to the desired repository destination, and then run the following:

```
# Clone the repository
git clone https://github.com/bartel-group/synth-assess.git

# Change to the repository directory
cd synth_assess

# Install the package
pip install .
```

# Modules

## Reaction generation and selectivity assessment (synth-assess.selectivity)
This module enables the user to generate reactions to form a target of interest and to identify the most thermodynamically selective reactions. Users may reconstruct the entire pipeline, including customization of data used in reaction generation and constraints on reaction generation (refer to solidstatesynth.selectivity.rxn_networks and solidstatesynth.selectivity.entries), but if you are simply seeking to generate reactions associated with a particular target using our settings, and to determine reaction selectivity, you can run the following code:

```
# import class to compute Γ 
from synth_assess.selectivity.rxn_metrics import GammaFromTarget
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

If this module is used, please cite [1] and [2]. Refer to examples for further details.

## Material generation from input chemical spaces (synth-assess.gen)

This module enables users to generate new materials in a specified chemical space using [Chemeleon](https://github.com/hspark1212/chemeleon/) and to compute material energetics using CHGNET.

If this module is used for material generation, please cite [3]. If this module is used for energy computation, please cite [4].

## Material synthesizability prediction (synth-assess.pred)
This module enables users to replicate our results for five synthesizability predictors applied to generative models. Two of the five models (PU-CGNF [5] and SynthNN [6]) take only formula as input while the other three (PU-CGCNN [7], SynCoTrain [8], and TSDNN [9]) require structural inputs as well. Note that we do not import these packages, but offer users the input files needed to make predictions on the materials used in this assessment.

For further information, please refer to [5-9].


```
[1] McDermott, M. J.; Dwaraknath, S. S.; Persson, K. A. A Graph-Based Network for Predicting Chemical Reaction Pathways in Solid-State Materials Synthesis. Nat. Commun. 2021, 12 (1), 3097. https://doi.org/10.1038/s41467-021-23339-x.

[2]	McDermott, M. J.; McBride, B. C.; Regier, C. E.; Tran, G. T.; Chen, Y.; Corrao, A. A.; Gallant, M. C.; Kamm, G. E.; Bartel, C. J.; Chapman, K. W.; Khalifah, P. G.; Ceder, G.; Neilson, J. R.; Persson, K. A. Assessing Thermodynamic Selectivity of Solid-State Reactions for the Predictive Synthesis of Inorganic Materials. ACS Cent. Sci. 2023, 9 (10), 1957–1975. https://doi.org/10.1021/acscentsci.3c01051.

[3]	Park, H.; Onwuli, A.; Walsh, A. Exploration of Crystal Chemical Space Using Text-Guided Generative Artificial Intelligence. Nat. Commun. 2025, 16 (1), 4379. https://doi.org/10.1038/s41467-025-59636-y.

[4] Deng, B.; Zhong, P.; Jun, K.; Riebesell, J.; Han, K.; Bartel, C. J.; Ceder, G. CHGNet as a Pretrained Universal Neural Network Potential for Charge-Informed Atomistic Modelling. Nat. Mach. Intell. 2023, 5 (9), 1031–1041. https://doi.org/10.1038/s42256-023-00716-3.

[5]	Jang, J.; Noh, J.; Zhou, L.; Gu, G. H.; Gregoire, J. M.; Jung, Y. Synthesizability of Materials Stoichiometry Using Semi-Supervised Learning. Matter 2024, 7 (6), 2294–2312. https://doi.org/10.1016/j.matt.2024.05.002.

[6]	Antoniuk, E. R.; Cheon, G.; Wang, G.; Bernstein, D.; Cai, W.; Reed, E. J. Predicting the Synthesizability of Crystalline Inorganic Materials from the Data of Known Material Compositions. Npj Comput. Mater. 2023, 9 (1), 155. https://doi.org/10.1038/s41524-023-01114-4.

[7] Jang, J.; Gu, G. H.; Noh, J.; Kim, J.; Jung, Y. Structure-Based Synthesizability Prediction of Crystals Using Partially Supervised Learning. J. Am. Chem. Soc. 2020, 142 (44), 18836–18843. https://doi.org/10.1021/jacs.0c07384.

[8]	Amariamir, S.; George, J.; Benner, P. SynCoTrain: A Dual Classifier PU-Learning Framework for Synthesizability Prediction. Digit. Discov. 2025, 4 (6), 1437–1448. https://doi.org/10.1039/D4DD00394B.

[9]	Gleaves, D.; Fu, N.; Dilanga Siriwardane, E. M.; Zhao, Y.; Hu, J. Materials Synthesizability and Stability Prediction Using a Semi-Supervised Teacher-Student Dual Neural Network. Digit. Discov. 2023, 2 (2), 377–391. https://doi.org/10.1039/D2DD00098A.
```

# Citation
If using this package, please consider citing the accompanying paper:

```
[1] Schlesinger, J.; Hjaltason, S.; Szymanski, N.; Bartel, C. Thermodynamic assessment of machine learning models for solid-state synthesis prediction. ** In preparation **. 2026 
  
```
