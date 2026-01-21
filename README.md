# solidstatesynth

A package for assessing a candidate target for solid-state synthesis. This package was developed in conjunction with the following work:

# Installation

```
# Clone the repository
git clone https://github.com/bartel-group/solidstatesynth.git

# Change to the repository directory
cd solidstatesynth

# Install the package
pip install .
```

# Modules

## Reaction generation and selectivity assessment (solidstatesynth.selectivity)
This module enables the user to generate reactions to form a target of interest and to identify the most thermodynamically selective reactions. This module enables the user the reconstruct the entire pipeline, including customization of data used in reaction generation and constraints on reaction generation (refer to solidstatesynth.selectivity.rxn_networks), but if you are simply seeking to generate selective reactions associated with a particular target using our settings, you can run the following code:

```
from solidstatesynth.selectivity.rxn_metrics import GammaFromTarget
```
To get all reactions for a given target and temperature (if temperature is unspecified, 1073 K is used):

```
all_rxns = GammaFromTarget(target, temperature).get_metrics(gen_data = None, is_gen = None)
```

To get the optimum rxn for a given target and temperature (if temperature is unspecified, 1073 K is used):
```
opt_rxn = GammaFromTarget(target, temperature).opt_rxn(gen_data = None, is_gen = None)
```
In both cases, if material is not in MP, additional data must be given-- for this purpose use is_gen = True and gen_data as input data (refer to /solidstatesynth/data/README.md for further details).

If this module is used, please consider citing the following:

```
[1] McDermott, M. J.; Dwaraknath, S. S.; Persson, K. A. A Graph-Based Network for Predicting Chemical Reaction Pathways in Solid-State Materials Synthesis. Nat. Commun. 2021, 12 (1), 3097. https://doi.org/10.1038/s41467-021-23339-x.

[2]	McDermott, M. J.; McBride, B. C.; Regier, C. E.; Tran, G. T.; Chen, Y.; Corrao, A. A.; Gallant, M. C.; Kamm, G. E.; Bartel, C. J.; Chapman, K. W.; Khalifah, P. G.; Ceder, G.; Neilson, J. R.; Persson, K. A. Assessing Thermodynamic Selectivity of Solid-State Reactions for the Predictive Synthesis of Inorganic Materials. ACS Cent. Sci. 2023, 9 (10), 1957–1975. https://doi.org/10.1021/acscentsci.3c01051.
```

- Material generation from input chemical spaces solidstatesynth.gen

serves to generate possible reactions to form a target of interest and to compute the selectivity associated with these reactions (solidstatesynth.selectivity). We also offer the mechanism to generate new materials from input chemical spaces (solidstatesynth.gen) and to compute material energetics using CHGNET and to predict synthesizability using different predictive models  (solidstatesynth.predict). Figures from the reference paper can also be regenerated using this package (solidstatesynth.plotting)

# Installation

```
# Clone the repository
git clone https://github.com/bartel-group/solidstatesynth.git

# Change to the repository directory
cd solidstatesynth

# Install the package
pip install .
```

This package enables the user the reconstruct the entire pipeline, including customization of data used in reaction generation and constraints on reaction generation (refer to solidstatesynth.selectivity.entries and solidstatesynth.selectivity.rxn_networks, but if you are simply seeking to generate selective reactions associated with a particular target using our settings, you can run the following code (with solidstatesynth installed).

```
from solidstatesynth.selectivity.rxn_metrics import GammaFromTarget
```
To get all reactions for a given target and temperature (if temperature is unspecified, 1073 K is used):

```
all_rxns = GammaFromTarget(target, temperature).get_metrics(gen_data = None, is_gen = None)
```

To get the optimum rxn for a given target and temperature (if temperature is unspecified, 1073 K is used):
```
opt_rxn = GammaFromTarget(target, temperature).opt_rxn(gen_data = None, is_gen = None)
```
In both cases, if material is not in MP, additional data must be given-- for this purpose use is_gen = True and gen_data as input data (refer to /solidstatesynth/data/README.md for further details).

If this package is used for selectivity assessment, please consider citing the following works:

```
[1] McDermott, M. J.; Dwaraknath, S. S.; Persson, K. A. A Graph-Based Network for Predicting Chemical Reaction Pathways in Solid-State Materials Synthesis. Nat. Commun. 2021, 12 (1), 3097. https://doi.org/10.1038/s41467-021-23339-x.

[2]	McDermott, M. J.; McBride, B. C.; Regier, C. E.; Tran, G. T.; Chen, Y.; Corrao, A. A.; Gallant, M. C.; Kamm, G. E.; Bartel, C. J.; Chapman, K. W.; Khalifah, P. G.; Ceder, G.; Neilson, J. R.; Persson, K. A. Assessing Thermodynamic Selectivity of Solid-State Reactions for the Predictive Synthesis of Inorganic Materials. ACS Cent. Sci. 2023, 9 (10), 1957–1975. https://doi.org/10.1021/acscentsci.3c01051.
```

If predictive models are used, please refer to the following works:

```
[3] Jang, J.; Gu, G. H.; Noh, J.; Kim, J.; Jung, Y. Structure-Based Synthesizability Prediction of Crystals Using Partially Supervised Learning. J. Am. Chem. Soc. 2020, 142 (44), 18836–18843. https://doi.org/10.1021/jacs.0c07384.

[4]	Jang, J.; Noh, J.; Zhou, L.; Gu, G. H.; Gregoire, J. M.; Jung, Y. Synthesizability of Materials Stoichiometry Using Semi-Supervised Learning. Matter 2024, 7 (6), 2294–2312. https://doi.org/10.1016/j.matt.2024.05.002.

[5]	Antoniuk, E. R.; Cheon, G.; Wang, G.; Bernstein, D.; Cai, W.; Reed, E. J. Predicting the Synthesizability of Crystalline Inorganic Materials from the Data of Known Material Compositions. Npj Comput. Mater. 2023, 9 (1), 155. https://doi.org/10.1038/s41524-023-01114-4.

[6]	Amariamir, S.; George, J.; Benner, P. SynCoTrain: A Dual Classifier PU-Learning Framework for Synthesizability Prediction. Digit. Discov. 2025, 4 (6), 1437–1448. https://doi.org/10.1039/D4DD00394B.

[7]	Gleaves, D.; Fu, N.; Dilanga Siriwardane, E. M.; Zhao, Y.; Hu, J. Materials Synthesizability and Stability Prediction Using a Semi-Supervised Teacher-Student Dual Neural Network. Digit. Discov. 2023, 2 (2), 377–391. https://doi.org/10.1039/D2DD00098A.
```

In addition, this work acknowledges contributions from [pydmclab](https://github.com/Bartel-Group/pydmclab/blob/main/README.md).
