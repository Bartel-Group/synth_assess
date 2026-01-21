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
This module enables the user to generate reactions to form a target of interest and to identify the most thermodynamically selective reactions. This module enables the user to reconstruct the entire pipeline, including customization of data used in reaction generation and constraints on reaction generation (refer to solidstatesynth.selectivity.rxn_networks), but if you are simply seeking to generate selective reactions associated with a particular target using our settings, you can run the following code:

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

If this module is used, please cite the following:

```
[1] McDermott, M. J.; Dwaraknath, S. S.; Persson, K. A. A Graph-Based Network for Predicting Chemical Reaction Pathways in Solid-State Materials Synthesis. Nat. Commun. 2021, 12 (1), 3097. https://doi.org/10.1038/s41467-021-23339-x.

[2]	McDermott, M. J.; McBride, B. C.; Regier, C. E.; Tran, G. T.; Chen, Y.; Corrao, A. A.; Gallant, M. C.; Kamm, G. E.; Bartel, C. J.; Chapman, K. W.; Khalifah, P. G.; Ceder, G.; Neilson, J. R.; Persson, K. A. Assessing Thermodynamic Selectivity of Solid-State Reactions for the Predictive Synthesis of Inorganic Materials. ACS Cent. Sci. 2023, 9 (10), 1957–1975. https://doi.org/10.1021/acscentsci.3c01051.
```

## Material generation from input chemical spaces (solidstatesynth.gen)

This module enables users to generate new materials in a specified chemical space using [Chemeleon](https://github.com/hspark1212/chemeleon/) and to compute material energetics using CHGNET.

If this module is used for material generation, please cite the following:
```
[3]	Park, H.; Onwuli, A.; Walsh, A. Exploration of Crystal Chemical Space Using Text-Guided Generative Artificial Intelligence. Nat. Commun. 2025, 16 (1), 4379. https://doi.org/10.1038/s41467-025-59636-y.![image](https://github.umn.edu/user-attachments/assets/7f47dcd2-a3d3-478b-8f32-18b4d24f0ea9)

```
If this module is used for energy computation, please cite the following:
```
[4] Deng, B.; Zhong, P.; Jun, K.; Riebesell, J.; Han, K.; Bartel, C. J.; Ceder, G. CHGNet as a Pretrained Universal Neural Network Potential for Charge-Informed Atomistic Modelling. Nat. Mach. Intell. 2023, 5 (9), 1031–1041. https://doi.org/10.1038/s42256-023-00716-3.![image](https://github.umn.edu/user-attachments/assets/f5a5ae3a-6456-49ef-acdb-9976fa2d668b)

```

## Material synthesizability prediction (solistatesynth.pred)
This module enables users to apply five synthesizability predictors to a material of interest. Two of the five models (PU-CGNF [5] and SynthNN [6]) take only formula as input while the other three (PU-CGCNN [7], SynCoTrain [8], and TSDNN [9]) require structural inputs as well.

If this module is used, please cite the following (depending which predictive models are used)
```
[5] Jang, J.; Gu, G. H.; Noh, J.; Kim, J.; Jung, Y. Structure-Based Synthesizability Prediction of Crystals Using Partially Supervised Learning. J. Am. Chem. Soc. 2020, 142 (44), 18836–18843. https://doi.org/10.1021/jacs.0c07384.

[6]	Jang, J.; Noh, J.; Zhou, L.; Gu, G. H.; Gregoire, J. M.; Jung, Y. Synthesizability of Materials Stoichiometry Using Semi-Supervised Learning. Matter 2024, 7 (6), 2294–2312. https://doi.org/10.1016/j.matt.2024.05.002.

[7]	Antoniuk, E. R.; Cheon, G.; Wang, G.; Bernstein, D.; Cai, W.; Reed, E. J. Predicting the Synthesizability of Crystalline Inorganic Materials from the Data of Known Material Compositions. Npj Comput. Mater. 2023, 9 (1), 155. https://doi.org/10.1038/s41524-023-01114-4.

[8]	Amariamir, S.; George, J.; Benner, P. SynCoTrain: A Dual Classifier PU-Learning Framework for Synthesizability Prediction. Digit. Discov. 2025, 4 (6), 1437–1448. https://doi.org/10.1039/D4DD00394B.

[9]	Gleaves, D.; Fu, N.; Dilanga Siriwardane, E. M.; Zhao, Y.; Hu, J. Materials Synthesizability and Stability Prediction Using a Semi-Supervised Teacher-Student Dual Neural Network. Digit. Discov. 2023, 2 (2), 377–391. https://doi.org/10.1039/D2DD00098A.
```

This work acknowledges contributions from [pydmclab](https://github.com/Bartel-Group/pydmclab/blob/main/README.md).
