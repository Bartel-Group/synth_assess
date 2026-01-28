# PU prediction models
This folder includes all the files and procedures to replicate the PU model predictions on the structure test set generated with Chemeleon. Each model will have to be downloaded from the respective GitHub page, but the generated materials have been formatted to fit the specifications of each model. Additionally, the models have their own descriptive README files for retraining or testing different datasets, but the following instructions will suffice for replicating these predictions. For any use of a PU model, cite using the included reference.

# PU-CGNF
## Model download
`git clone https://github.com/kaist-amsg/Synthesizability-stoi-CGNF.git`

## Prerequisites
* python
* numpy
* pytorch
* pymatgen

## Usage
1. Copy **id_prop.csv** from **synth_assess/pred/PU-CGNF/** into **Synthesizability-stoi-CGNF/**
2. Run the prediction in **Synthesizability-stoi-CGNF/** with `python predict_PU_learning.py --bag 100 --data id_prop.csv --embedding cgcnn_hd_rcut4_nn8.element_embedding.json --modeldir ./models`
3. Results are in the file **test_results_ensemble_100models.csv**

## Reference
* Jang, J.; Noh, J.; Zhou, L.; Gu, G. H.; Gregoire, J. M.; Jung, Y. Synthesizability of Materials Stoichiometry Using Semi-Supervised Learning. Matter 2024, 7 (6), 2294–2312. https://doi.org/10.1016/j.matt.2024.05.002.

# SynthNN
## Model download
`git clone https://github.com/antoniuk1/SynthNN.git`

## Prerequisites
* python
* pymatgen
* tensorflow

## Usage
1. Copy **input_formulas.txt** from **synth_assess/pred/SynthNN/** into **SynthNN/**
2. Open the jupyter notebook **SynthNN_predict.ipynb** and run the first cell for imports
3. Run the third cell labeled "Option 2: Read in formulas from text file."
4. Results are in the file **output_formula_preds.txt**

## Reference
* Antoniuk, E. R.; Cheon, G.; Wang, G.; Bernstein, D.; Cai, W.; Reed, E. J. Predicting the Synthesizability of Crystalline Inorganic Materials from the Data of Known Material Compositions. Npj Comput. Mater. 2023, 9 (1), 155. https://doi.org/10.1038/s41524-023-01114-4.

# PU-CGCNN
## Model download
`git clone https://github.com/kaist-amsg/Synthesizability-PU-CGCNN.git`

## Prerequisites
* python
* numpy
* pytorch
* pymatgen

## Usage
1. Copy the folder **synth_assess/pred/PU-CGCNN/cif_files/** with the files **atom_init.json** and **id_prop.csv** into **Synthesizability-PU-CGCNN/**
2. Make a new directory **Synthesizability-PU-CGCNN/saved_crystal_graph/** - The folder must exist even if it is empty to generate the crystal graphs
3. Unzip the structure files in **synth_assess/data/data/cifs_relaxed.zip** and copy all the CIFs directly into **Synthesizability-PU-CGCNN/cif_files/** so that the directory has **atom_init.json**, **id_prop.csv**, and ~22,000 CIFs
4. Generate the crystal graphs by running `python generate_crystal_graph.py --cifs ./cif_files --n 12 --r 8 --f ./saved_crystal_graph` in **Synthesizability-PU-CGCNN/**
5. Run the prediction model in **Synthesizability-PU-CGCNN/** with `python predict_PU_learning.py --bag 100 --graph ./saved_crystal_graph --cifs ./cif_files --modeldir ./trained_models`
6. Results are in the file **test_results_ensemble_100models.csv**

## Reference
* Jang, J.; Gu, G. H.; Noh, J.; Kim, J.; Jung, Y. Structure-Based Synthesizability Prediction of Crystals Using Partially Supervised Learning. J. Am. Chem. Soc. 2020, 142 (44), 18836–18843. https://doi.org/10.1021/jacs.0c07384.

# SynCoTrain
## Model download
`git clone https://github.com/BAMeScience/SynCoTrainMP.git`

## Prerequisites
Create a new environment with **SynCoTrainMP/condaEnvs/sync.yml** and install the repository

## Usage
1. Copy the file **chemeleon_filtered_structures_3.pkl** into **SynCoTrainMP/schnet_pred/data/**
2. Run the prediction model in **SynCoTrainMP/** with `python schnet_pred/predict_schnet.py --input_file chemeleon_filtered_structures_3`
3. Results are in the file that looks like **SynCoTrainMP/schnet_pred/results/chemeleon_filtered_structures_3_ ... .csv** 

## Reference
* Amariamir, S.; George, J.; Benner, P. SynCoTrain: A Dual Classifier PU-Learning Framework for Synthesizability Prediction. Digit. Discov. 2025, 4 (6), 1437–1448. https://doi.org/10.1039/D4DD00394B.

# TSDNN
## Model download
`git clone https://github.com/usccolumbia/tsdnn.git`

## Prerequisites
* python
* scikit-learn
* pytorch
* torchvision
* pymatgen
* conda-forge

## Usage
1. Copy the folder **synth_assess/pred/tsdn/root_dir/** with the files **atom_init.json** and **data_test.csv** into **tsdnn/data/**
2. Unzip the structure files in **synth_assess/data/data/cifs_relaxed.zip** and copy all the CIFs directly into **tsdnn/data/root_dir** so that the directory has **atom_init.json**, **data_test.csv**, and ~22,000 CIFs
3. Run the prediction model in **tsdnn/** with `python python predict.py checkpoints/pre-trained/synthesizability.pth.tar data/root_dir`
4. Results are in the file **tsdnn/results/predictions/predictions_0.csv**

## Reference
* Gleaves, D.; Fu, N.; Dilanga Siriwardane, E. M.; Zhao, Y.; Hu, J. Materials Synthesizability and Stability Prediction Using a Semi-Supervised Teacher-Student Dual Neural Network. Digit. Discov. 2023, 2 (2), 377–391. https://doi.org/10.1039/D2DD00098A.

