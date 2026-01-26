This file describes how to run the synthesizability prediction for each model tested. The models's github pages are self-explanatory, but the necessary steps for only pretrained predictions are outlined below. Each model's section in this file will first describe how to arrange the downloaded repos to test on the given structures, and then there is documentation from the repo on how to run the predictor once all the data is in place.







# Synthesizability-stoi-CGNF

#### Place id_prop.csv in the folder and run the command in [2]


## Prerequisites
Python3<br> Numpy<br> Pytorch<br> Pymatgen<br>


## Usage
### [1] Define a customized data format and prepare atomic embedding vector file for generation of CGNF
To input crystal structures to Synthesizability-stoi-CGNF, you will need to define a customized dataset and pre-generate CGNF as pickle files for bootstrap aggregating in semi-supervised learning. Note that this is required for both training and predicting.
Following files should be needed to generate CGNF.
#### id_prop.csv: a CSV file with two columns for positive data(synthesizable) and unlabeled data(not-yet-synthesized). The first column recodes a inorganic composition (The formula string format of Composition class in Pymatgen package is recommended), and the second column recodes the value (1 = positive, 0 = unlabeled) according to whether they were synthesized already or not.
#### cgcnn_hd_rcut4_nn8.element_embedding.json: a JSON file containing atomic embedding vectors for generation of CGNF


### [2] Predict synthesizability of new crystals with pre-trained models
`python predict_PU_learning.py --bag 100 --data id_prop.csv --embedding cgcnn_hd_rcut4_nn8.element_embedding.json --modeldir ./models`<br>

Load composition information from 'id_prop_test.csv' file for test materials and pre-trained models from 'models' folder.<br>
Predict synthesizability of crystal composition in id_prop_test.csv file using the loaded models.<br>
Result of bootstrap aggregating is saved as 'test_results_ensemble_100models.csv'




# SynthNN

#### Replace input_formulas.txt in the repo with the provided input_formulas.txt. Open the jupyter notebook SynthNN_predict.ipynb and run the first cell for imports, and then run the cell labeled "Option 2: Read in formulas from text file." The results are shown in the file output_formula_preds.txt. Additinoal information is given below.


## Prerequisites
Requirements:
- Python
- [Pymatgen](https://pymatgen.org/installation.html)
- [Tensorflow](https://www.tensorflow.org/install)


### Predict Synthesizability
Predicting the synthesizability of a material composition with a pre-trained version of SynthNN can be done with SynthNN_predict.ipynb.
We recommend referring to the below performance metrics when choosing a decision threshold to label a material as synthesizable or not. The below table indicates the performance of
SynthNN of a dataset with a 20:1 ratio of unsynthesized:synthesized examples. Note, a threshold value of '0.10' means that any material with a SynthNN output greater than 0.10 is taken to be synthesizable, which leads to low precision but high recall.
Threshold | Precision | Recall | 
| :---: | :---: | :---: |
0.10 | 0.239 | 0.859 |
0.20 | 0.337 | 0.783 |
0.30 | 0.419 | 0.721 |
0.40 | 0.491 | 0.658 |
0.50 | 0.563 | 0.604 |
0.60 | 0.628 | 0.545 |
0.70 | 0.702 | 0.483 |
0.80 | 0.765 | 0.404 |
0.90 | 0.851 | 0.294 |




# Synthesizability-PU-CGCNN

#### To generate the crystal graphs in [1], move the provided folders cif_files and saved_crystal_graph into the model directory. Copy the cifs from cifs_relaxed (unzip first) into cif_files and run the command in [1]. After generating the crystal graphs, run the prediction in [2].


## Prerequisites
Python3<br> Numpy<br> Pytorch<br> Pymatgen<br>


## Usage
### [1] Define a customized dataset and generate crystal graphs
To input crystal structures to Synthesizability-PU-CGCNN, you will need to define a customized dataset and pre-generate crystal graph as pickle files for bootstrap aggregating in partially supervised learning. Note that this is required for both training and predicting.
If you want to use cif data in the folder named as “cif_files”, following files should be needed to generate crystal graph.
#### 1) id_prop.csv: a CSV file with two columns for positive data(synthesizable) and unlabeled data(not-yet-synthesized). The first column recodes a unique ID for each crystal, and the second column recodes the value (1 = positive, 0 = unlabeled) according to whether they were synthesized already or not.
#### 2) atom_init.json: a JSON file that stores the initialization vector for each element.
#### 3) ID.cif: a CIF file that recodes the crystal structure, where ID is the unique ID for the crystal.
ex) If you want to generate crystal graph with cutoff radius 8A, maximum 12 neighbors:<br>
`python generate_crystal_graph.py --cifs ./cif_files --n 12 --r 8 --f ./saved_crystal_graph`<br>
Then, you will obtain preloaded crystal graph files in folder “saved_crystal_graph”<br>


### [2] Predict synthesizability of new crystals with pre-trained models
`python predict_PU_learning.py --bag 100 --graph ./saved_crystal_graph --cifs ./cif_files --modeldir ./trained_models`<br>

Load crystal graph information from 'saved_crystal_graph' folder and pre-trained models from 'trained_models' folder.<br>
Predict synthesizability of crystal structures in 'cif_files' folder (with id_prop.csv file) using the loaded models.<br>
Result of bootstrap aggregating is saved as 'test_results_ensemble_100models.csv'




# TSDNN

#### Copy the root_dir folder into the /data directory and copy all the cifs from cifs_relaxed into root_dir. Run the prediction command using the pretrained model.


##  Prerequisites

This package requires:

- [PyTorch](http://pytorch.org)
- [scikit-learn](http://scikit-learn.org/stable/)
- [pymatgen](http://pymatgen.org)

If you are new to Python, the easiest way of installing the prerequisites is via [conda](https://conda.io/docs/index.html). After installing [conda](http://conda.pydata.org/), run the following command to create a new [environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) named `cgcnn` and install all prerequisites:

```bash
conda upgrade conda
conda create -n tsdnn python=3 scikit-learn pytorch torchvision pymatgen -c pytorch -c conda-forge
```

Alternatively, you can import our conda environment from the `environment.yml file:

```bash
conda create -n tsdnn
conda install -f environment.yml
```

*Note: this code is tested for PyTorch v1.0.0+ and is not compatible with versions below v0.4.0 due to some breaking changes.

This creates a conda environment for running TSDNN. Before using TSDNN, activate the environment by:

```bash
conda activate tsdnn
```


### Define a customized dataset 

The structure of the `root_dir` should be:
```
root_dir
├── data_test.csv
├── atom_init.json
├── id0.cif
├── id1.cif
├── ...
```

### Predict material properties with a pre-trained TSDNN model

Before predicting the material properties, you will need to:

- [Define a customized dataset](#define-a-customized-dataset) at `root_dir` for all the crystal structures that you want to predict.
- Obtain a pre-trained TSDNN model (example found in checkpoints/pre-trained/pre-train.pth.tar).

Then, in directory `synth-tsdnn`, you can predict the properties of the crystals in `root_dir`:

```bash
python python predict.py checkpoints/pre-trained/synthesizability.pth.tar data/root_dir
```

After predicting, you will get one file in `synth-tsdnn` directory:

- `predictions.csv`: stores the `ID` and predicted value for each crystal in test set.




# SynCoTrain

#### Follow the steps by copying the crystal data file into the schnet_pred/data/ folder and running the command below.


## Installation
It is recommended to create a virtual environment with mamba and miniforge to install the different packages easily. Start by installing mamba according to the instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

Start by cloning this repository in your preferred path as shown below:
```bash
cd /path/to/parent/directory
git clone git@github.com:BAMeScience/SynCoTrainMP.git
```
Next, navigate to the cloned directory. You can create the appropriate mamba environment there based on the `sync.yml` file:
```bash
cd SynCoTrain
mamba env create -f condaEnvs/sync.yml
mamba activate sync
```
This might take a while, as all the required packages are being installed. Please note that you may need to change ther exact version of dgl and cudatoolkit based on your current setup. You can check your current cuda version using the `nvidia-smi` command. Then, you can search for a compatible dgl with cuda using the command `mamba search dgl --channel conda-forge`. Pick a version of dgl earlier than 2.0.0 which is compatible with your cudatoolkit.

Once the packages are installed, you may activate the `sync` conda environment and install this repository with the following commands:
```bash
pip install -e .
```


## Predicting Synthesizability of Oxides

If you are only interested in predicting the synthesizability of oxides, there’s no need to train the model from scratch. The current version of SynCoTrain comes pre-trained for synthesizability prediction of oxide crystals, using SchNet as the classifier.

### How to Predict Synthesizability for Your Data

1. **Prepare Your Data**: Save your crystal data as a pickled DataFrame and place it in the `schnet_pred/data` directory. For example:
```
schnet_pred/data/chemeleon_filtered_structures_3.pkl
```
2. **Run Prediction**: Use the following command to feed your data into the model:

```bash
python schnet_pred/predict_schnet.py --input_file chemeleon_filtered_structures_3
```
3. **View Results**: The prediction results will be saved in the following location:
```bash
schnet_pred/results/chemeleon_filtered_structures_3_predictions.csv
```
The result will be saved in `schnet_pred/results/chemeleon_filtered_structures_3_predictions.csv`.
