# dipolar_cycloaddition_dataset
The repository contains the code and auxiliary data associated to the 1,3-dipolar cycloaddition reaction dataset construction project.

## Requirements

1. python 3.10
2. autode 1.2
3. rdkit 2020.03
4. ipykernel 6.9
5. pandas 1.4
6. pebble 4.6
7. xtb 6.3
8. tqdm 4.64
9. pip 2.22
10. rdchiral 1.1

### Conda environment
To set up the conda environment:
```
conda env create --name <env-name> --file environment.yml
```

## Generating the search space and reaction SMILES

The Jupyter notebooks used to generate the search space for both dipoles and dipolarophiles is included in the `dataset_construction` directory. The full list of generated species and the extracted samples can be found in the `full_list_dipoles_dipolarophiles` and `sample_list` sub-directories respectively. The `dataset_construction` directory furthermore contains Python scripts to combine dipoles and dipolarophiles into valid reaction SMILES. The first script generates reaction SMILES for the synthetic and biofragment-based dipolarophiles separately:
```
python construct_dataset_finalized.py
```

By default, the script generates two `.csv` files based on the samples defined in the `sample_list` sub-directory. The resulting files which were used as input for the reaction profile generation workflow are included in the `final_datasets` sub-directory. In the same subdirectory, a combined `.csv` file can also be found in which the two sets of reaction SMILES have been appended.

To generate the azide test reaction SMILES, another Python script in the `dataset_construction` directory needs to be executed:
```
python construct_dataset_azides_finalized.py
```

The outputted `.csv` file used in the automated workflow is also included in the `final_datasets` sub-directory 

## High-throughput reaction profile computation

