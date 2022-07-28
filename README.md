# dipolar_cycloaddition_dataset
This repository contains the code and auxiliary data associated to the 1,3-dipolar cycloaddition reaction dataset construction project. Code is provided "as-is". Minor edits may be required to tailor the scripts for different computational systems.

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

Additionally, in order to execute the autodE high-throughput reaction profile computation workflow, Gaussian09/Gaussian16 needs to be accessible. More information about the autodE module can be found [here](https://github.com/duartegroup/autodE).

### Conda environment
To set up the conda environment:
```
conda env create --name <env-name> --file environment.yml
```

## Generating the search space and reaction SMILES

The Jupyter notebooks used to generate the search space for both dipoles and dipolarophiles are included in the `dataset_construction` directory. The full list of generated species and the extracted samples used throughout the workflow can be found in the `full_list_dipoles_dipolarophiles` and `sample_list` sub-directories respectively. The `dataset_construction` directory furthermore contains Python scripts to combine dipoles and dipolarophiles into valid reaction SMILES. The first script generates reaction SMILES for the synthetic and biofragment-based dipolarophiles separately:
```
python construct_dataset_finalized.py
```

By default, the script generates two `.csv` files based on the samples defined in the `sample_list` sub-directory. Note that the reaction SMILES are ordered based on the number of electrons present in the reacting system, which facilitates the construction of batches of autodE reaction profile computations with a similar computational cost (see below). The resulting files which were used as input for the reaction profile generation workflow are included in the `final_datasets` sub-directory. In the same subdirectory, a combined `.csv` file can also be found in which the two sets of reaction SMILES have been appended.

To generate the azide test reaction SMILES, another Python script in the `dataset_construction` directory needs to be executed:
```
python construct_dataset_azides_finalized.py
```

The outputted `.csv` file used in the automated workflow is also included in the `final_datasets` sub-directory 

## High-throughput reaction profile computation

Input files for high-throughput reaction profile computations, based on the reaction SMILES outputted in the previous section, can be generated with the help of the `initialize.py` script in the `high_throughput_reaction_profiles` directory as follows:
```
python initialize.py --data-file <path to the input .csv file> --num-input-files <total number of files to be generated> --autodE-folder <name of the autodE folder to be generated> [--n-cores <number of cores per computation>] [--DFT-theory <functional/low_basis_set/high_basis_set/dispersion_correction>] [--free-energy] [--complexes]
```
The script batches up the reaction SMILES in individual Python files, each containing a linear sequence of autodE reaction profile computations, and places them in the `autodE_folder` (the number of files can be defined in the `num_input_files` command line option). Additionally, the script generates sub-directories for each individual reaction SMILES in the `autodE_folder` (the sub-directories are numbered based on the indices in the input file). By default, 4 cores are used per computation, no free energy corrections nor complexes are computed and the PBE0-D3(BJ)/def2tzvp//PBE0-D3(BJ)/def2svp level of theory is used. These default options can be adjusted with the help of the corresponding command line flags.

Once the input files have been generated, they can be executed with the help of a job array (potentially [with triples](https://supercloud.mit.edu/job-arrays-llsub-triples-3-steps)).

To resume an interrupted workflow, a re-initialization can be performed which will generate new input files (with an 'r' prefix) containing only the autodE  computations which haven't been run yet:
```
python re_initialize.py --data-file <path to the input .csv file> --num_input-files <total number of files to be generated> --autodE-folder <name of the autodE folder to be generated> [--n-cores <number of cores per computation>] [--DFT-theory <functional/low_basis_set/high_basis_set/dispersion_correction>] [--free-energy] [--complexes]
```

Finally, the `high_throughput_reaction_profiles` directory also contains a script to extract the relevant output from the `autodE_folder`. This script can be executed as follows:
```
python extract_output.py --data-file <path to the input .csv file> --output-folder <autodE_folder> [--complexes]
```

This script will copy all final .xyz files as well as the `energies.csv` to a new directory (`xyz_folder_<output_folder>`). Additionally, it creates a .csv file containing the successfully computed reaction SMILES together with the activition energies/enthalpies/free energies as well as the reaction energies/enthalpies/free energies (`output_<output_folder>.csv`). 

Note that many reactant dipoles require postprocessing to make their geometry compatible with the selected TS conformer. To facilitate this postprocessing, the scripts below have also been included in this repository.

## Post-processing reaction SMILES to ensure stereo-compatibility of the dipoles

The main script in the `postprocess_reaction_profiles` directory, `postprocess_reaction_profiles.py`, checks for compatibility between the reactant dipole conformer and the conformation taken on by this species in the TS. It takes as an input the full data set `.xyz` file, as well as the output folder obtained after execution of the `extract_output.py` script from the previous section. The output is a new folder, `lowest_energy_conformers_RR`, containing directories and xyz-files for each individual dipole which might need to be replaced (vide infra), as well as a `.csv` file indicating for each dipole whether or not it needs to be run and whether the geometry was constrained. 

More specifically, the script extracts the dipole geometry from the TS `.xyz` file, optimizes the geometry, checks rotatability of the individual dihedral angles over the (partially) double bonds connecting the 3 reactive centers in the dipole and then performs a (constrained) randomize and relax conformer search. If constraints were applied during the search, i.e., if at least one of the bonds was not rotatable, then the lowest energy conformer is always saved. If this is not the case, then the newly generated conformer is only saved when it is lower in energy than the originally registered conformer.   

The script can be executed as follows:
```
python postprocess_reaction_profiles.py --data-file <path to the input .csv file> --xyz-folder <xyz_folder> [--n-cores <number of cores to be used>]
```

Once a `lowest_energy_conformers_RR` directory has been generated, a script making heavy use of autodE can be called to geometrically and energetically refine the alternative RR conformers generated during postprocessing. This script is very similar to the `initialize.py` script and can be executed as follows:
```
python initialize_dipole_calc.py --data_file <path to .csv file with the postprocessing output> --num_input_files <total number of files to be generated> [--autodE_folder <name of the autodE folder>] [--n_cores <number of cores per computation>] [--DFT_theory <functional/low_basis_set/high_basis_set/dispersion_correction>] [--free_energy]
```

The script batches up the dipoles in individual Python files, each containing a linear sequence of autodE geometry optimizations and frequency calculations, and places them in the selected `autodE_folder` containing the alternative dipole conformer `.xyz` files (the default name is once again `lowest_energy_conformers_RR`).  Again, the number of files can be defined in the `num_input_files` command line option. By default, 4 cores are used per computation, no free energy corrections nor complexes are computed and the PBE0-D3(BJ)/def2tzvp//PBE0-D3(BJ)/def2svp level of theory is used. These default options can be adjusted with the help of the corresponding command line flags.

Finally, the `finalize_reaction_profiles.py` script can be called to combine the preliminary reaction profiles with the corrected dipole conformers and recompute the activation and reaction energies. The input of this script consists of the full data set `.xyz` file, the output `.csv` file from `postprocess_reaction_profiles.py` (vide supra), the preliminary, i.e., uncorrected, reaction profile directory and the `lowest_energy_conformers_RR` directory which is obtained after execution of the two preceding scripts in this section. As output, the script generates a finalized reaction profile folder as well as a `.csv` file with all the reaction SMILES and their computed activation and reaction energies. It can be executed as follows:
```
python finalize_reaction_profiles.py --data-file <data set .csv file> --preliminary-profile-folder <folder containing reaction profiles before postprocessing> [--output-postprocessing <input .csv file obtained after postprocessing steps>] [--output-folder-postprocessing <folder containing the optimized dipole .xyz-files and energy values>] [--finalized-profile-folder <folder in which to store the finalized reaction profiles>]
```

## Final dataset

The `.csv` file, containing the curated final dataset, can be found in the main directory (`output_reaction_profiles_finalized.csv`). The geometries and energies for each of the reactions can be found at xxx.
## References

If (parts of) this workflow are used as part of a publication please cite the associated paper:

xxx

Additionally, since the workflow makes heavy use of autodE, please also cite the paper in which this code was originally presented:
```
@article{autodE,
  doi = {10.1002/anie.202011941},
  url = {https://doi.org/10.1002/anie.202011941},
  year = {2021},
  publisher = {Wiley},
  volume = {60},
  number = {8},
  pages = {4266--4274},
  author = {Tom A. Young and Joseph J. Silcock and Alistair J. Sterling and Fernanda Duarte},
  title = {{autodE}: Automated Calculation of Reaction Energy Profiles -- Application to Organic and Organometallic Reactions},
  journal = {Angewandte Chemie International Edition}
}
```
