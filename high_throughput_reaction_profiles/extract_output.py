from argparse import ArgumentParser
import os
import csv
import pandas as pd
import shutil
import numpy as np

hartree = 627.5094740631

parser = ArgumentParser()
parser.add_argument(
    "--data-file",
    type=str,
    required=True,
    help="input .csv file containing the reaction data",
)
parser.add_argument(
    "--output-folder",
    type=str,
    default="reaction_profiles",
    help="folder containing computed reaction profiles",
)
parser.add_argument(
    "--complexes",
    dest="complexes",
    action="store_true",
    help="whether extract energies of complexes need to be extracted",
)


def get_rxn_smiles_and_targets(filename):
    ''' Read in the input .csv-file, rename rxn_id column and initialize the output columns '''
    df = pd.read_csv(filename)
    df.rename(columns={"Unnamed: 0": "rxn_id"}, inplace=True)

    for output_column in ["E_r", "E_act", "H_r", "H_act", "G_r", "G_act"]:
        df[f"{output_column}"] = df["rxn_smiles"].apply(lambda x: None)

    return df


def extract_energies(row):
    ''' Extract energies from a single line in an energies.csv-file '''
    energy = float(row[-1])
    try:
        enthalpy = float(row[-1]) + float(row[-2])
        free_energy = float(row[-1]) + float(row[-3])
    except Exception:
        return np.array([energy, None, None])

    return np.array([energy, enthalpy, free_energy])


def extract_output(idx, complexes: bool = False):
    ''' Extract and process energy values for each of the species associated with a reaction '''
    r_alt_name = None

    path = os.path.join(os.getcwd(), str(idx))
    file = os.path.join(path, "output/energies.csv")
    if not os.path.isfile(file):
        path = os.path.join(os.getcwd(), f"r{str(idx)}")
        file = os.path.join(path, "output/energies.csv") 

    try:
        with open(file, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            reaction_data = [row for row in csv_reader]
            for row in reaction_data[2:]:
                if row[0].startswith('r0'):
                    r0_energy_array = extract_energies(row)
                elif row[0].startswith('r1'):
                    r1_energy_array = extract_energies(row)
                elif row[0].startswith('p0'):
                    product_energy_array = extract_energies(row)
                elif row[0].startswith('TS'):
                    ts_energy_array = extract_energies(row)
                elif row[0].endswith('_alt'):
                    r_alt_energy_array = extract_energies(row)
                elif row[0].endswith("_reactant"):
                    r_complex_energy_array = extract_energies(row)
                elif row[0].endswith("_product"):
                    p_complex_energy_array = extract_energies(row)

        if r_alt_name is not None:
            if r_alt_name.startswith('r0'):
                reactant_energy_array = r_alt_energy_array + r1_energy_array
            elif r_alt_name.startswith('r1'):
                reactant_energy_array = r_alt_energy_array + r0_energy_array
        else:
            reactant_energy_array = r0_energy_array + r1_energy_array

        reaction_energy_array = (product_energy_array - reactant_energy_array) * hartree

        try:
            activation_energy_array = (ts_energy_array - reactant_energy_array) * hartree
        except UnboundLocalError:
            print(f"No Barrier found for {idx}!")
            activation_energy_array = np.array([None, None, None])

        if complexes:
            complexation_energy_array = (r_complex_energy_array - reactant_energy_array) * hartree

            G_r_complexes = (p_complex_energy_array[-1] - r_complex_energy_array[-1]) * hartree
            try:
                G_act_complexes = (ts_energy_array[-1] - r_complex_energy_array[-1]) * hartree
            except UnboundLocalError:
                G_act_complexes = None

    except Exception:
        print(f"File for {idx} does not exist!")
        if complexes:
            (
                reaction_energy_array, 
                activation_energy_array, 
                complexation_energy_array,
                G_r_complexes,
                G_act_complexes,
            ) = (np.array([None, None, None]), np.array([None, None, None]), np.array([None, None, None]), None, None)
        else:
            reaction_energy_array, activation_energy_array = np.array([None, None, None]), np.array([None, None, None])

    if None not in activation_energy_array:
        print(f"Reaction profile found for {idx}!")

    if complexes:
        return [
            reaction_energy_array,
            activation_energy_array,
            complexation_energy_array,
            G_r_complexes,
            G_act_complexes,
        ]
    else:
        return [reaction_energy_array, activation_energy_array]


def get_xyz(idx, xyz_folder):
    ''' Copy .xyz-files from output folder to xyz-folder; if failed calculation, then just return None '''
    path = os.path.join(os.getcwd(), str(idx))
    folder = os.path.join(path, "output")
    file = os.path.join(path, "output/energies.csv")
    if not os.path.isfile(file):
        path = os.path.join(os.getcwd(), f"r{str(idx)}")
        file = os.path.join(path, "output/energies.csv")
        folder = os.path.join(path, "output")
        if not os.path.isfile(file):
            return None

    xyz_sub_folder = os.path.join(xyz_folder, str(idx))
    if not os.path.isdir(xyz_sub_folder):
        os.mkdir(xyz_sub_folder)
    file_names = os.listdir(folder)
    for file_name in file_names:
        full_file_name = os.path.join(folder, file_name)
        if full_file_name.endswith(".xyz") and "imag_mode" not in full_file_name:
            shutil.copy(full_file_name, xyz_sub_folder)
        elif full_file_name.endswith("energies.csv"):
            shutil.copy(full_file_name, xyz_sub_folder)


if __name__ == "__main__":
    args = parser.parse_args()
    xyz_folder = f'xyz_folder_{args.output_folder.strip("/")}'
    if not os.path.isdir(xyz_folder):
        os.mkdir(xyz_folder)
    dataset = get_rxn_smiles_and_targets(args.data_file)
    original_path = os.getcwd()
    os.chdir(os.path.join(os.getcwd(), args.output_folder))
    if args.complexes:
        dataset["output"] = dataset["rxn_id"].apply(lambda x: extract_output(x, True))
    else:
        dataset["output"] = dataset["rxn_id"].apply(lambda x: extract_output(x))

    dataset["rxn_id"].apply(
        lambda x: get_xyz(x, os.path.join(original_path, xyz_folder))
    )

    dataset["E_r"] = dataset["output"].apply(lambda x: x[0][0])
    dataset["E_act"] = dataset["output"].apply(lambda x: x[1][0])

    dataset["H_r"] = dataset["output"].apply(lambda x: x[0][1])
    dataset["H_act"] = dataset["output"].apply(lambda x: x[1][1])

    dataset["G_r"] = dataset["output"].apply(lambda x: x[0][2])
    dataset["G_act"] = dataset["output"].apply(lambda x: x[1][2])

    if args.complexes:
        dataset["E_complexation"] = dataset["output"].apply(lambda x: x[2][0])
        dataset["H_complexation"] = dataset["output"].apply(lambda x: x[2][1])
        dataset["G_complexation"] = dataset["output"].apply(lambda x: x[2][2])
        dataset["G_r_complexes"] = dataset["output"].apply(lambda x: x[3])
        dataset["G_act_complexes"] = dataset["output"].apply(lambda x: x[4])

    os.chdir(original_path)

    # clean up potential errors
    dataset = dataset[dataset["G_act"].notnull()]
    dataset = dataset[dataset["G_act"] > 0.001]
    dataset = dataset[dataset["G_act"] <= 100]

    if args.complexes:
        if "E_target" in dataset.columns:
            dataset[
                [
                    "rxn_id",
                    "rxn_smiles",
                    "solvent",
                    "temp",
                    "E_target",
                    "E_r",
                    "E_act",
                    "H_r",
                    "H_act",
                    "G_r",
                    "G_act",
                    "E_complexation",
                    "H_complexation",
                    "G_complexation",
                    "G_r_complexes",
                    "G_act_complexes",
                ]
            ].to_csv(f'output_{args.output_folder.strip("/")}.csv')
        else:
            dataset[
                [
                    "rxn_id",
                    "rxn_smiles",
                    "solvent",
                    "temp",
                    "E_r",
                    "E_act",
                    "H_r",
                    "H_act",
                    "G_r",
                    "G_act",
                    "E_complexation",
                    "H_complexation",
                    "G_complexation",
                    "G_r_complexes",
                    "G_act_complexes",
                ]
            ].to_csv(f'output_{args.output_folder.strip("/")}.csv')
    else:
        if "E_target" in dataset.columns:
            dataset[
                [
                    "rxn_id",
                    "rxn_smiles",
                    "solvent",
                    "temp",
                    "E_target",
                    "E_r",
                    "E_act",
                    "H_r",
                    "H_act",
                    "G_r",
                    "G_act",
                ]
            ].to_csv(f'output_{args.output_folder.strip("/")}.csv', index=False)
        else:
            dataset[
                [
                    "rxn_id",
                    "rxn_smiles",
                    "solvent",
                    "temp",
                    "E_r",
                    "E_act",
                    "H_r",
                    "H_act",
                    "G_r",
                    "G_act",
                ]
            ].to_csv(f'output_{args.output_folder.strip("/")}.csv', index=False)
