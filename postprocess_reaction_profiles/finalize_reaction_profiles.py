from argparse import ArgumentParser
import pandas as pd
import os
import shutil
import autode as ade
from autode.mol_graphs import find_cycles
from autode.input_output import xyz_file_to_atoms
from autode.geom import calc_rmsd
import numpy as np
import csv

hartree = 627.5094740631

parser = ArgumentParser()

parser.add_argument(
    "--data-file",
    type=str,
    required=True,
    help="data set file (.csv)",
)
parser.add_argument(
    "--preliminary-profile-folder",
    type=str,
    required=True,
    help="folder containing reaction profiles before postprocessing",
)
parser.add_argument(
    "--output-postprocessing",
    type=str,
    default="output_postprocessing_dipole_conformers.csv",
    help="input .csv file obtained after postprocessing steps",
)
parser.add_argument(
    "--output-folder-postprocessing",
    type=str,
    default="lowest_energy_conformers_RR",
    help="folder containing the optimized dipole .xyz-files and energy values",
)
parser.add_argument(
    "--finalized-profile-folder",
    type=str,
    default="reaction_profiles_finalized",
    help="folder in which to store the finalized reaction profiles",
)


def get_dipole_smiles_dict(csv_file):
    df = pd.read_csv(csv_file)
    df["dipole_reactant"] = df["rxn_smiles"].apply(
        lambda x: [
            smiles
            for smiles in x.split(">")[0].split(".")
            if ("+" in smiles and "-" in smiles)
        ][0]
    )
    dipole_smiles_dict = dict(zip(df.index, df.dipole_reactant))

    return dipole_smiles_dict


def correct_reaction_profiles(
    id, to_run, constrained, preliminary_dir, dipole_output_dir, final_dir
):
    """
    Copy initial output folder to final location and decide whether the alternative dipole conformer needs
    to be included as well (if so, append line to energies.csv and copy .xyz-file)
    """
    shutil.copytree(
        os.path.join(preliminary_dir, str(id)), os.path.join(final_dir, str(id))
    )
    if to_run == True:
        # first get information about the alternative dipole conformer
        try:
            line_alt, energy_alt, name, xyz_alt = get_information_alt_dipole_conf(
                dipole_output_dir, id
            )
            xyz_original = os.path.join(os.path.join(final_dir, str(id)), f"{name}.xyz")
        except:
            print(
                f"An error has taken place during information extraction from alternative dipole conformer for {id}"
            )
            # abort
            return None
        # if constrained and different geometry, paste the line at the end of the original .csv file and copy the xyz-file
        if constrained == True:
            # First do a final check that cyclization has not occured in the alternative dipole conformer;
            # if it did -> remove profile and abort
            cyclization = check_cyclization(xyz_original, xyz_alt)
            if cyclization:
                print(id, cyclization)
                shutil.rmtree(os.path.join(final_dir, str(id)))
                return None
            # if constrained and alternative conformer valid, then you first need to check whether the geometries
            # are really different
            is_unique = check_rmsd(xyz_original, xyz_alt)
            if is_unique:
                include_alt_dipole_conf(final_dir, id, line_alt, xyz_alt)
        else:
            # determine the energy for the original dipole conformer and compare with alternative dipole
            with open(
                os.path.join(os.path.join(final_dir, str(id)), "energies.csv"), "r"
            ) as f:
                line_original = [
                    line for line in f.readlines() if line.startswith(name)
                ][0]
                energy_original = float(line_original.split(",")[-1]) + float(
                    line_original.split(",")[-3]
                )
            # if new energy more than than 1 kJ/0.239 kcal below original and no cyclization in the alternative dipole conformer -> replace
            cyclization = check_cyclization(xyz_original, xyz_alt)
            if (energy_alt - energy_original) * hartree < -0.239 and not cyclization:
                include_alt_dipole_conf(final_dir, id, line_alt, xyz_alt)
            else:
                pass


def get_information_alt_dipole_conf(dipole_output_dir, id):
    """Read in relevant information from the generated RR conformer of the reactant dipole"""
    with open(
        os.path.join(os.path.join(dipole_output_dir, str(id)), "energies.csv"), "r"
    ) as f:
        line_alt = f.readlines()[0]
        energy_alt = float(line_alt.split(",")[-1]) + float(line_alt.split(",")[-3])
        name = line_alt.split(",")[0].rstrip("_alt")
        xyz_alt = os.path.join(
            os.path.join(dipole_output_dir, str(id)), f'{line_alt.split(",")[0]}.xyz'
        )

    return line_alt, energy_alt, name, xyz_alt


def include_alt_dipole_conf(final_dir, id, line_alt, xyz_alt):
    """Append energetic information to energies.csv of full reaction profile and copy .xyz file"""
    with open(os.path.join(os.path.join(final_dir, str(id)), "energies.csv"), "a") as f:
        f.write(line_alt)
    shutil.copy(xyz_alt, os.path.join(final_dir, str(id)))


def get_molecule_no_hydrogens(molecule):
    """Remove all hydrogen atoms from  a list of molecules"""
    molecule.atoms = [atom for atom in molecule.atoms if atom.label != "H"]
    return molecule


def check_rmsd(xyz_original, xyz_alt, threshold_rmsd=0.05):
    """Generate autodE species for the original and alternative dipole conformers and perform an RMSD comparison to verify uniqueness"""
    mol_original = get_molecule_no_hydrogens(
        ade.Species(
            name=xyz_original.rstrip(".xyz"),
            atoms=xyz_file_to_atoms(xyz_original),
            charge=0,
            mult=1,
        )
    )
    mol_alt = get_molecule_no_hydrogens(
        ade.Species(
            name=xyz_alt.rstrip(".xyz"),
            atoms=xyz_file_to_atoms(xyz_alt),
            charge=0,
            mult=1,
        )
    )
    rmsd = calc_rmsd(mol_alt.coordinates, mol_original.coordinates)
    if rmsd < threshold_rmsd:
        return False
    else:
        return True


def check_cyclization(xyz_original, xyz_alt):
    """Check if cyclization has taken place during generation of alternative conformer"""
    ade_mol_original = ade.Molecule(xyz_original, name="original")
    ade_mol_alt = ade.Molecule(xyz_alt, name="alt")

    if len(find_cycles(ade_mol_original.graph)) == len(find_cycles(ade_mol_alt.graph)):
        return False
    else:
        return True


def get_rxn_smiles_and_targets(filename):
    """Read in the input .csv-file, rename rxn_id column and initialize the output columns"""
    df = pd.read_csv(filename)
    df["rxn_id"] = df.index

    for output_column in ["E_r", "E_act", "H_r", "H_act", "G_r", "G_act"]:
        df[f"{output_column}"] = df["rxn_smiles"].apply(lambda x: None)

    return df


def extract_energies(row):
    """Extract energies from a single line in an energies.csv-file"""
    energy = float(row[-1])
    try:
        enthalpy = float(row[-1]) + float(row[-2])
        free_energy = float(row[-1]) + float(row[-3])
    except Exception:
        return np.array([energy, None, None])

    return np.array([energy, enthalpy, free_energy])


def extract_output(main_path, idx):
    """Extract and process energy values for each of the species associated with a reaction"""
    r_alt_name = None

    path = os.path.join(main_path, str(idx))
    file = os.path.join(path, "energies.csv")

    if not os.path.isfile(file):
        path = os.path.join(os.getcwd(), f"r{str(idx)}")
        file = os.path.join(path, "energies.csv")

    try:
        with open(file, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            reaction_data = [row for row in csv_reader]
            for row in reaction_data[2:]:
                if row[0].startswith("r0") and not row[0].endswith("_alt"):
                    r0_energy_array = extract_energies(row)
                elif row[0].startswith("r1") and not row[0].endswith("_alt"):
                    r1_energy_array = extract_energies(row)
                elif row[0].startswith("p0"):
                    product_energy_array = extract_energies(row)
                elif row[0].startswith("TS"):
                    ts_energy_array = extract_energies(row)
                elif row[0].endswith("_alt"):
                    r_alt_energy_array = extract_energies(row)
                    r_alt_name = row[0]

        if r_alt_name is not None:
            if r_alt_name.startswith("r0"):
                reactant_energy_array = r_alt_energy_array + r1_energy_array
            elif r_alt_name.startswith("r1"):
                reactant_energy_array = r_alt_energy_array + r0_energy_array
        else:
            reactant_energy_array = r0_energy_array + r1_energy_array

        reaction_energy_array = (product_energy_array - reactant_energy_array) * hartree

        try:
            activation_energy_array = (
                ts_energy_array - reactant_energy_array
            ) * hartree
        except UnboundLocalError:
            print(f"No Barrier found for {idx}!")
            activation_energy_array = np.array([None, None, None])

    except Exception:
        print(f"File for {idx} does not exist!")
        reaction_energy_array, activation_energy_array = np.array(
            [None, None, None]
        ), np.array([None, None, None])

    return [reaction_energy_array, activation_energy_array]


def clean_up_incorrect_datapoints(dataset, finalized_profile_folder):
    """Remove datapoints involving negative barriers and barriers above 100 kcal/mol"""
    # collect the ids of these problematic datapoints and remove the corresponding folder
    incorrect_ids = (
        dataset[dataset["G_act"] < 0.001]["rxn_id"].tolist()
        + dataset[dataset["G_act"] >= 100]["rxn_id"].tolist()
    )
    for idx in incorrect_ids:
        shutil.rmtree(os.path.join(finalized_profile_folder, str(idx)))
    # remove those ids from the dataset dataframe
    dataset = dataset[dataset["G_act"] > 0.001]
    dataset = dataset[dataset["G_act"] <= 100]

    return dataset


if __name__ == "__main__":
    args = parser.parse_args()
    # initialize
    df_dipoles_alt = pd.read_csv(args.output_postprocessing)
    os.makedirs(args.finalized_profile_folder, exist_ok=True)
    dipole_smiles_dict = get_dipole_smiles_dict(args.data_file)

    # first correct the reaction profiles that need correcting
    df_dipoles_alt.apply(
        lambda x: correct_reaction_profiles(
            x["rxn_id"],
            x["to_run"],
            x["constrained"],
            args.preliminary_profile_folder,
            args.output_folder_postprocessing,
            args.finalized_profile_folder,
        ),
        axis=1,
    )

    # then, extract energy values
    dataset = get_rxn_smiles_and_targets(args.data_file)
    path = os.path.join(os.getcwd(), args.finalized_profile_folder)

    dataset["output"] = dataset["rxn_id"].apply(lambda x: extract_output(path, x))
    dataset["filter"] = dataset["output"].apply(
        lambda x: None not in x[0] and None not in x[1]
    )
    dataset = dataset[dataset["filter"]]

    dataset["E_r"] = dataset["output"].apply(lambda x: x[0][0])
    dataset["E_act"] = dataset["output"].apply(lambda x: x[1][0])

    dataset["H_r"] = dataset["output"].apply(lambda x: x[0][1])
    dataset["H_act"] = dataset["output"].apply(lambda x: x[1][1])

    dataset["G_r"] = dataset["output"].apply(lambda x: x[0][2])
    dataset["G_act"] = dataset["output"].apply(lambda x: x[1][2])

    # clean up potential errors that have slipped through the cracks
    dataset = clean_up_incorrect_datapoints(dataset, path)

    # write output
    dataset.reset_index(inplace=True)
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
    ].to_csv(f'output_{args.finalized_profile_folder.strip("/")}.csv')
