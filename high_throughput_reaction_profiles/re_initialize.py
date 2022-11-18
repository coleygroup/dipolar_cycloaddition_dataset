from argparse import ArgumentParser
import os
import pandas as pd
import math

from lib import ReactionDataPoint, InputFile

parser = ArgumentParser()

parser.add_argument(
    "--data-file",
    type=str,
    required=True,
    help="input .csv file containing the reaction data",
)
parser.add_argument(
    "--num-input-files",
    type=int,
    default=288,
    help="number of calculations to run in parallel",
)
# autodE calculation
parser.add_argument(
    "--autodE-folder",
    type=str,
    default="reaction_profiles",
    help="folder for DFT calculation",
)
parser.add_argument(
    "--DFT-theory",
    type=str,
    default=None,
    help="level of theory for the DFT calculation"
    "functional/basis_set_low/basis_set_high/disp_corr_keyword",
)
parser.add_argument(
    "--n-cores", type=int, default=8, help="number of cores for autodE calculations"
)
parser.add_argument(
    "--free-energy",
    dest="free_energy",
    action="store_true",
    help="To compute free energies",
)
parser.add_argument(
    "--complexes",
    dest="complexes",
    action="store_true",
    help="whether or not to compute complexes (not compatible with --free_energy,"
    "but freq computation is done by default during sp-step)",
)


def not_yet_run(idx, failed_calc_list):
    """Determines if a data point has already been run through index matching"""
    if idx in failed_calc_list:
        return False
    else:
        return True


if __name__ == "__main__":
    args = parser.parse_args()

    output = f"output_{args.autodE_folder.strip('/')}.csv"

    failed_calculations = [
        archive.split(".tar.gz")[0]
        for archive in os.listdir(args.autodE_folder)
        if archive.endswith(".tar.gz")
    ]

    df_input = pd.read_csv(args.data_file)
    df_input.reset_index(inplace=True)

    df_output = pd.read_csv(output)

    output_list = df_output["rxn_id"].tolist()

    df_input["not_yet_run"] = df_input["index"].apply(
        lambda x: str(x) not in failed_calculations and x not in output_list
    )

    df_to_do = df_input[df_input["not_yet_run"]]

    name = os.path.splitext(args.data_file)[0]

    df_to_do.rename(columns={"index": "rxn_id"}, inplace=True)

    df_to_do["data_point"] = df_to_do.apply(
        lambda x: ReactionDataPoint(
            f"r{str(x['rxn_id'])}",
            x["rxn_smiles"],
            x["solvent"],
            x["temp"],
            os.path.join(args.autodE_folder, f"r{str(x['rxn_id'])}"),
            args.free_energy,
            args.complexes,
        ),
        axis=1,
    )

    data_point_list = df_to_do["data_point"].tolist()

    os.chdir(args.autodE_folder)

    num_calc = math.ceil(len(data_point_list) / args.num_input_files)

    for i in range(0, args.num_input_files):
        file = InputFile(i + 1, args.n_cores, args.DFT_theory, args.autodE_folder, relaunch=True)
        for j in range(0, num_calc):
            try:
                file.write_reaction_point(
                    data_point_list[j * args.num_input_files + i],
                    j * args.num_input_files + i,
                )
            except IndexError:
                break
