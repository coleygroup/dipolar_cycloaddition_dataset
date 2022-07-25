from argparse import ArgumentParser
import os
import pandas as pd
import math

from lib import DipoleDataPoint, InputFile


parser = ArgumentParser()

parser.add_argument(
    "--data_file",
    type=str,
    required=False,
    help="input .csv file containing the reaction data",
)
parser.add_argument(
    "--num_input_files", type=int, default=288, help="number of calculations to run in parallel"
)

# autodE calculation
parser.add_argument(
    "--autodE_folder",
    type=str,
    default="lowest_energy_conformers_RR",
    help="folder for the DFT calculations",
)
parser.add_argument(
    "--DFT_theory",
    type=str,
    default=None,
    help="level of theory for the DFT calculation"
    "functional/basis_set_low/basis_set_high/disp_corr_keyword",
)
parser.add_argument(
    "--n_cores", type=int, default=4, help="number of cores for autodE calculations"
)
parser.add_argument(
    "--free_energy",
    dest="free_energy",
    action="store_true",
    help="To compute free energies",
)


if __name__ == "__main__":
    args = parser.parse_args()

    name = os.path.splitext(args.data_file)[0]

    if not os.path.isdir(args.autodE_folder):
        os.mkdir(args.autodE_folder)

    df = pd.read_csv(args.data_file)

    df = df[df["to_run"] != False]

    df["data_point"] = df.apply(
        lambda x: DipoleDataPoint(
            str(x["rxn_id"]),
            x["name"],
            "water",
            "298.15",
            os.path.join(args.autodE_folder, str(x["rxn_id"])),
            args.free_energy,
        ),
        axis=1,
    )

    data_point_list = df["data_point"].tolist()

    os.chdir(args.autodE_folder)

    num_calc = math.ceil(len(data_point_list) / args.num_input_files)

    for i in range(0, args.num_input_files):
        file = InputFile(i, args.n_cores, args.DFT_theory, args.autodE_folder)
        for j in range(0, num_calc):
            try:
                file.write_dipole_datapoint(data_point_list[j * args.num_input_files + i])
            except IndexError:
                break

        file.close_input_file()
