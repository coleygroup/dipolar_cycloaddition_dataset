class DipoleDataPoint:
    '''
    Class storing the settings for a single reaction

    Attributes:
        rxn_id (int): The identifier of the reaction
        xyz_file (str): The name of the xyz-file
        solvent (str): The name of the solvent (optional)
        temp (str): The temperature at which the reaction takes place
        folder_name (str): The name of the folder in which the output of the calculation will be placed
        free_energy (bool): Whether or not thermal corrections need to be computed

    Methods:
        make_folder: Makes the folder in which the reaction profile calculation will be performed
        process_smiles: Removes the atom-mapping from the reaction SMILES (necessary for autodE)
        clean_up_folder: Cleans up the reaction profile folder after calculation has finished (moving intermediate files to archive)
    '''

    def __init__(
        self,
        rxn_id: int = None,
        xyz_file: str = None,
        solvent: str = None,
        temp: float = 298.15,
        folder_name: str = None,
        free_energy: bool = False,
    ):

        self.rxn_id = rxn_id
        self.xyz_file = xyz_file
        self.solvent = solvent
        self.temp = temp
        self.folder_name = folder_name
        self.free_energy = free_energy

    def clean_up_folder(self):
        ''' Cleans up the reaction profile folder after calculation has finished (moving intermediate files to archive) '''

        return (
            f'\tfiles_to_tar = [file for file in os.listdir() if not file.endswith(".csv") and not file.endswith(".xyz")]\n'
            f'\tsubprocess.run(["tar", "-zcf", "intermediate.tar"] + files_to_tar) \n'
            f'\ttime.sleep(2)\n'
            f'\tsubprocess.run(["rm"] + files_to_tar) \n \n'
        )

    def __str__(self):
        ''' Prints output corresponding to the reaction data point '''
        if self.solvent == "gas":
            if self.free_energy :
                return (
                    f'\tdipole_{self.rxn_id} = ade.Molecule("{self.xyz_file.strip(" ")}", name="{self.xyz_file.strip(" ").split(".")[0]}")\n'
                    f'\tdipole_{self.rxn_id}.optimise(method=G16)\n\tdipole_{self.rxn_id}.calc_thermo(temp={self.temp})\n'
                    f'\tdipole_{self.rxn_id}.single_point(method=G16)\n\tdipole_{self.rxn_id}.print_xyz_file()\n\tprint_energies_to_csv(dipole_{self.rxn_id}))\n\n'
                )
            else:
                return (
                    f'\tdipole_{self.rxn_id} = ade.Molecule("{self.xyz_file.strip(" ")}", name="{self.xyz_file.strip(" ").split(".")[0]}")\n'
                    f'\tdipole_{self.rxn_id}.optimise(method=G16)\n\tdipole_{self.rxn_id}.single_point(method=G16)\n'
                    f'\tdipole_{self.rxn_id}.print_xyz_file()\n\tprint_energies_to_csv(dipole_{self.rxn_id})\n\n'
                )

        elif self.solvent != "gas":
            if self.free_energy:
                return (
                    f'\tdipole_{self.rxn_id} = ade.Molecule("{self.xyz_file.strip(" ")}", name="{self.xyz_file.strip(" ").split(".")[0]}",'
                    f'solvent_name="{self.solvent}")\n\tdipole_{self.rxn_id}.optimise(method=G16)\n\tdipole_{self.rxn_id}.calc_thermo(method=G16,temp={self.temp})\n'
                    f'\tdipole_{self.rxn_id}.single_point(method=G16)\n\tdipole_{self.rxn_id}.print_xyz_file()\n\tprint_energies_to_csv(dipole_{self.rxn_id})\n\n'
                )
            else:
                return (
                    f'\tdipole_{self.rxn_id} = ade.Molecule("{self.xyz_file.strip(" ")}", name="{self.xyz_file.strip(" ").split(".")[0]}",'
                    f'solvent_name={self.solvent})\n\tdipole_{self.rxn_id}.optimise(method=G16)\n\tprint(dipole_{self.rxn_id}.__dict__)\n'
                    f'\tdipole_{self.rxn_id}.print_xyz_file()\n\tprint_energies_to_csv(dipole_{self.rxn_id})\n\n'
                )


class InputFile:
    '''
    Class which sets up an autodE inputfile

    Attributes:
        idx (int): The identifier of the inputfile
        ncores (int): The number of cores that will be used during the calculation
        level_of_theory (str): The DFT level of theory (func/basis1/basis2/disp_corr) 
        autode_folder (str): The overarching autodE folder in which the input files need to be stored

    Methods:
        determine_keywords: Determines and formats all the keywords necessary for the autodE calculation
        initialize_input_file: Writes the common lines to the autodE input script
        close_input_file: Closes the input file after all the individual data points have been written
        write_reaction_point: Writes a reaction to the input file
    '''

    def __init__(
        self,
        idx: int = None,
        n_cores: int = 4,
        level_of_theory: str = None,
        autode_folder: str = None,
    ):

        self.idx = idx
        self.n_cores = n_cores
        self.autode_folder = autode_folder

        self.functional = None

        if level_of_theory is not None:
            self.functional = level_of_theory.split("/")[0]
            self.basis_set_low = level_of_theory.split("/")[1]
            self.basis_set_high = level_of_theory.split("/")[2]
            if len(level_of_theory.split("/")) == 4:
                self.disp_corr = level_of_theory.split("/")[3]

        self.file = open(f"{idx}.py", "w")

        self.initialize_input_file()

    def write_dipole_datapoint(self, data_point):
        ''' Writes a reaction to the input file '''
        self.file.write(f"\nos.chdir(os.path.join(path, str({data_point.rxn_id})))\n\n")
        self.file.write("try: \n")
        self.file.write(str(data_point))
        self.file.write(data_point.clean_up_folder())
        self.file.write("except Exception as e: \n")
        self.file.write(f'    print(f"Failure for dipole {data_point.rxn_id}: {{e}}") \n')

    def close_input_file(self):
        ''' Closes the input file after all the individual data points have been written '''
        self.file.close()

    def initialize_input_file(self):
        ''' Writes the common lines to the autodE input script '''
        if self.functional is not None:
            low_opt, grad, opt, opt_ts, hess, sp = self.determine_keywords()

        self.file.write("import autode as ade\n")
        self.file.write("import subprocess\n")
        self.file.write("import time\n")
        self.file.write("import os\n\n")
        self.file.write(f"ade.Config.n_cores = {self.n_cores}\n")
        self.file.write(f"ade.Config.hmethod_conformers = False\n\n")
        self.file.write(f"path = os.path.join(os.getcwd(),'lowest_energy_conformers_RR')\n\n")

        self.file.write("def print_energies_to_csv(mol):\n")
        self.file.write("\twith open('energies.csv', 'w') as f:\n")
        self.file.write("\t\tf.write(f'{mol.name}, {mol.energies.first_potential}, {mol.g_cont}, {mol.h_cont}, {mol.energies.last_potential}')\n\n")

        if self.functional is not None:
            self.file.write("G16 = ade.methods.G16()\n\n")
            self.file.write(
                f"G16.keywords.low_opt = ade.SinglePointKeywords({low_opt})\n"
            )
            self.file.write(
                f"G16.keywords.grad = ade.SinglePointKeywords({grad})\n"
            )
            self.file.write(
                f"G16.keywords.opt = ade.SinglePointKeywords({opt})\n"
            )
            self.file.write(
                f"G16.keywords.opt_ts = ade.SinglePointKeywords({opt_ts})\n"
            )
            self.file.write(
                f"G16.keywords.hess = ade.SinglePointKeywords({hess})\n"
            )
            self.file.write(
                f"G16.keywords.sp = ade.SinglePointKeywords({sp})\n\n"
            )


    def determine_keywords(self):
        ''' Determines and formats all the keywords necessary for the autodE calculation '''
        ts_str = (
            "Opt=(TS, CalcFC, NoEigenTest, MaxCycles=100, MaxStep=10, "
            "NoTrustUpdate, RecalcFC=30)"
        )
        iop_str = "IOp(2/9=2000)"

        try:
            low_opt = [
                self.functional,
                self.basis_set_low,
                "Opt=(loose, maxcycles=10)",
                self.disp_corr,
                iop_str,
            ]
            grad = [
                self.functional,
                self.basis_set_low,
                "Force(NoStep)",
                self.disp_corr,
                iop_str,
            ]
            opt = [self.functional, self.basis_set_low, "Opt", self.disp_corr, iop_str]
            opt_ts = [
                self.functional,
                self.basis_set_low,
                "Freq",
                self.disp_corr,
                ts_str,
                iop_str,
            ]
            hess = [
                self.functional,
                self.basis_set_low,
                "Freq",
                self.disp_corr,
                iop_str,
            ]
            sp = [self.functional, self.basis_set_high, self.disp_corr, iop_str]
        except AttributeError:
            low_opt = [
                self.functional,
                self.basis_set_low,
                "Opt=(loose, maxcycles=10)",
                iop_str,
            ]
            grad = [self.functional, self.basis_set_low, "Force=(NoStep)", iop_str]
            opt = [self.functional, self.basis_set_low, "Opt", iop_str]
            opt_ts = [self.functional, self.basis_set_low, "Freq", ts_str, iop_str]
            hess = [self.functional, self.basis_set_low, "Freq", iop_str]
            sp = [self.functional, self.basis_set_high, iop_str]

        return low_opt, grad, opt, opt_ts, hess, sp


