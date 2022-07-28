import os
from rdkit import Chem


class ReactionDataPoint:
    """
    Class storing the settings for a single reaction

    Attributes:
        rxn_id (int): The identifier of the reaction
        rxn_smiles (str): The reaction SMILES string associated with the reaction
        solvent (str): The name of the solvent (optional)
        temp (str): The temperature at which the reaction takes place
        folder_name (str): The name of the folder in which the output of the calculation will be placed
        free_energy (bool): Whether or not thermal corrections need to be computed
        complexes (bool): Whether or not to include complexes in the computation

    Methods:
        make_folder: Makes the folder in which the reaction profile calculation will be performed
        process_smiles: Removes the atom-mapping from the reaction SMILES (necessary for autodE)
        clean_up_folder: Cleans up the reaction profile folder after calculation has finished (moving intermediate files to archive)
    """

    def __init__(
        self,
        rxn_id: int = None,
        rxn_smiles: str = None,
        solvent: str = None,
        temp: float = 298.15,
        folder_name: str = None,
        free_energy: bool = False,
        complexes: bool = False,
    ):

        self.rxn_id = rxn_id
        self.rxn_smiles = self.process_smiles(rxn_smiles)
        self.solvent = solvent
        self.temp = temp
        self.folder_name = folder_name
        self.free_energy = free_energy
        self.complexes = complexes

        self.make_folder()

    def make_folder(self):
        """Makes the folder in which the reaction profile calculation will be performed"""
        if not os.path.isdir(self.folder_name):
            os.mkdir(self.folder_name)

    def process_smiles(self, rxn_smiles):
        """Removes the atom-mapping from the reaction SMILES (necessary for autodE)"""
        if ":" in rxn_smiles:
            r_smiles, p_smiles = rxn_smiles.split(">>")
            r_mol = Chem.MolFromSmiles(r_smiles)
            p_mol = Chem.MolFromSmiles(p_smiles)
            for atom in r_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in p_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            return f"{Chem.MolToSmiles(r_mol)}>>{Chem.MolToSmiles(p_mol)}"
        else:
            return rxn_smiles

    def clean_up_folder(self):
        """Cleans up the reaction profile folder after calculation has finished (moving intermediate files to archive)"""
        path = os.path.join(os.getcwd(), self.rxn_id)
        if self.free_energy and self.complexes:
            dir_list = [
                "conformers",
                "reactants_and_products",
                "single_points",
                "thermal",
                "complexes",
                "transition_states",
            ]
        elif self.free_energy:
            dir_list = [
                "conformers",
                "reactants_and_products",
                "single_points",
                "thermal",
                "transition_states",
            ]
        elif self.complexes:
            dir_list = [
                "conformers",
                "reactants_and_products",
                "single_points",
                "complexes",
                "transition_states",
            ]
        else:
            dir_list = [
                "conformers",
                "reactants_and_products",
                "single_points",
                "transition_states",
            ]

        dir_list = list(map(lambda x: os.path.join(path, x), dir_list))

        return (
            f'    subprocess.run(["tar", "-zcf", "{os.path.join(path, "intermediate.tar")}"] + {dir_list}) \n'
            f'    subprocess.run(["rm", "-r"] + {dir_list}) \n \n'
        )

    def __str__(self):
        """Prints output corresponding to the reaction data point"""
        if self.solvent == "gas":
            if self.free_energy and self.complexes:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f"temp={self.temp}) \n \n    rxn_{self.rxn_id}.calculate_reaction_profile(with_complexes=True, free_energy=True) \n \n"
                )
            elif self.free_energy and not self.complexes:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f"temp={self.temp}) \n \n    rxn_{self.rxn_id}.calculate_reaction_profile(free_energy=True) \n \n"
                )
            elif self.complexes and not self.free_energy:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f"temp={self.temp}) \n \n    rxn_{self.rxn_id}.calculate_reaction_profile(with_complexes=True) \n \n"
                )
            else:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f"temp={self.temp}) \n \n    rxn_{self.rxn_id}.calculate_reaction_profile() \n \n"
                )

        elif self.solvent != "gas":
            if self.free_energy and self.complexes:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f'temp={self.temp}, solvent_name="{self.solvent}") \n \n    rxn_{self.rxn_id}.calculate_reaction_profile(with_complexes=True, free_energy=True) \n \n'
                )

            elif self.free_energy and not self.complexes:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f'temp={self.temp}, solvent_name="{self.solvent}") \n \n    rxn_{self.rxn_id}.calculate_reaction_profile(free_energy=True) \n \n'
                )
            elif self.complexes and not self.free_energy:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f'temp={self.temp}, solvent_name="{self.solvent}") \n \n    rxn_{self.rxn_id}.calculate_reaction_profile(with_complexes=True) \n \n'
                )
            else:
                return (
                    f'    rxn_{self.rxn_id} = ade.Reaction(r"{self.rxn_smiles}", name="{self.folder_name}",'
                    f'temp={self.temp}, solvent_name="{self.solvent}") \n \n    rxn_{self.rxn_id}.calculate_reaction_profile() \n \n'
                )


class InputFile:
    """
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
    """

    def __init__(
        self,
        idx: int = None,
        n_cores: int = 16,
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

    def write_reaction_point(self, data_point, reaction_num):
        """Writes a reaction to the input file"""
        self.file.write("os.chdir(path) \n \n")
        self.file.write("try: \n")
        self.file.write(str(data_point))
        self.file.write(data_point.clean_up_folder())
        self.file.write("except Exception as e: \n")
        self.file.write(f'    print(f"Failure for reaction {reaction_num}: {{e}}") \n')
        self.file.write(
            f'    path_to_tar = os.path.join(os.path.join(path, "{self.autode_folder}"),str({reaction_num})) \n'
        )
        self.file.write(
            f'    tar_file = os.path.join(os.path.join(path, "{self.autode_folder}"),"{reaction_num}.tar.gz") \n \n'
        )
        self.file.write(
            f'    subprocess.run(["tar", "-zcf", f"{{tar_file}}",'
            f'f"{{path_to_tar}}"]) \n'
        )
        self.file.write(f'    subprocess.run(["rm", "-r", f"{{path_to_tar}}"]) \n \n')

    def close_input_file(self):
        """Closes the input file after all the individual data points have been written"""
        self.file.close()

    def initialize_input_file(self):
        """Writes the common lines to the autodE input script"""
        if self.functional is not None:
            low_opt, grad, opt, opt_ts, hess, sp = self.determine_keywords()

        self.file.write("import autode as ade \n")
        self.file.write("import subprocess \n")
        self.file.write("import os \n \n")
        self.file.write(f"ade.Config.n_cores = {self.n_cores} \n")
        self.file.write(f"ade.Config.hmethod_conformers = False \n \n")
        self.file.write(f"path = os.getcwd() \n \n")
        if self.functional is not None:
            self.file.write(
                f"ade.Config.G16.keywords.low_opt = ade.SinglePointKeywords({low_opt}) \n"
            )
            self.file.write(
                f"ade.Config.G16.keywords.grad = ade.SinglePointKeywords({grad}) \n"
            )
            self.file.write(
                f"ade.Config.G16.keywords.opt = ade.SinglePointKeywords({opt}) \n"
            )
            self.file.write(
                f"ade.Config.G16.keywords.opt_ts = ade.SinglePointKeywords({opt_ts}) \n"
            )
            self.file.write(
                f"ade.Config.G16.keywords.hess = ade.SinglePointKeywords({hess}) \n"
            )
            self.file.write(
                f"ade.Config.G16.keywords.sp = ade.SinglePointKeywords({sp}) \n \n"
            )

    def determine_keywords(self):
        """Determines and formats all the keywords necessary for the autodE calculation"""
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
