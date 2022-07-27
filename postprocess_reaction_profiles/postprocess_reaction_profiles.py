from argparse import ArgumentParser
from rdkit import Chem
import os
import pandas as pd
import numpy as np
import autode as ade
import subprocess
from autode.conformers import conf_gen
from pebble import ProcessPool
from lib import create_logger
import shutil
from autode.mol_graphs import find_cycles
import tempfile

hartree = 627.5094740631

parser = ArgumentParser()
parser.add_argument(
    "--data-file",
    type=str,
    required=True,
    help="data set file (.csv)",
)
parser.add_argument(
    "--xyz-folder", 
    type=str, 
    required=True,
    help="name (not path) of the folder with the (preliminary) reaction profiles"
)
parser.add_argument(
    "--n-cores", type=int, default=48, help="number of cores to be used"
)

def find_dipole(rxn_smiles):
    ''' Locate the dipole within the reaction SMILES '''
    reactant_list = rxn_smiles.split('>')[0].split('.')
    for reactant in reactant_list:
        if '+' in reactant and '-' in reactant:
            return reactant


def parse_coordinates(lines):
    ''' Turn lines from xyz-file into a list of coordinates '''
    coordinates = []
    for line in lines:
        _,x,y,z = line.split()
        coordinates.append(list(map(lambda x: float(x), [x,y,z])))
    return np.array(coordinates)


def extract_ts_geometry(r_smiles, p_smiles, rxn_id, home_dir, work_dir, geom_dir):
    ''' Wrapper around the TSGeomExtractor object (can take care of potential failures etc.) '''
    os.chdir(home_dir)
    geom_path = os.path.join(home_dir, geom_dir)
    return TSGeomExtractor(r_smiles, p_smiles, rxn_id, work_dir, geom_path)


def remove_atom_map_nums(smiles):
    ''' Remove atom-mapping from SMILES string'''
    mol = Chem.MolFromSmiles(smiles)
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
    return Chem.MolToSmiles(mol)


def get_ordering(mapped_mol, unmapped_mol):
    ''' Get the ordering of the atoms in the xyz-file '''
    ordering = unmapped_mol.GetSubstructMatch(mapped_mol)
    return ordering


def get_map_num_dict(smiles):
    ''' Get a dict {map_num: Idx} for a mapped SMILES string'''
    map_num_dict = {}
    for atom in Chem.MolFromSmiles(smiles).GetAtoms():
        map_num_dict[atom.GetAtomMapNum()] = atom.GetIdx() 
    return map_num_dict


def save_rotation_xyzs(home_path, work_dir, id, bond_num):
    ''' Get xyz-files of rotation profile id and write a subset of them '''
    with open(os.path.join(work_dir, f'xtbscan.log'), 'r') as f:
        lines = f.readlines()
        energy_indices = [lines.index(line) for line in lines if line.startswith(' e')]

    os.chdir(os.path.join(home_path, 'rotation_profiles'))        
    os.makedirs(f'profile_{id}_{bond_num}', exist_ok=True)
    for i in range(6):
        try:
            with open(f'profile_{id}_{bond_num}/{i}.xyz', 'w') as f:
                for index in range(energy_indices[i*10]-1, energy_indices[i*10+1]-1):
                    f.write(lines[index])
        except Exception as e:
            print(e)
            print(f'Saving of rotation xyzs failed for {id}')
    subprocess.run(['tar','cvzf', f'profile_{id}_{bond_num}.tar.gz', f'profile_{id}_{bond_num}'])
    shutil.rmtree(f'profile_{id}_{bond_num}')


def assess_bond_rotatability_from_scan(outfile):
    ''' 
    Get rotation barrier and bias energies from GFN2-XTB output files and determine whether rotation 
    has taken place succesfully and with a low barrier

    Parameters:
        outfile (Path): The path to the output file

    Returns:
        bool: True or False
    '''
    # get_data from output-file to decide whether rotation was succesful
    with open(outfile, 'r') as f:
        lines = f.readlines()
        energy_values = [float(line.split()[2]) for line in lines if line.startswith('unbiased energy:')]
        bias_energy_values = [float(line.split()[2]) for line in lines if line.startswith('    bias energy:')]

    # decide whether rotation has been succesful and low in energy:
    # (1) are energies the same once rotation is complete? 
    # (2) does the bias energy exceed 20 kcal/mol at any point? 
    # (3) is the barrier reasonable (> 2 kcal/mol and < 20 kcal/mol)?
    if abs(energy_values[0] - energy_values[-1]) * hartree > 5 or max(bias_energy_values) * hartree > 20 or \
     (max(energy_values) - min(energy_values)) * hartree > 20 or (max(energy_values) - min(energy_values)) * hartree < 2:
        return False
    else:
        return True


def write_scan_input_file(settings_path, dihedral_atoms, init_angle):
    ''' Write input file for GFN2-XTB scan '''
    with open(settings_path, 'w') as f:
        f.write('$constrain\n')
        f.write(' force constant=0.1\n')
        f.write(f' dihedral: {dihedral_atoms[0]+1},{dihedral_atoms[1]+1},{dihedral_atoms[2]+1},{dihedral_atoms[3]+1},{init_angle}\n')
        f.write('$end\n')
        f.write('$scan\n')
        f.write(f'1:{init_angle},{init_angle+360},60\n')
        f.write('$end')


def determine_rotatability(geometry_extractor, home_dir):
    ''' 
    Scan the dihedral angle with GFN2-xTB and determine whether bonds are rotatable

    Parameters:
        geometry_extractor: object containing geometries and data about individual datapoints

    Returns:
        unrotatable_bonds (list): list of bonds that are unrotatable
    '''
    unrotatable_bonds = []
    work_dir = geometry_extractor.work_dir

    # If no bond needs to be checked, return empty list
    if not geometry_extractor.dihedral_angles_to_track:
        return unrotatable_bonds

    # Iterate through bonds and perform scan
    for i, bond in enumerate(geometry_extractor.dihedral_angles_to_track):
        os.makedirs(os.path.join(work_dir, f'scan_{geometry_extractor.rxn_id}_{i}'), exist_ok=True)
        os.chdir(os.path.join(work_dir, f'scan_{geometry_extractor.rxn_id}_{i}'))
        input_xyz_scan = os.path.join(os.path.join(work_dir, f'scan_{geometry_extractor.rxn_id}_{i}'), 
                        f'input_scan_{geometry_extractor.rxn_id}_{i}.xyz')
        shutil.copy(geometry_extractor.reactive_opt_xyz_path, input_xyz_scan)
        outfile = f'{os.path.join(os.path.join(work_dir, f"scan_{geometry_extractor.rxn_id}_{i}"), str(geometry_extractor.rxn_id))}_{i}.out'
        settings_path = os.path.join(os.path.join(work_dir, f"scan_{geometry_extractor.rxn_id}_{i}"), 'scan.inp')
        init_angle = geometry_extractor.ade_ts_mol.dihedral(bond[0], bond[1], bond[2], bond[3]).to('degrees')

        write_scan_input_file(settings_path, bond, init_angle)

        with open(outfile, "w") as out:
            subprocess.run(f'xtb {input_xyz_scan} --alpb water --opt --input {settings_path}',
            shell=True, stdout=out, stderr=out,)

        rotatable = assess_bond_rotatability_from_scan(outfile)
        save_rotation_xyzs(home_dir, os.path.join(work_dir, f'scan_{geometry_extractor.rxn_id}_{i}'), geometry_extractor.rxn_id, i)
        if not rotatable:
            unrotatable_bonds.append(bond)

        shutil.copy(outfile, os.path.join(home_dir, f'output_folder/{outfile.split("/")[-1]}'))

    os.chdir(work_dir)

    for i in range(len(geometry_extractor.dihedral_angles_to_track)):
        shutil.rmtree(os.path.join(work_dir, f'scan_{geometry_extractor.rxn_id}_{i}'))

    return unrotatable_bonds


def get_dihedral_angle_atoms(atom1, atom2):
    ''' Determine neighbors of the bonding partners to follow during dihedral angle determination '''
    try:
        neighbor_atom1 = [atom for atom in atom1.GetNeighbors() if (atom.GetIdx() != atom2.GetIdx() and atom.GetAtomMapNum() != 0)][0]
        neighbor_atom2 = [atom for atom in atom2.GetNeighbors() if (atom.GetIdx() != atom1.GetIdx() and atom.GetAtomMapNum() != 0)][0]
        return [neighbor_atom1, atom1, atom2, neighbor_atom2]
    except IndexError:
        # for triple bonds and O-sites, this will generally fail
        return None


def generate_and_select_conformers(geometry_extractor, geom_constraints, home_dir):
    ''' Generate new (RR) conformers (with geometry constraints applied) with autodE and select lowest energy one '''
    dipole_mol = geometry_extractor.ade_ts_mol
    dipole_mol.name = f'{geometry_extractor.min_xyz_path.split("/")[-1].split(".")[0]}_alt'
    distance_constraints = {}

    # constrain dihedral angles by fixing the bond lengths of 3 consecutive bonds 
    # as well as the distance from the first until the last atom of the series
    for dihedral in geom_constraints:
        distance_constraints[(dihedral[0], dihedral[1])] = dipole_mol.distance(dihedral[0], dihedral[1])
        distance_constraints[(dihedral[1], dihedral[2])] = dipole_mol.distance(dihedral[1], dihedral[2])
        distance_constraints[(dihedral[2], dihedral[3])] = dipole_mol.distance(dihedral[2], dihedral[3])
        distance_constraints[(dihedral[0], dihedral[3])] = dipole_mol.distance(dihedral[0], dihedral[3])

    # generate constrained conformers
    conformer_gen_dir = os.path.join(os.getcwd(), 'siman_trash')
    os.makedirs(conformer_gen_dir, exist_ok=True)
    os.chdir(conformer_gen_dir)
    for n in range(1000):
        conformer = conf_gen.get_simanl_conformer(dipole_mol, dist_consts=distance_constraints, conf_n=n)
        dipole_mol.conformers.append(conformer)

    # go through pruning steps and optimize with XTB
    dipole_mol.conformers.prune_on_energy(e_tol=1E-10)
    dipole_mol.conformers.prune_on_rmsd()
    dipole_mol.conformers.optimise(method=ade.methods.XTB())
    dipole_mol.conformers.prune(remove_no_energy=True)
    dipole_mol.conformers.prune_diff_graph(dipole_mol.graph)
    
    # only replace the original compatible conformer if one of the newly generated ones is lower in energy
    if dipole_mol.energy > dipole_mol.conformers.lowest_energy.energy:
        dipole_mol._set_lowest_energy_conformer()

    # copy to final output directory and clean up
    path = os.path.join(os.path.join(home_dir, 'lowest_energy_conformers_RR'), str(geometry_extractor.rxn_id))
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    shutil.rmtree(conformer_gen_dir)

    # print the final geometry of the conformer
    dipole_mol.print_xyz_file()
    os.chdir(home_dir)

    return dipole_mol.name, (dipole_mol.energy.to('kcalmol') - geometry_extractor.ade_r_mol.energy.to('kcalmol'))


class TSGeomExtractor:
    '''
    Class responsible for extracting geometry of dipole from TS xyz-file and preparing it for dihedral scan

    Attributes:
        r_smiles (str): the SMILES string of the dipole
        p_smiles (str): the SMILES string of the full addition product
        rxn_id (int): the reaction ID of the reaction considered
        geom_path (str): the path to the directory containing all the xyz-files

    Methods:
        get_reactant_xyz_file_name: determines the path to the dipole xyz-file
        write_xyz_file: writes xyz-file for the dipole geometry extracted from the TS-xyz file
        extract_lines_from_ts_xyz: extracts lines from TS-xyz file which correspond to the reactant dipole
        optimize_geom: optimizes the TS and reactant geometry of the dipole from its xyz-files with autodE (returns opt xyz-files and autodE mol-objects)
        get_dihedral_angles_to_track: determine dihedral angles which should be checked for rotatability

    '''
    def __init__(self, r_smiles, p_smiles, rxn_id, work_dir, geom_path):
        # set the continue_workflow flag to True 
        self.continue_workflow = True

        # if geom_path does not exist, break off immediately
        if not os.path.isdir(geom_path):
            self.continue_workflow == False
        else:
            # initialize
            self.r_smiles = r_smiles
            self.p_smiles = p_smiles
            self.rxn_id = rxn_id

            self.work_dir = os.path.join(work_dir, f'xtb_trash/{self.rxn_id}')
            os.makedirs(self.work_dir, exist_ok=True)
            os.chdir(self.work_dir)

            # store all relevant properties of the dipole in the reactant geometry
            self.r_mol = Chem.AddHs(Chem.MolFromSmiles(self.r_smiles))
            self.min_xyz_path = self.get_reactant_xyz_file_name(geom_path)
            self._r_mol_unmapped = Chem.AddHs(Chem.MolFromSmiles(remove_atom_map_nums(self.r_smiles)))
            self._r_mol_map_nums = [atom.GetAtomMapNum() for atom in Chem.MolFromSmiles(self.r_smiles).GetAtoms()]
            self._r_mol_ordering = get_ordering(self.r_mol, self._r_mol_unmapped)

            # store all relevant properties of the dipole in the TS geometry
            self.ts_mol = Chem.AddHs(Chem.MolFromSmiles(self.p_smiles))
            self._ts_mol_unmapped = Chem.AddHs(Chem.MolFromSmiles(remove_atom_map_nums(self.p_smiles)))
            self._ts_mol_ordering = get_ordering(self.ts_mol, self._ts_mol_unmapped)
            self._map_num_dict = get_map_num_dict(self.p_smiles)
            self._retained_lines = self.extract_lines_from_ts_xyz(geom_path)
            self._reactive_xyz_path = os.path.join(self.work_dir, f'reactive_dipole_geom_{self.rxn_id}.xyz')

            # optimize and store dihedral angles to track
            self.write_xyz()
            self.reactive_opt_xyz_path, self.ade_r_mol, self.ade_ts_mol = self.optimize_geom()
            if self.ade_r_mol is not None:
                self.dihedral_angles_to_track = self.get_dihedral_angles_to_track()
            else:
                self.continue_workflow = False
                self.dihedral_angles_to_track = None


    def get_reactant_xyz_file_name(self, geom_path):
        ''' Determine the path to the dipole xyz-file, copy it the work directory and return the final path '''
        file_list = [file for file in os.listdir(os.path.join(geom_path, str(self.rxn_id))) if (file.startswith('r') and file.endswith('.xyz') and '_opt' not in file)]
        for file_name in file_list:
            with open(os.path.join(os.path.join(geom_path, str(self.rxn_id)),file_name), 'r') as f:
                lines = f.readlines()
                if int(lines[0]) == self.r_mol.GetNumAtoms() and len([line for line in lines if not line.startswith('H')]) - 2 == self.r_mol.GetNumHeavyAtoms() \
                and len([line for line in lines if line.startswith('C')]) == len([atom for atom in self.r_mol.GetAtoms() if atom.GetSymbol() == 'C']) \
                and len([line for line in lines if line.startswith('N')]) == len([atom for atom in self.r_mol.GetAtoms() if atom.GetSymbol() == 'N']) \
                and len([line for line in lines if line.startswith('O')]) == len([atom for atom in self.r_mol.GetAtoms() if atom.GetSymbol() == 'O']):
                    original_path = os.path.join(os.path.join(geom_path, str(self.rxn_id)),file_name)
                    destination_path = os.path.join(self.work_dir, file_name)
                    shutil.copy(original_path, destination_path)

                    return destination_path

    def write_xyz(self):
        ''' Writes xyz-file for the dipole geometry extracted from the TS-xyz file '''
        retained_lines = self._retained_lines
        with open(f'{self._reactive_xyz_path}', 'w') as f:
            f.write(f'{len(retained_lines)} \n')
            f.write(f'{self._reactive_xyz_path}\n')
            for line in retained_lines:
                f.write(line)             

    def extract_lines_from_ts_xyz(self, geom_path):
        ''' Extracts lines from TS-xyz file which correspond to the reactant dipole '''
        folder_path = os.path.join(geom_path, str(self.rxn_id))
        file = [file for file in os.listdir(folder_path) if file.startswith('TS')][0]
        # First, translate from map numbers to (mapped) TS atom indices
        indices_ts_mol_to_retain = [self._map_num_dict[map_num] for map_num in self._r_mol_map_nums]
        # Second, translate atom indices from mapped TS to unmapped TS
        line_indices_to_retain = [self._ts_mol_ordering[index] + 2 for index in indices_ts_mol_to_retain]

        # Get those lines from TS xyz-file
        retained_lines = []
        with open(os.path.join(folder_path,file), 'r') as f:
            lines = f.readlines()
            for index in line_indices_to_retain:
                retained_lines.append(lines[index])
        # Then, reorder the lines so that they match the ordering of the reactant
        retained_lines_reordered = [i for i in range(len(retained_lines))]
        for i, line in enumerate(retained_lines):
            retained_lines_reordered[self._r_mol_ordering[i]] = line

        # Get the associated coordinates
        heavy_atom_coordinates = parse_coordinates(retained_lines_reordered)

        # Finally, add the H-atom lines
        H_coordinate_lines = [line for line in lines if line.startswith('H')]
        H_coordinates = parse_coordinates(H_coordinate_lines)

        for heavy_atom_coordinate in heavy_atom_coordinates:
            for i, H_coordinate in enumerate(H_coordinates):
                if np.linalg.norm(H_coordinate - heavy_atom_coordinate) < 1.30:
                    retained_lines_reordered.append(H_coordinate_lines[i])

        return retained_lines_reordered

    def optimize_geom(self):
        ''' Optimize the TS and reactant geometry of the dipole from its xyz-files with autodE (returns opt xyz-files and autodE mol-objects) '''
        # set up ade_r_mol object
        ade_r_mol = ade.Molecule(self.min_xyz_path, name=self.min_xyz_path.rstrip(".xyz"), solvent_name='water')
        ade_r_mol.optimise(method=ade.methods.XTB())
        ade_r_mol.reset_graph()

        # if ade_r_mol does not have the same number of cycles as suggested by the SMILES connectivity, abort
        ade_r_smiles_mol = ade.Molecule(name='r_smiles_mol', smiles= self.r_smiles)
        if not len(find_cycles(ade_r_mol.graph)) == len(find_cycles(ade_r_smiles_mol.graph)):
            return None, None, None

        # set up ade_ts_mol object
        ade_ts_mol = ade.Molecule(self._reactive_xyz_path, name=self._reactive_xyz_path.rstrip(".xyz"), solvent_name='water')
        ade_ts_mol.optimise(method=ade.methods.XTB())
        ade_ts_mol.reset_graph()

        # if connectivity of ade_r_mol and ade_ts_mol do not match, check if the atoms are the same, and if so, work with the non-optimized geometry
        # else, abort
        if len(find_cycles(ade_r_mol.graph)) == len(find_cycles(ade_ts_mol.graph)):
            return f'{ade_ts_mol.name}_optimised_xtb.xyz', ade_r_mol, ade_ts_mol
        else: 
            if len(ade_r_mol.atoms) == len(ade_ts_mol.atoms) and [atom.label for atom in ade_r_mol.atoms] == [atom.label for atom in ade_ts_mol.atoms]:
                return self._reactive_xyz_path, ade_r_mol, ade_ts_mol
            else:
                print(self.rxn_id, ade_ts_mol.atoms, ade_r_mol.atoms)
                return None, None, None

    def get_dihedral_angles_to_track(self):
        ''' Determine dihedral angles which should be checked for rotatability '''
        dihedral_angles_to_track = []

        # First, find the positive and negative poles of the dipole
        for atom in self.r_mol.GetAtoms():
            if atom.GetFormalCharge() == -1:
                neg_charged_atom = atom
            elif atom.GetFormalCharge() == 1:
                pos_charged_atom = atom

        # Then, find the final atom
        for bond in pos_charged_atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                if bond.GetBeginAtom().GetFormalCharge() == 0:
                    final_atom = bond.GetBeginAtom()
                elif bond.GetEndAtom().GetFormalCharge() == 0:
                    final_atom = bond.GetEndAtom()
            # If there is a triple bond, or the bond is aromatic (i.e., in a ring), then there is nothing to be constrained 
            # -> break off the execution of the function
            if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE or bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                return dihedral_angles_to_track
        
        # Finally, add the two bonds, i.e., final_atom -> pos_pole and pos_pole -> min_pole, to the dihedral_angles_to_track list
        for atom1, atom2 in [(final_atom, pos_charged_atom), (pos_charged_atom, neg_charged_atom)]:
            try:
                dihedral_angles_to_track.append(list(map(lambda x: self._r_mol_ordering[x.GetIdx()], get_dihedral_angle_atoms(atom1, atom2))))
            except TypeError:
                pass

        return dihedral_angles_to_track


def get_input_pool(df, home_dir, tmp_dir, geom_dir):
    ''' Extract relevant input from each row in the dataframe and return in list format '''
    df['dipole_smiles'] = df['rxn_smiles'].apply(lambda x: find_dipole(x))
    df['product_smiles'] = df['rxn_smiles'].apply(lambda x: x.split('>')[-1])

    dipole_smiles_list = df.dipole_smiles.values.tolist()
    product_smiles_list = df.product_smiles.values.tolist()
    rxn_id_list = df.index.tolist()
    home_dir_list = [home_dir for i in range(len(rxn_id_list))]
    tmp_dir_list = [tmp_dir for i in range(len(rxn_id_list))]
    geom_dir_list = [geom_dir for i in range(len(rxn_id_list))]

    input_list = list(zip(dipole_smiles_list, product_smiles_list, rxn_id_list, home_dir_list, tmp_dir_list, geom_dir))

    return input_list


def get_compatible_conformer(input_tuple):
    ''' Generate a stereo-compatible conformer for the dipole '''
    # extract all info necessary to perform 
    dipole_smiles, product_smiles, rxn_id, home_dir, tmp_dir, geom_dir = input_tuple
    geometry_extractor = extract_ts_geometry(dipole_smiles, product_smiles, rxn_id, home_dir, tmp_dir, geom_dir)

    if not geometry_extractor.continue_workflow:
        return dipole_smiles, rxn_id, None, None, None

    geom_constraints = determine_rotatability(geometry_extractor, home_dir)
    name, energy_difference = generate_and_select_conformers(geometry_extractor, geom_constraints, home_dir)

    if len(geom_constraints) == 0:
        # If no geometry constraints, then simply check if the generated XTB conformer is lower
        # in energy (by at least 0.1 kcal/mol) than the originally registered dipole reactant conformer
        if energy_difference < -0.1:
            to_run, constrained = True, False
        else:
            to_run, constrained = False, False
    else:
        # If geometry constraints, then always run the alternative conformer
        to_run, constrained = True, True

    return name, rxn_id, energy_difference, to_run, constrained


def get_all_compatible_conformers(csv_file, geom_dir, num_cores):
    ''' Main function which controls the process pool to obtain stereo-compatible conformers for the all the dipoles '''
    # Set up the logger
    logger_name = f"{csv_file.split('/')[-1].rstrip('.csv')}_dipole_geom_compatibility"
    logger = create_logger(logger_name)
    pwd = os.getcwd()
    tmp_dir = tempfile.gettempdir()

    # Prepare the input
    df = pd.read_csv(csv_file)
    input_list = get_input_pool(df, pwd, tmp_dir, geom_dir)[5738:5740]

    # Get everything related to the output set up
    os.makedirs(os.path.join(pwd, 'rotation_profiles'), exist_ok=True)
    os.makedirs(os.path.join(pwd, 'output_folder'), exist_ok=True)
    os.makedirs(os.path.join(tmp_dir, 'xtb_trash'), exist_ok=True)
    os.makedirs(os.path.join(pwd, 'lowest_energy_conformers_RR'), exist_ok=True)

    output_list = []

    # Process data points in parallell with the help of a process pool
    with ProcessPool(max_workers=num_cores) as pool:
        future = pool.map(get_compatible_conformer, input_list, timeout = 10000)
        iterator = future.result()

        while True:
            try:
                name, rxn_id, energy_difference, to_run, constrained = next(iterator)
                if energy_difference != None:
                    # A positive energy difference means that the original conformer was lower in energy
                    logger.info(f'Energy difference for {name}, rxn_id {rxn_id} : {energy_difference}')
                    output_list.append([name, rxn_id, energy_difference, to_run, constrained])
                else:
                    logger.info(f'Workflow got interrupted for rxn_id {rxn_id} due to issues with geometries')
            
            except StopIteration:
                break
            except TimeoutError as error:
                logger.info(f'get_compatible_conformers call took more than {error.args} seconds')
                raise
            except Exception as error:
                logger.info(f'Failure: {error}')
                os.chdir(pwd)
        
        pool.close()
        pool.join()

    with open(os.path.join(pwd, 'output_postprocessing_dipole_conformers.csv'), 'w') as f:
        f.write(',name,rxn_id,energy_difference,to_run,constrained\n')
        for i, output in enumerate(output_list):
            f.write(f'{i}, {output[0]}.xyz,{output[1]},{output[2]},{output[3]},{output[4]}\n')

if __name__ == '__main__':
    args = parser.parse_args()
    get_all_compatible_conformers(args.data_file, args.xyz_folder, args.n_cores)
