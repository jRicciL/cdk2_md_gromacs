import subprocess
import os 
import shutil
from pathlib import Path
import warnings
from typing import List, Tuple
import inspect
import numpy as np
from decimal import Decimal


SOLVENTS_LIB_DIR = '../data/cosolvent_libs/'

def run_pdb2pqr(input_pdb: str, 
                output_basename: str,
                ph: float = 7.4,
                verbose: str = 'error'):

    pdb2pqr_command = \
    f"""
    pdb2pqr30 --ff='AMBER' --ffout='AMBER'\
        --with-ph={ph} -o={ph}\
        --drop-water --keep-chain\
        --log-level='INFO'\
        --titration-state-method='propka'\
        --pdb-output {output_basename}.pdb\
        {input_pdb} {output_basename}.pqr
    """

    ps = subprocess.Popen(pdb2pqr_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    if verbose:
        print(output)
        print(error)


def run_obabel(input_ligand: str,
           output_name: str,
           ph: float = 7.0,
           partial_charges: str = 'gasteiger',
           use_amber_reduce: bool = True,
           verbose: bool = True):

    if use_amber_reduce:
        amber_reduce_command = 'reduce -NUClear -OH -ROTNH3 -ALLALT ' +\
                                 f'{input_ligand} > {output_name}.amber_reduce.pdb'
        input_ligand = f'{output_name}.amber_reduce.pdb'
        p = subprocess.Popen(amber_reduce_command, 
                             shell  = True, 
                             stdout = subprocess.PIPE,
                             stderr = subprocess.STDOUT,
                             encoding = 'UTF-8')
        p.communicate()

    # Save mol2
    obabel_command = 'obabel ' + \
                     f'-ipdb {input_ligand} ' + \
                     f'-omol2 -O {output_name} ' + \
                     f'--partialcharge "{partial_charges}" ' + \
                     ('' if use_amber_reduce else f'-p {ph} ')
                     # f'--minimize --steps 1500 --sd --ff "GAFF" '

    ps = subprocess.Popen(obabel_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    ps = subprocess.Popen(obabel_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    if verbose:
        print(output)
        print(error)



def run_get_charge(filename: str, round_value: bool = True) -> int:
    get_charge_command = "sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/p' " + \
                            f"{filename}" + \
                            " | sed '1d;$d' | awk '{ SUM += $NF} END {print SUM}'"
    ps = subprocess.Popen(get_charge_command, shell = True, 
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                          encoding='UTF-8')
    output = ps.communicate()[0]
    charge = float(output.replace('\n', ''))
    if round_value:
        LIG_NET_CHARGE = round(charge)
    else:
        LIG_NET_CHARGE = charge
    return LIG_NET_CHARGE


def correct_non_integer_charges(mol2_file: str):
    LIG_FLOAT_CHARGE = run_get_charge(
                            filename = mol2_file, 
                            round_value = False)
    
    print('Current charge:', LIG_FLOAT_CHARGE)
    with open(mol2_file, 'r') as f:
        lines = f.readlines()
        # Get atom files
        atom_lines = []
        start = lines.index('@<TRIPOS>ATOM\n')
        end   = lines.index('@<TRIPOS>BOND\n')
        atom_lines = lines[start + 1: end]

        # get charges
        charges = np.array([float(l.split()[-1]) 
                            for l in atom_lines])

        closest_int = round(LIG_FLOAT_CHARGE)
        frac = round(float(Decimal(LIG_FLOAT_CHARGE - closest_int)), 6)
            
        # Check if the charge is positive or negative
        if closest_int < 0:
            max_idx   = np.argmax(charges)
            new_value = round(Decimal(charges[max_idx] - frac), 6)
        else:
            max_idx = np.argmin(charges)
            new_value = round(Decimal(charges[max_idx] - frac), 6)
        # Update charges
        charges[max_idx] = new_value
        charges_str = [str(round(Decimal(c), 6)) + '\n' 
                       for c in charges]
        charges_str = [s if s[0] == '-' else ' ' + s 
                       for s in charges_str]
        # update atom_lines
        new_atom_lines = [l[:-10] + c 
                          for l, c in zip(atom_lines, charges_str)]
        first_block = lines[:start + 1]
        end_block   = lines[end: ]
        final_lines = first_block + new_atom_lines + end_block

    with open(mol2_file, 'w') as f:
        f.writelines(final_lines)



def run_parmchk2(mol2_filename: str, 
                 frcmod_filename: str,
                 verbose: bool = True):
    parmchk2_command = f'parmchk2 -i {mol2_filename} -f mol2 -o {frcmod_filename} -s gaff'

    ps = subprocess.Popen(parmchk2_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    if verbose:
        print(output)
        print(error)


def run_antechamber(mol2_filename: str, 
                    lig_resname: str,
                    lig_net_charge: int,
                    output_dir: str = '.',
                    verbose: bool = True):
    antechamber_command = f"antechamber -i {mol2_filename}" + \
                f" -fi mol2 -o {output_dir}/LIG.mol2" + \
                 " -fo mol2 -c bcc -s 2 -at gaff -eq 2" + \
                f" -rn {lig_resname}" + \
                f" -nc {lig_net_charge}"

    ps = subprocess.Popen(antechamber_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    if verbose:
        print(output)
        print(error)
        

    # Move the intermediate files
    move_temp_files_command = f'mv -v ANTECHAMBER* ATOMTYPE.INF sqm* {output_dir}'
    ps = subprocess.Popen(move_temp_files_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    o, e = ps.communicate()

    # Check partial charges
    correct_non_integer_charges(f"{output_dir}/LIG.mol2")

    # Run PARMCHK2
    run_parmchk2(mol2_filename = f'{output_dir}/LIG.mol2', 
                 frcmod_filename = f'{output_dir}/LIG.frcmod')

    


def run_leap_lig_lib(tmp_dir: str, 
                     lig_basename: str = 'LIG',
                     verbose: bool = True):
    frcmod = f'{tmp_dir}/{lig_basename}.frcmod'
    mol2   = f'{tmp_dir}/{lig_basename}.mol2'
    lib    = f'{tmp_dir}/{lig_basename}.lib'
    
    assert os.path.isfile(frcmod)
    assert os.path.isfile(mol2)

    leap_lig_lib = \
    f"""
    source leaprc.gaff
    loadamberparams {frcmod}
    LIG = loadmol2 {mol2}
    saveoff LIG {lib}
    quit
    """

    with open(f'{tmp_dir}/leap_lig_lib.in', 'w') as f:
        f.write(leap_lig_lib)

    # Run tleap
    tleap_command_lib = f'tleap -f {tmp_dir}/leap_lig_lib.in'
    ps = subprocess.Popen(tleap_command_lib, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()
    # Move temp files
    os.rename('leap.log', f'{tmp_dir}/leap_lig_lib.in.log')

    if verbose:
        print(output)
        print(error)


def _add_tleap_head_lines(input_ligand_basename: str = None,
                          input_waters_pdb: str = None,
                          solvent_type: str = None,
                          tmp_dir: str = None):

    SOLVENTS_DICT = {'WAT'  : 'TIP3PBOX', # Water
                     'ETA20': 'ETAWAT20', # Ethanol 20%
                     'ISO20': 'ISOWAT20', # Isopropanol 20% 
                     'MAM20': 'MAMWAT20'  # Acetamide 20%
                }

    _ligand_lines = ''
    _crys_wat_lines = ''
    _solvent_type_lines = ''

    if input_ligand_basename != None:
        _ligand_lines = inspect.cleandoc(f"""
        # Ligand parameters
        loadOff {tmp_dir}/{input_ligand_basename}.lib 
        loadamberparams {tmp_dir}/{input_ligand_basename}.frcmod
        ligand  = loadmol2 {tmp_dir}/{input_ligand_basename}.mol2
        """)
        
    if input_waters_pdb != None:
        _crys_wat_lines = inspect.cleandoc(f"""
        # Cocrystalized waters
        waters = loadpdb {input_waters_pdb}
        """)

    if solvent_type != 'TIP3PBOX':
        assert os.path.isfile(f'{SOLVENTS_LIB_DIR}/{solvent_type}.off'), f'Please provide the `{solvent_type}`.off file'
        _solvent_type_lines = f'loadoff {SOLVENTS_LIB_DIR}/{solvent_type}.off'

    head_lines = inspect.cleandoc(
        f"""
        source leaprc.protein.ff14SB
        source leaprc.gaff
        source leaprc.water.tip3p
        loadamberparams frcmod.ionsjc_tip3p
        {_ligand_lines}
        {_crys_wat_lines}
        {_solvent_type_lines}
        """)

    return head_lines


def run_tleap_prepare_system(input_protein_pdb: str,
                             output_basename: str, 
                             tmp_dir: str = '.', 
                             input_ligand_basename: str = None,
                             input_waters_pdb: str = None,
                             verbose: bool = True):

    _head_lines = _add_tleap_head_lines(
                            input_ligand_basename = input_ligand_basename,
                            input_waters_pdb = input_waters_pdb,
                            solvent_type = 'TIP3PBOX'
                            )

    _system_components = ['protein']
    if input_ligand_basename != None:
        _system_components.append('ligand')
    if input_waters_pdb != None:
        _system_components.append('waters')
    _system_line = 'system = combine { ' + \
                    ' '.join(_system_components) + \
                   ' }'

    leap_prepare_system = \
    f"""
    {_head_lines}

    protein = loadpdb {tmp_dir}/{input_protein_pdb}
    
    {_system_line}

    # Sanity check of the system
    check system

    # Check charges
    charge protein
    charge system

    # Save coords and params
    savepdb system {tmp_dir}/{output_basename}.pdb
    saveamberparm system {tmp_dir}/{output_basename}.prmtop {tmp_dir}/{output_basename}.rst7
    quit
    """

    with open(f'{tmp_dir}/leap_prep_{output_basename.lower()}.in', 'w') as f:
        f.write(leap_prepare_system)

    # Run tleap
    tleap_command_pl = f'tleap -f {tmp_dir}/leap_prep_{output_basename.lower()}.in'
    ps = subprocess.Popen(tleap_command_pl, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    # Move temp files
    os.rename('leap.log', f'{tmp_dir}/leap_prep_{output_basename}.in.log')

    if verbose:
        print(output)
        print(error)


def run_tleap_prepare_pl_complex(output_basename: str, 
                              tmp_dir: str, 
                              input_waters_pdb: str = None,
                              verbose: bool = True):
    _add_crys_wat_line = ''
    _system_line = 'system = combine {protein ligand}'

    if input_waters_pdb != None:
        _add_crys_wat_line = f'waters = loadpdb {input_waters_pdb}'
        _system_line = 'system = combine {protein waters ligand}'

    leap_prepare_system = \
    f"""
    source leaprc.protein.ff14SB
    source leaprc.gaff
    source leaprc.water.tip3p
    loadOff {tmp_dir}/LIG.lib 
    loadamberparams {tmp_dir}/LIG.frcmod

    protein = loadpdb {tmp_dir}/prot.TEMP.pdb
    ligand  = loadmol2 {tmp_dir}/LIG.mol2
    {_add_crys_wat_line}
    {_system_line}

    savepdb system {tmp_dir}/{output_basename}.pdb
    # Sanity check of the system
    check system

    saveamberparm system {tmp_dir}/{output_basename}.prmtop {tmp_dir}/{output_basename}.rst7
    # Protein charge
    charge protein
    # Ligand charge
    charge ligand
    # System charge
    charge system
    quit
    """

    with open(f'{tmp_dir}/leap_prep_{output_basename.lower()}.in', 'w') as f:
        f.write(leap_prepare_system)

    # Run tleap
    tleap_command_pl = f'tleap -f {tmp_dir}/leap_prep_{output_basename.lower()}.in'
    ps = subprocess.Popen(tleap_command_pl, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    # Move temp files
    os.rename('leap.log', f'{tmp_dir}/leap_prep_{output_basename}.in.log')

    if verbose:
        print(output)
        print(error)




def get_sytem_charges(tleap_log_file: str, 
                      mol_names: List = ['system']):
    assert os.path.isfile(tleap_log_file)
                          
    def get_charge(line: str):
        charge = float(line.strip().split(' ')[-1])
        return round(charge)
                          
    charges = {}
    with open(tleap_log_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            for mol in mol_names:
                if f'charge {mol}' in line:
                    mol_charge = get_charge(lines[i + 1])
                    charges[mol] = mol_charge
                    
    if not bool(charges):
        warnings.warn(f"No molecule definitions found in {tleap_log_file}") 
    
    return charges



def run_minimize_hydrogens(input_basename: str,
                           output_basename: str, 
                           tmp_dir: str,
                           verbose: bool = True):
    assert (input_basename != output_basename)

    sander_mdrun_min_H = \
    f"""\
    Taken from  biobb_amber module from the BioBB library
    Type of mdin: min_vacuo
    &cntrl
        imin = 1 
        ncyc = 250 
        ntb  = 0 
        igb  = 0 
        cut  = 12 
        maxcyc = 500 
        ntpr = 5 
        ntr  = 1 
        restraintmask = ":*&!@H=" 
        restraint_wt = 500.0 
    &end
    """

    with open(f'{tmp_dir}/sander_mdrun_min_H.mdin', 'w') as f:
        f.write(sander_mdrun_min_H)

    # Run Sander
    sander_minH_command = \
    f"""
    sander -O \
    -i {tmp_dir}/sander_mdrun_min_H.mdin \
    -p {tmp_dir}/{input_basename}.prmtop \
    -c {tmp_dir}/{input_basename}.rst7 \
    -ref {tmp_dir}/{input_basename}.rst7 \
    -r {tmp_dir}/{output_basename}.rst7 \
    -o {tmp_dir}/sander.minH.log \
    -inf {tmp_dir}/sander.minH.mdinfo \
    -x {tmp_dir}/{output_basename}.x 
    """

    ps1 = subprocess.Popen(sander_minH_command, 
                           shell  = True, 
                           stdout = subprocess.PIPE,
                           stderr = subprocess.STDOUT,
                           encoding = 'UTF-8')
    output1, error1 = ps1.communicate()

    if verbose:
        print(output1)
        print(error1)
    
    
    # Convert output files to PDB
    amber_pdb_command = 'ambpdb ' + \
            f'-p {tmp_dir}/{input_basename}.prmtop ' +\
            f'-c {tmp_dir}/{output_basename}.rst7 ' +\
            f'>  {tmp_dir}/{output_basename}.pdb'

    ps2 = subprocess.Popen(amber_pdb_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output2, error2 = ps2.communicate()




def run_leap_solv(input_basename: str, 
                  output_basename: str, 
                  box_padding: float,
                  solvent_type: str,
                  tmp_dir: str = '.',
                  input_ligand_basename: str = None,
                  verbose: bool = True):
    
    _head_lines = _add_tleap_head_lines(
                            input_ligand_basename = input_ligand_basename,
                            input_waters_pdb = None,
                            solvent_type = solvent_type,
                            tmp_dir = tmp_dir
                            )

    leap_solv = f"""
    {_head_lines}

    system = loadpdb {tmp_dir}/{input_basename}.pdb

    solvateOct system {solvent_type} {box_padding}
    savepdb system {tmp_dir}/{output_basename}.pdb

    # Sanity check of the system
    check system

    charge system

    saveamberparm system {tmp_dir}/{output_basename}.prmtop {tmp_dir}/{output_basename}.rst7
    quit
    """

    with open(f'{tmp_dir}/leap_prep_{output_basename}.in', 'w') as f:
        f.write(leap_solv)

    # Run tleap
    tleap_command_solv = f'tleap -f {tmp_dir}/leap_prep_{output_basename}.in'
    ps = subprocess.Popen(tleap_command_solv, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()
    # Move temp files
    os.rename('leap.log', f'{tmp_dir}/leap_prep_{output_basename}.in.log')

    if verbose:
        print(output)
        print(error)

    

def run_count_water_molecules(PDB_FILE_PATH: str):
    assert os.path.isfile(PDB_FILE_PATH)

    # Count the number of water moleculer
    count_wat_command = f'grep WAT {PDB_FILE_PATH} | wc -l'
    ps = subprocess.Popen(count_wat_command, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()
    N_WATS = round(float(output.replace('\n', '')) / 3)
    
    return N_WATS


def get_tleap_neutralization_line(ION_MOLAR: float, 
                                  N_WATS: int,
                                  sys_CHARGE: int) -> str:
    if (ION_MOLAR > 0) and (ION_MOLAR < 1):
        _WAT_MOL_MASS = 18.02 # g/mol
        _WAT_DENSITY  = 0.997 # kg/m**3
        # Number of water mols in 1 L
        _N_WAT_MOLS_1L = _WAT_DENSITY / _WAT_MOL_MASS * 1000
        N_IONS = int(N_WATS / _N_WAT_MOLS_1L * ION_MOLAR)

        #print(N_IONS, ION_MOLAR, N_WATS, sys_CHARGE)

        # Number of ions to get the required MOLAR concentration
        N_Na = N_IONS
        N_Cl = N_IONS
        # Add extra ions to neutralize the system
        if sys_CHARGE > 0:
            N_Cl += sys_CHARGE
        else:
            N_Na += abs(sys_CHARGE)
        _neutralize_line = f'addionsRand system Na+ {N_Na} Cl- {N_Cl}'

    else:
        ION_NEUTRALIZE =  'Cl-' if sys_CHARGE > 0 else 'Na+'
        _neutralize_line = f'addions system {ION_NEUTRALIZE} 0'

    return _neutralize_line


def run_leap_neutralization(input_basename: str, 
                            output_basename: str,
                            tmp_dir: str,
                            ion_concentration: float,
                            solvent_type: str,
                            input_ligand_basename: str = None,
                            verbose: bool = True):

    _head_lines = _add_tleap_head_lines(
                            input_ligand_basename = input_ligand_basename,
                            input_waters_pdb = None,
                            solvent_type = solvent_type,
                            tmp_dir = tmp_dir
                            )
    
    sys_charges = get_sytem_charges(
        tleap_log_file = f'{tmp_dir}/leap_prep_{input_basename}.in.log',
        mol_names = ['system']
    )
    sys_CHARGE = sys_charges['system']

    n_waters = run_count_water_molecules(f'{tmp_dir}/{input_basename}.pdb')

    _neutralize_line = get_tleap_neutralization_line(ion_concentration, 
                                                     n_waters, sys_CHARGE)
    
    leap_neutral = \
    f"""
    {_head_lines}

    system =  loadpdb {tmp_dir}/{input_basename}.pdb
    {_neutralize_line}
    setBox system vdw
    check system
    charge system

    savepdb system {tmp_dir}/{output_basename}.pdb
    saveamberparm system {tmp_dir}/{output_basename}.prmtop {tmp_dir}/{output_basename}.rst7
    quit
    """

    with open(f'{tmp_dir}/leap_prep_{output_basename}.in', 'w') as f:
        f.write(leap_neutral)

    # Run tleap
    tleap_command_neutral = f'tleap -f {tmp_dir}/leap_prep_{output_basename}.in'
    ps = subprocess.Popen(tleap_command_neutral, 
                          shell  = True, 
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT,
                          encoding = 'UTF-8')
    output, error = ps.communicate()

    # Move temp files
    os.rename('leap.log', f'{tmp_dir}/leap_prep_{output_basename}.in.log')

    if verbose:
        print(output)
        print(error)



def run_center_system_to_origin(filename: str, 
                                topology: str,
                                out_filename: str = None,
                                overwrite: bool = True):
    import pytraj
    system = pytraj.load(
                filename = filename, 
                top = topology)
    x, y, z = pytraj.center_of_geometry(traj = system)[0]

    print(f'Moving system to origin ({x:.2f}, {y:.2f}, {z:.2f}) -> (0,0,0)')
    system = pytraj.center(traj = system, 
                           center = 'origin')
    pytraj.write_traj(
                filename = out_filename, 
                traj = system, 
                overwrite = overwrite)
    if overwrite:
        pass



def create_posre_prot(input_pdb: str, 
                      posre_basename: str,
                      out_dir: str,
                      verbose: bool = True):
    
    POSRE_FILENAME = f"{posre_basename}.prot_posre.itp"
    
    gmx_posre_command = f'gmx pdb2gmx -f {input_pdb}' +\
                        f' -o ./TEMP.gro' +\
                        f' -p ./TEMP.top' +\
                        f' -i {out_dir}/{POSRE_FILENAME}' +\
                         ' -ff amber99sb-ildn -water tip3p -ignh'
    
    ps = subprocess.Popen(gmx_posre_command, 
                           shell  = True, 
                           stdout = subprocess.PIPE,
                           stderr = subprocess.STDOUT,
                           encoding = 'UTF-8')
    output, error = ps.communicate()
    os.remove('./TEMP.gro')
    os.remove('./TEMP.top')
    if verbose:
        print(output)
        print(error)



def include_posre_in_top(input_top_file: str, 
                         out_top_file: str, 
                         posre_filename: str,
                         molecule_index: int = 1,
                         posre_sufix: str = 'prot_posre'
                       ):
    POSRES_LINES = [
        '; Include Position restraint file for protein\n',
        '#ifdef POSRES\n',
        f'#include "{posre_filename}"\n',
        '#endif\n',
        '\n','\n'
    ]

    with open(input_top_file) as in_topfile:
        top_lines = in_topfile.readlines()
        mol_type_str = '[ moleculetype ]'
        mol_type_idxs = [i for (i, l) in enumerate(top_lines) 
                         if l.startswith(mol_type_str)]

    # Insert POSRES_LINES between the last lines of the protein definition
    # and the lines of the solvent definitions
    second_molecule_type_idx = mol_type_idxs[molecule_index]
    FIRST_BLOCK  = top_lines[:second_molecule_type_idx]
    SECOND_BLOCK = top_lines[second_molecule_type_idx:]

    with open(out_top_file, 'w') as out_topfile:
        out_topfile.writelines(
            FIRST_BLOCK + POSRES_LINES + SECOND_BLOCK
        )


def create_posre_file(input_gro: str, 
                      out_posre_filename: str,
                      tmp_dir: str = '.',
                      mol_selection = 'protein',
                      rest_energies = [1000, 1000, 1000],
                      specify_index_mol = None,
                      verbose: bool = True):
    
    SELECTIONS = {
        'protein': "2\n q",
        'ligand': "r LIG & ! a H*\n q",
        'MG': "r MG\n q"
    }
    
    # assert mol_selection in SELECTIONS.keys(), \
    # 'Only `protein` and `ligand` sel available for `mol_selection`'

    # Create the .gro file
    if specify_index_mol:
        GRO_SELECTION = str(specify_index_mol)
    elif mol_selection == 'ligand':
        GRO_SELECTION = '13'
    elif mol_selection == 'protein':
        GRO_SELECTION = '1'
    
    # Create the gro file
    gmx_gro_command = f'echo "{GRO_SELECTION}" | gmx trjconv' +\
                      f' -f {input_gro} -s {input_gro}' +\
                      f' -o {out_posre_filename}.{mol_selection}.gro'
    ps_gro = subprocess.Popen(gmx_gro_command, 
                           shell  = True, 
                           stdout = subprocess.PIPE,
                           stderr = subprocess.STDOUT,
                           encoding = 'UTF-8')
    output, error = ps_gro.communicate()

    
    # make_ndx
    gmx_make_ndx = f'echo "' + SELECTIONS[mol_selection] + '" |' +\
                        f' gmx make_ndx' +\
                        f' -f {out_posre_filename}.{mol_selection}.gro' +\
                        f' -o {out_posre_filename}.ndx'
    ps_ndx = subprocess.Popen(gmx_make_ndx, 
                           shell  = True, 
                           stdout = subprocess.PIPE,
                           stderr = subprocess.STDOUT,
                           encoding = 'UTF-8')
    output, error = ps_ndx.communicate()
    
    
    SELECTIONS_GPRS = {
        'protein': "2",
        'ligand': "LIG_&_!H*",
        'MG': "3"
    }
    
    gmx_genrestr = f'echo "' + SELECTIONS_GPRS[mol_selection] + '" |' +\
                        f' gmx genrestr' +\
                        f' -f {input_gro}' +\
                        f' -n {out_posre_filename}.ndx' +\
                        f' -o {out_posre_filename}' +\
                        f' -fc ' +\
                        ' '.join([str(e) for e in rest_energies])
    ps_rest = subprocess.Popen(gmx_genrestr, 
                           shell  = True, 
                           stdout = subprocess.PIPE,
                           stderr = subprocess.STDOUT,
                           encoding = 'UTF-8')
    output, error = ps_rest.communicate()
    
    os.remove(f'{out_posre_filename}.ndx')
    os.remove(f'{out_posre_filename}.{mol_selection}.gro')
    
    if verbose:
        print(output)
        print(error)