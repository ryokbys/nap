"""
Module for handling in.pmd file and related information.

Usage:
  inpmd.py [options] INPMD_NAME

Options:
  -h, --help  Show this help message and exit.
"""
import copy
from docopt import docopt

_entry_to_varname = {
    'num_nodes_x': 'nx',
    'num_nodes_y': 'ny',
    'num_nodes_z': 'nz',
    'io_format': 'ciofmt',
    'print_level': 'iprint',
    'time_interval': 'dt',
    'num_iteration': 'nstp',
    'min_iteration': 'minstp',
    'num_out_energy': 'nerg',
    'flag_out_pmd': 'ifpmd',
    'flag_out_pos': 'ifpmd',
    'combine_out_pos': 'lcomb_pos',
    'combine_out_pmd': 'lcomb_pos',
    'num_out_pmd': 'npmd',
    'num_out_pos': 'npmd',
    'flag_sort': 'ifsort',
    'force_type': 'force_list',
    'cutoff_radius': 'rc',
    'cutoff_buffer': 'rbuf',
    'flag_damping': 'ifdmp',
    'damping_coeff': 'dmp',
    'converge_eps': 'eps_conv',
    'converge_num': 'n_conv',
    'initial_temperature': 'tinit',
    'final_temperature': 'tfin',
    'flag_multi_temp': 'lmultemps',
    'temperature_control': 'ctctl',
    'temperature_target': 'ttgt',
    'temperature_relax_time': 'trlx',
    'temperature_limit': 'tlimit',
    'flag_temp_dist': 'ltdst',
    'factor_direction': 'fmv',
    'stress_control': 'cpctl',
    'pressure_target': 'ptgt',
    'stress_target': 'stgt',
    'stress_relax_time': 'srlx',
    'pressure_relax_time': 'srlx',
    'flag_compute_stress': 'lstrs0',
    'cell_fix': 'lcellfix',
    'zload_type': 'czload_type',
    'final_strain': 'strfin',
    'boundary': 'boundary',
    'max_num_neighbors': 'nnmax',
    'allow_reallocation': 'lrealloc',
    'remove_translation': 'nrmtrans',
}

_default_params = {
    'num_nodes_x': 1,
    'num_nodes_y': 1,
    'num_nodes_z': 1,
    'io_format': 'ascii',
    'print_level': 1,
    'time_interval': 1.0,
    'num_iteration': 0,
    'min_iteration': 0,
    'num_out_energy': 10,
    'flag_out_pmd': 3,
    'num_out_pmd': 1,
    'combine_out_pos': True,
    'flag_sort': 1,
    'force_type': None,
    'cutoff_radius': 5.0,
    'cutoff_buffer': 0.0,
    'flag_damping': 0,
    'damping_coeff': 0.95,
    'converge_eps': 1e-4,
    'converge_num': 3,
    'initial_temperature': -10.0,
    'final_temperature': -10.0,
    'temperature_control': 'none',
    'flag_multi_temp': False,
    'temperature_target': [300.0, 100.0, 300.0,
                           300.0, 300.0, 300.0,
                           300.0, 300.0, 300.0],
    'temperature_relax_time': 100.0,
    'temperature_limit': 1.0e+5,
    'remove_translation': 1,
    'flag_temp_dist': False,
    'factor_direction':[[1.0, 1.0, 1.0],
                        [1.0, 0.0, 1.0],
                        [1.0, 1.0, 1.0],
                        [1.0, 0.0, 1.0],
                        [1.0, 1.0, 1.0],
                        [1.0, 0.0, 1.0],
                        [1.0, 1.0, 1.0],
                        [1.0, 0.0, 1.0],
                        [1.0, 1.0, 1.0]],
    'stress_control': 'none',
    'pressure_target': 0.0,
    'stress_target': [[0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0]],
    'stress_relax_time': 20.0,
    'cell_fix': [[False,False,False],
                 [False,False,False],
                 [False,False,False]],
    'flag_compute_stress': True,
    'zload_type': 'none',
    'final_strain': 0.0,
    'boundary': 'ppp',
    'max_num_neighbors': 70,
    'allow_reallocation': True,
}

_int_keys = [
    'num_nodes_x', 'num_nodes_y', 'num_nodes_z',
    'num_iteration', 'num_out_energy', 'flag_out_pmd',
    'num_out_pmd', 'flag_damping',
    'converge_num', 'min_iteration', 'flag_sort',
    'print_level', 'max_num_neighbors', 'remove_translation'
]
_float_keys = [
    'time_interval', 'cutoff_radius', 'cutoff_buffer',
    'damping_coeff', 'initial_temperature',
    'final_temperature', 'temperature_limit',
    'temperature_relax_time', 'pressure_target',
    'stress_relax_time', 'shear_stress',
    'converge_eps', 'final_strain'
]
_str_keys = [
    'io_format', 'force_type', 'temperature_control',
    'stress_control',
    'zload_type', 'boundary',
    'overlay_type',
]
_bool_keys = [
    'allow_reallocation', 'flag_temp_dist', 'flag_compute_stress',
    'flag_multi_temp', 'combine_out_pos', 'combine_out_pmd'
]

def get_default():
    return _default_params

def correct_fortran_double(arg):
    """
    Correct fortran-type double precision variable, 1.0d-3,
    to python float, 1.0e-3.
    """

    if type(arg) != str:
        raise ValueError('arg is not a string.')
    if arg[0].isdigit() and 'd' in arg:
        arg = arg.replace('d','e')
    return arg

def read_inpmd(fname='in.pmd'):
    """
    Read in.pmd and return parameters as a dictionary.
    """
    with open(fname,'r') as f:
        lines = f.readlines()
    inputs = copy.copy(_default_params)
    mode = None
    for line in lines:
        if line[0] in ('#','!'):
            continue
        data = line.split()
        data = [ correct_fortran_double(d) for d in data ]
        if len(data) == 0:
            mode = None
            continue
        key = data[0]
        if key in _entry_to_varname.keys():
            mode = None
        if mode is not None:
            if mode == 'factor_direction':
                inputs[mode][facdir_inc] = [ float(x) for x in data ]
                facdir_inc += 1
            elif mode == 'stress_target':
                inputs[mode].append([ float(x) for x in data ])
            elif mode == 'cell_fix':
                inputs[mode].append([ x in ('T','True','.true.') for x in data ])
            else:
                continue
        else:
            #...Special entries
            if key == 'force_type':
                inputs[key] = data[1:]
                mode = None
            elif key == 'temperatuer_target':
                if inputs['flag_multi_temp']:
                    ifmv = int(data[1])
                    if not 1 <= ifmv <= 9:
                        ValueError('temperature_target ifmv wrong.')
                    inputs[key][ifmv-1] = float(data[2])
                else:
                    inputs[key][0] = float(data[1])
                mode = None
            elif key == 'factor_direction':
                mode = key
                facdir_dim = int(data[1])
                facdir_num = int(data[2])
                facdir_inc = 0
                # inputs[key] = []
                continue
            elif key == 'stress_target':
                mode = key
                inputs[key] = []
                continue
            #...Non special entries
            elif key in _int_keys:
                inputs[key] = int(data[1])
                mode = None
            elif key in _float_keys:
                inputs[key] = float(data[1])
                mode = None
            elif key in _str_keys:
                inputs[key] = data[1]
                mode = None
            elif key in _bool_keys:
                inputs[key] = data[1] in ('T', 'True', '.true.')
                mode = None

    return inputs

def inputs_to_vars(inputs={}):
    """
    Convert inputs, which is a dictionary with keys of entries,
    to vars, whose names are used in pmd source code.
    """
    vs = {}
    for k,v in inputs.items():
        if k not in _entry_to_varname.keys():
            continue
        vname = _entry_to_varname[k]
        vs[vname] = v
    return vs

def main():
    import os,sys
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])))

    print("Under construction...")
    return None

if __name__ == '__main__':
    main()
