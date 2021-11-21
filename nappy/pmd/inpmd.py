"""
Module for handling in.pmd file and related information.
"""
from __future__ import print_function

import copy

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
    'temperature_control': 'ctctl',
    'temperature_target': 'ttgt',
    'temperature_relax_time': 'trlx',
    'flag_temp_dist': 'ltdst',
    'factor_direction': 'fmv',
    'stress_control': 'cpctl',
    'pressure_target': 'ptgt',
    'stress_target': 'stgt',
    'stress_relax_time': 'srlx',
    'pressure_relax_time': 'srlx',
    'flag_compute_stress': 'lstrs0',
    'zload_type': 'czload_type',
    'final_strain': 'strfin',
    'boundary': 'boundary',
    'max_num_neighbors': 'nnmax',
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
    'flag_out_pmd': 2,
    'num_out_pmd': 1,
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
    'temperature_target': [300.0, 100.0,],
    'temperature_relax_time': 100.0,
    'flag_temp_dist': 'F',
    'factor_direction':[[1.0, 1.0, 1.0],
                        [1.0, 0.0, 1.0]],
    'stress_control': 'none',
    'pressure_target': 0.0,
    'stress_target': [[0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0]],
    'stress_relax_time': 20.0,
    'flag_compute_stress': 'T',
    'zload_type': 'none',
    'final_strain': 0.0,
    'boundary': 'ppp',
    'max_num_neighbors': 200,
}

_int_keys = [
    'num_nodes_x','num_nodes_y','num_nodes_z',
    'num_iteration','num_out_energy','flag_out_pmd',
    'num_out_pmd','flag_damping',
    'converge_num','min_iteration','flag_sort',
    'print_level','max_num_neighbors',
]
_float_keys = [
    'time_interval','cutoff_radius','cutoff_buffer',
    'damping_coeff','initial_temperature',
    'final_temperature',
    'temperature_relax_time','pressure_target',
    'stress_relax_time','shear_stress',
    'converge_eps','final_strain'
]
_str_keys = [
    'io_format','force_type','temperature_control',
    'stress_control','flag_temp_dist',
    'flag_compute_stress','zload_type','boundary',
    'overlay_type'
]

def get_default():
    return _default_params

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
        if len(data) == 0:
            mode = None
            continue
        key = data[0]
        if key in _entry_to_varname.keys():
            mode = None
        if mode is not None:
            if mode == 'factor_direction':
                inputs[mode].append([ float(x) for x in data ])
            elif mode == 'stress_target':
                inputs[mode].append([ float(x) for x in data ])
            else:
                continue
        else:
            #...Special entries
            if key == 'force_type':
                inputs[key] = data[1:]
                mode = None
            elif key == 'temperatuer_target':
                if type(inputs[key]) is list:
                    pass
                else:
                    inputs[key] = []
                inputs[key].append(float(data[2]))
                mode = None
            elif key == 'factor_direction':
                mode = key
                inputs[key] = []
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
