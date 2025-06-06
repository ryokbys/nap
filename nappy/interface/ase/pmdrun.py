"""
ASE interface to pmd.

Changes:
  <2024.05.04 RK>  Changed to call pmd via f2py and dynamic library instead of via subprocess.
"""
import os
import subprocess
import numpy as np
from ase.calculators.calculator import FileIOCalculator,Calculator

import nappy
from nappy.interface.ase.pmdio import get_fmvs
from nappy.napsys import NAPSystem

try:
    import nappy.pmd.pmd_wrapper as pw
except:
    pass

__author__  = "Ryo KOBAYASHI"
__version__ = "240504"
__LICENSE__ = "MIT"

CALC_END_MARK = "Job finished "
some_changes = ['positions', 'numbers', 'cell',]

def str2char(string,nlen):
    slen = len(string)
    if slen > nlen:
        raise ValueError('slen > nlen !')
    c = np.empty(nlen,dtype='c')
    c[:] = ' '
    c[0:slen] = [ string[i] for i in range(slen) ]
    c[slen:] = [ ' ' for i in range(nlen-slen) ]
    return c


class PMD_ASE_Calculator(Calculator):
    """
    Class for PMD calculation for ASE.

    calc = PMD(label='pmd')
    """

    implemented_properties = ['energy', 'forces',
                              'relaxed_scaled_positions',
                              'stress',
                              'num_step_relax',
                              'relaxed_cell']
    command = 'pmd'

    default_parameters = {
        'num_nodes_x': -1,
        'num_nodes_y': -1,
        'num_nodes_z': -1,
        'num_omp_threads': 4,
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
        # 'mass': {'W':183.14, 'H':1.008,},
        'zload_type': 'none',
        'final_strain': 0.0,
        'boundary': 'ppp',
    }

    def __init__(self, restart=None,
                 label='pmd', atoms=None, 
                 command='pmd > out.pmd',
                 dimension=(True,True,True),
                 specorder=None, **kwargs):
        """Construct PMD-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (in.label, erg.label, ...). [Default: 'pmd']
        command: str
            Command string to execute pmd. [Default: 'pmd > out.pmd']
        force_type: str or tuple/list
            Force fields to be used, which must be set. [Default: None]
        specorder: list
            Order of species. This is probably very important, since the order of species
            is usually fixed in pmd whereas not in ASE atoms object. [Default: None]

        Examples
        ========
        Use default values:

        >>> h = Atoms('H', calculator=PMD(label='pmd',force_type='NN'))
        >>> e = h.get_potential_energy()

        """

        Calculator.__init__(self, restart=restart,
                            label=label, atoms=atoms,
                            directory='.', **kwargs)

        if label != 'pmd':
            raise RuntimeError('label must be pmd.')

        if self.parameters['force_type'] is None:
            raise RuntimeError('force_type must be specified.')

        if command is None:
            self.command = self.label+' > out.'+self.label
        elif '>' in command:
            self.command = command.split('>')[0] +' > out.'+self.label
        else:
            self.command = command +' > out.'+self.label

        self.dimension = dimension

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=some_changes):
        import nappy.pmd
        Calculator.calculate(self, atoms, properties, system_changes)
        # # Write in.pmd to be read from pmd executable
        # self.write_input(self.atoms, properties, system_changes)

        # olddir = os.getcwd()
        # try:
        #     os.chdir(self.directory)
        #     errorcode = subprocess.call(self.command, shell=True)
        # finally:
        #     os.chdir(olddir)

        # if errorcode:
        #     raise RuntimeError('%s returned an error: %d' %
        #                        (self.name, errorcode))
        # self.read_results()
        
        #...The code below is for calling pmd via f2py and dynamic library
        if atoms is not None:
            self.nsys = nappy.io.from_ase(atoms, get_forces=False)

        #...Call pmd via f2py and dynamic library
        pmd = nappy.pmd.PMD(self.nsys)
        pmd.set_params(force_type = self.parameters['force_type'],
                       cutoff_radius = self.parameters['cutoff_radius'],)
        pmd.run(nstp=0, iprint=0)

        self.results['energy'] = pmd.get_potential_energy()
        self.results['forces'] = pmd.get_forces().T
        self.results['stresses'] = pmd.get_stress()
        self.results['stress'] = pmd.get_stress().sum() /3
        return None

    def relax(self, atoms=None, properties=['energy'],
              system_changes=some_changes,
              flag_damping=1,damping_coeff=0.95,
              converge_eps=1.0e-3,num_iteration=100,min_iteration=5,
              converge_num=3,time_interval=0.5,
              initial_temperature=10.0, stress_control='none',
              pressure_target=0.0,
              stress_target=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]):
        """
        Relax atom positions by running damped MD in pmd instead of using
        optimize module in ASE.
        """
        self.set(flag_damping=flag_damping,
                 damping_coeff=damping_coeff,
                 converge_eps=converge_eps,
                 num_iteration=num_iteration,
                 min_iteration=min_iteration,
                 num_out_energy=num_iteration,
                 converge_num=converge_num,
                 time_interval=time_interval,
                 initial_temperature=initial_temperature,
                 stress_control=stress_control,
                 stress_target=stress_target,
                 flag_sort=1)
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)

        olddir = os.getcwd()
        try:
            os.chdir(self.directory)
            errorcode = subprocess.call(self.command, shell=True)
        finally:
            os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))
        self.read_results(relax=True)

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input parameters to in.pmd file."""

        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        if self.label == 'pmd':
            infname = 'in.pmd'
            # write_pmd(atoms,fname='pmdini',specorder=self.specorder)
            #nsys = NAPSystem.from_ase_atoms(atoms,specorder=self.specorder)
            nsys = nappy.io.from_ase(atoms,specorder=self.specorder)
            nappy.io.write_pmd(nsys,fname='pmdini')
            
        with open(infname,'w') as f:
            fmvs,ifmvs = get_fmvs(atoms)
            for i in range(len(fmvs)):
                for ii in range(3):
                    if fmvs[i][ii] > 0.1 and not self.dimension[ii]:
                        fmvs[i][ii] = 0.0
            f.write(get_input_txt(self.parameters,fmvs))
        
    def read_results(self,relax=False):
        """
        Only erg.pmd and frc.pmd are to be read.
        """
        outfname= 'out.'+self.label
        ergfname= 'erg.'+self.label
        frcfname= 'frc.'+self.label
        strfname= 'strs.'+self.label
        if not os.path.exists(outfname):
            raise RuntimeError(outfname+' does not exists.')
        if not os.path.exists(ergfname):
            raise RuntimeError(ergfname+' does not exists.')
        if not os.path.exists(frcfname):
            raise RuntimeError(frcfname+' does not exists.')
        if not os.path.exists(strfname):
            print('Warning: '+strfname+' does not exists.')

        self.results={ k : None for k in self.implemented_properties}

        fout= open(outfname,'r')
        lines= fout.readlines()
        if CALC_END_MARK not in lines[-1]:
            raise RuntimeError(self.label+' seems to stop somewhere..')
        if relax:
            relax_converged = False
            num_step_relax = -1
            for line in lines:
                if 'Damped MD converged with' in line:
                    relax_converged = True
                    num_step_relax = int(line.split()[4])
                    break
            if not relax_converged:
                print('')
                print('** Warning: pmd relaxation does not' +\
                    ' seem to be converged**')
                print('')
            self.results['num_step_relax'] = num_step_relax
        fout.close()

        with open(ergfname,'r') as f:
            erg = float(f.readline().split()[0])
            self.results['energy'] = erg

        with open(frcfname,'r') as f:
            num= int(f.readline().split()[0])
            frcs= np.zeros((num,3))
            for i in range(num):
                data= [ float(x) for x in f.readline().split() ]
                frcs[i,0:3] = data[0:3]
            self.results['forces'] = frcs

        if os.path.exists(strfname):
            try:
                with open(strfname,'r') as f:
                    strs = np.array([ float(x) for x in f.readline().split() ])
                self.results['stress'] = strs
            except:
                self.results['srress'] = None
        
        if relax:
            posfile = 'pmdfin'
            #nsys = NAPSystem(fname=posfile,specorder=self.specorder)
            nsys = nappy.io.read(fname=posfile,specorder=self.specorder)
            #tmpatoms = read_pmd(fname=posfile,specorder=self.specorder)
            tmpatoms = nsys.to_ase_atoms()
            self.results['relaxed_scaled_positions'] \
                = tmpatoms.get_scaled_positions()
            self.results['relaxed_cell'] = tmpatoms.get_cell()

    def get_relaxed_scaled_positions(self):
        return self.results['relaxed_scaled_positions']

    def get_relaxed_cell(self):
        return self.results['relaxed_cell']

def get_input_txt(params,fmvs):
    txt = ''

    order=['num_nodes_x','num_nodes_y','num_nodes_z','num_omp_threads','',
           'io_format','print_level','',
           'time_interval','num_iteration','min_iteration','num_out_energy','',
           'flag_out_pmd','num_out_pmd','flag_sort','',
           'force_type','cutoff_radius','cutoff_buffer','',
           'flag_damping','damping_coeff','converge_eps','converge_num','',
           'initial_temperature','final_temperature',
           'temperature_control','temperature_target',
           'temperature_relax_time','flag_temp_dist','',
           'factor_direction','',
           'stress_control','pressure_target','stress_target',
           'stress_relax_time','',
           'overlay','overlay_type','',
           'final_strain','',
           'boundary']

    int_keys=['num_nodes_x','num_nodes_y','num_nodes_z','num_omp_threads',
              'num_iteration','num_out_energy','flag_out_pmd',
              'num_out_pmd','flag_damping','print_level',
              'converge_num','min_iteration','flag_sort']
    float_keys=['time_interval','cutoff_radius','cutoff_buffer',
                'damping_coeff','initial_temperature',
                'final_temperature',
                'temperature_relax_time','pressure_target',
                'stress_relax_time','shear_stress',
                'converge_eps','final_strain']
    str_keys=['io_format','force_type','temperature_control',
              'stress_control','flag_temp_dist',
              'boundary',
              'overlay_type']

    for key in order:
        
        # special keys first
        if key == '':
            txt += '\n'
        elif key not in params:
            continue
        elif key == 'force_type':
            vals = params[key]
            txt += '{0:25s} '.format('force_type')
            if isinstance(vals,str):
                txt += '  {0:s}'.format(vals)
            elif type(vals) in (list,tuple):
                for v in vals:
                    txt += '  {0:s}'.format(v)
            txt += '\n'
        elif key == 'temperature_target':
            vals = params[key]
            for i,v in enumerate(vals):
                #...Format of temperature target has been changed...
                # txt += '{0:25s} {1:2d} {2:6.1f}\n'.format(key,i+1,v)
                txt += f'{key:25s} {v:6.1f}\n'
        elif key == 'factor_direction':
            # vals = params[key]
            vals= fmvs
            txt += '{0:25s} 3 {1:d}\n'.format(key,len(vals))
            for i,v in enumerate(vals):
                txt += '  {0:6.2f} {1:6.2f} {2:6.2f}\n'.format(v[0],v[1],v[2])
        elif key == 'stress_target':
            vals = params[key]
            txt += '{0:25s}\n'.format(key)
            for i in range(3):
                v = vals[i]
                txt += '  {0:6.2f} {1:6.2f} {2:6.2f}\n'.format(v[0],v[1],v[2])
        # elif key == 'mass':
        #     masses = params[key]
        #     for k,v in masses.items():
        #         txt += '{0:25s} {1:3s}  {2:10.4f}\n'.format(key,k,v)
        elif key == 'converge_eps':
            txt += '{0:25s} {1:10.1e}\n'.format(key,params[key])
        elif key == 'overlay':
            ols = params[key]
            for k,v in ols.items():
                txt += '{0:25s} {1:3s}  {2:6.2f} {3:6.2f}\n'.format(key,k,v[0],v[1])
        elif key in int_keys:
            txt += '{0:25s} {1:3d}\n'.format(key,params[key])
        elif key in float_keys:
            if key == 'cutoff_radius' or key == 'damping_coeff':
                txt += '{0:25s} {1:8.4f}\n'.format(key,params[key])
            else:
                txt += '{0:25s} {1:6.1f}\n'.format(key,params[key])
        elif key in str_keys:
            txt += '{0:25s} {1:s}\n'.format(key,params[key])
        else:
            raise RuntimeError('Input parameter '+key+' is not defined.')

    return txt
