"""
ASE interface to pmd.
"""

__author__  = "Ryo KOBAYASHI"
__version__ = "160605"
__LICENSE__ = "MIT"

import os
import copy
import subprocess
import numpy as np
import ase.calculators.calculator
from ase.calculators.calculator import FileIOCalculator,Calculator

from pmdio import read_pmd,write_pmd,get_fmvs

CALC_END_MARK = "time  "
some_changes = ['positions', 'numbers', 'cell',]

class PMD(FileIOCalculator):
    """
    Class for doint PMD/SMD calculations.

    calc = PMD(label='pmd')
    """

    implemented_properties = ['energy', 'forces',
                              'relaxed_scaled_positions',
                              'stress']
    command = 'pmd'

    default_parameters = {
        'num_nodes_x': 1,
        'num_nodes_y': 1,
        'num_nodes_z': 1,
        'io_format': 'ascii',
        'time_interval': 1.0,
        'num_iteration': 0,
        'min_iteration': 0,
        'num_out_energy': 10,
        'flag_out_pmd': 1,
        'num_out_pmd':  1,
        'flag_sort':  1,
        'force_type': None,
        'cutoff_radius': 5.0,
        'cutoff_buffer': 0.0,
        'flag_damping': 0,
        'damping_coeff': 0.95,
        'converge_eps': 1e-4,
        'converge_num': 3,
        'initial_temperature': -10.0,
        'final_temperature': -10.0,
        'temperature_control': 'None',
        'temperature_target': [300.0, 100.0,],
        'temperature_relax_time': 100.0,
        'factor_direction':[[1.0, 1.0, 1.0],
                            [1.0, 0.0, 1.0]],
        'flag_isobaric': 0,
        'pressure_target': 0.0,
        'vol_mass_coeff': 1e-3,
        'vol_change_damping': 1.0,
        'shear_stress': 0.0,
        'mass': [28.0855,4.0,]
    }

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='pmd', atoms=None, 
                 command='pmd > out.pmd', 
                 specorder=None, **kwargs):
        """Construct PMD-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'pmd'.

        Examples
        ========
        Use default values:

        >>> h = Atoms('H', calculator=PMD(label='pmd',force_type='NN'))
        >>> e = h.get_potential_energy()

        """

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

        if not label in ['pmd']:
            raise RuntimeError('label must be pmd.')

        if self.parameters['force_type'] is None:
            raise RuntimeError('force_type must be specified.')

        if command is None:
            self.command = self.label+' > out.'+self.label
        elif '>' in command:
            self.command = command.split('>')[0] +' > out.'+self.label
        else:
            self.command = command +' > out.'+self.label

        self.specorder= specorder

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=some_changes):
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
        self.read_results()


    def relax(self, atoms=None, properties=['energy'],
              system_changes=some_changes,
              flag_damping=2,
              converge_eps=1.0e-3,nsteps=100,min_iteration=5,
              converge_num=3,time_interval=2.0,
              initial_temperature=10.0):
        """
        Relax atom positions by running damped MD in pmd instead of using
        optimize module in ASE.
        """
        self.set(flag_damping=flag_damping,
                 converge_eps=converge_eps,
                 num_iteration=nsteps,
                 min_iteration=min_iteration,
                 converge_num=converge_num,
                 time_interval=time_interval,
                 initial_temperature=initial_temperature,
                 flag_sort=1)
        Calculator.calculate(self, atoms, properties, system_changes)
        # print 'type(atoms),type(self.atoms)= ',type(atoms),type(self.atoms)
        # print 'self.atoms = ',self.atoms
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
        """Write input parameters to files-file."""

        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        if self.label == 'pmd':
            infname = 'in.pmd'
            #self.write_pmd(atoms)
            write_pmd(atoms,fname='pmd0000',specorder=self.specorder)
            
        with open(infname,'w') as f:
            fmvs,ifmvs = get_fmvs(atoms)
            f.write(get_input_txt(self.parameters,fmvs))
        
        
    def read_results(self,relax=False):
        """
        Only erg.pmd and frc.pmd are to be read.
        """
        import glob
        outfname= 'out.'+self.label
        ergfname= 'erg.'+self.label
        frcfname= 'frc.'+self.label
        strfname= 'str.'+self.label
        if not os.path.exists(outfname):
            raise RuntimeError(outfname+' does not exists.')
        if not os.path.exists(ergfname):
            raise RuntimeError(ergfname+' does not exists.')
        if not os.path.exists(frcfname):
            raise RuntimeError(frcfname+' does not exists.')
        if not os.path.exists(strfname):
            print 'Warning: '+frcfname+' does not exists.'

        fout= open(outfname,'r')
        lines= fout.readlines()
        if not CALC_END_MARK in  lines[-1]:
            raise RuntimeError(self.label+' seems to stop somewhere..')
        if relax:
            relax_converged = False
            for line in lines:
                if 'damped MD converged with' in line:
                    relax_converged = True
            if not relax_converged:
                print ''
                print '** Warning: pmd relaxation does not seem to be converged**'
                print ''
        fout.close()

        self.results={}
        
        with open(ergfname,'r') as f:
            erg = float(f.readline().split()[0])
            self.results['energy'] = erg

        with open(frcfname,'r') as f:
            num= int(f.readline().split()[0])
            frcs= np.zeros((num,3))
            for i in range(num):
                data= [ float(x) for x in f.readline().split() ]
                frcs[i,:] = data[:]
            self.results['forces'] = frcs

        if os.path.exists(strfname):
            with open(strfname,'r') as f:
                strs = np.array([ float(x) for x in f.readline().split() ])
            self.results['stress'] = strs
        
        if relax:
            posfile = max(glob.glob('pmd[0-9]*'), key=os.path.getctime)
            tmpatoms = read_pmd(fname=posfile,specorder=self.specorder)
            self.results['relaxed_scaled_positions'] = tmpatoms.get_scaled_positions()

    def get_relaxed_scaled_positions(self):
        # print 'self.results = ',self.results
        # relaxed_spos = self.get_property('relaxed_scaled_positions',self.atoms)
        return self.results['relaxed_scaled_positions']
        
def get_input_txt(params,fmvs):
    txt = ''

    order=['num_nodes_x','num_nodes_y','num_nodes_z','',
           'io_format','',
           'time_interval','num_iteration','min_iteration','num_out_energy','',
           'flag_out_pmd','num_out_pmd','flag_sort','',
           'force_type','cutoff_radius','cutoff_buffer','',
           'flag_damping','damping_coeff','converge_eps','converge_num','',
           'initial_temperature','final_temperature',
           'temperature_control','temperature_target',
           'temperature_relax_time','',
           'factor_direction','',
           'flag_isobaric','pressure_target',
           'vol_mass_coeff','vol_change_damping','shear_stress','',
           'mass','',]

    int_keys=['num_nodes_x','num_nodes_y','num_nodes_z',
              'num_iteration','num_out_energy','flag_out_pmd',
              'num_out_pmd','flag_damping','flag_isobaric',
              'converge_num','min_iteration','flag_sort']
    float_keys=['time_interval','cutoff_radius','cutoff_buffer',
                'damping_coeff','initial_temperature',
                'final_temperature',
                'temperature_relax_time','pressure_target',
                'vol_mass_coeff','vol_change_damping','shear_stress',
                'converge_eps']
    str_keys=['io_format','force_type','temperature_control']

    for key in order:
        
        # special keys first
        if key is '':
            txt += '\n'
        elif key is 'temperature_target':
            vals = params[key]
            for i,v in enumerate(vals):
                txt += '{0:25s} {1:2d} {2:6.1f}\n'.format(key,i+1,v)
        elif key is 'factor_direction':
            #vals = params[key]
            vals= fmvs
            txt += '{0:25s} 3 {1:d}\n'.format(key,len(vals))
            for i,v in enumerate(vals):
                txt += '  {0:6.2f} {1:6.2f} {2:6.2f}\n'.format(v[0],v[1],v[2])
        elif key is 'mass':
            vals = params[key]
            for i,v in enumerate(vals):
                txt += '{0:25s} {1:2d} {2:10.4f}\n'.format(key,i+1,v)
        elif key is 'converge_eps':
            txt += '{0:25s} {1:10.1e}\n'.format(key,params[key])
        elif key in int_keys:
            txt += '{0:25s} {1:3d}\n'.format(key,params[key])
        elif key in float_keys:
            txt += '{0:25s} {1:6.1f}\n'.format(key,params[key])
        elif key in str_keys:
            txt += '{0:25s} {1:s}\n'.format(key,params[key])
        else:
            raise RuntimeError('Input parameter '+key+' is not defined.')

    return txt
