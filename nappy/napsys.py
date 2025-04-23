#!/usr/bin/env python
"""
System information used in pmd, including cell information, lattice constant,
number of atoms, atom species, atom positions.

Users can use `napsys.py` as a converter program as described below.
Available formats are,
  pmd, POSCAR, dump, xsf

Usage:
  napsys.py convert [options] INFILE OUTFILE
  napsys.py analyze [options] INFILE

Options:
  -h, --help  Show this help message and exit.
  --in-format=INFMT
              Format of the input file. [default: None]
  --out-format=OUTFMT
              Format of the output file. [default: None]
  --specorder=SPECORDER
              Order of species. [default: None]
  --scale=SCALE
              Scale the cell. [default: None]
  --shift=SHIFT
              Shift original atom positions within a periodic system. [default: 0.0,0.0,0.0]
  --periodic-copy=COPIES
              Number of copies of a,b,and c directions. [default: 1,1,1]
  --cycle-coord=NCYCLE
              Number of cyclic shift of coordinates, e.g. x-y-z to y-z-x. [default: 0]
  --charges=CHARGES
              Charges of species. Charges should be numbers separated by comma.
              If it is None, no charge is set.
              If it is set and the output format accepts, charges are written.
              [default: None]
"""
import copy
from docopt import docopt
import numpy as np

import nappy
from nappy.atom import get_symbol_from_number, get_number_from_symbol
from nappy.util import pbc

from icecream import ic

__author__ = "RYO KOBAYASHI"
__version__ = "240323"

#...constants
# FILE_FORMATS = ('pmd','POSCAR','dump','xsf','lammps',
#                 'cube','CHGCAR','pdb')
DEFAULT_LABELS = ('x','y','z','vx','vy','vz','fx','fy','fz','sid')
DEFAULT_DTYPES = {'x':'float64', 'y':'float64', 'z':'float64',
                  'vx':'float32', 'vy':'float32', 'vz':'float32',
                  'fx':'float32', 'fy':'float32', 'fz':'float32',
                  'sid':'int8', 'chg':'float16', 'ekin':'float32',
                  'epot':'float32'}
# _file_formats = ('pmd','POSCAR','dump','xsf','lammps',
#                  'cube','CHGCAR')
# _default_labels = ('pos','vel','frc','sid')


class NAPSystem(object):
    """
    NAPSystem class that contains cell information and atoms, and provides some functionalities.
    Atom information is stored as a pandas DataFrame object.
    """
    def __init__(self, fname=None, format=None, specorder=[], ase_atoms=None):
        self.alc = 1.0
        self.a1 = np.zeros(3)
        self.a2 = np.zeros(3)
        self.a3 = np.zeros(3)
        self.specorder = specorder.copy()
        #...DataFrame for all the atoms including atomic data
        self.init_atoms()

        specorder_good = False
        for s in self.specorder:
            if len(s) > 0:
                specorder_good = True
                break
        if not specorder_good:
            self.specorder = None
        
        if fname is not None:
            # self.read(fname=fname,format=format)
            raise ValueError('Initializing NAPSystem with file is obsolete.\n'
                             +'Use nappy.io.read() instead.')
        elif ase_atoms is not None:
            # self.load_ase_atoms(ase_atoms)
            raise ValueError('Initializing NAPSystem with ase_atoms is obsolete.\n'
                             +'Use nappy.io.from_ase(atoms) instead.')
        else:
            pass

        return None

    def __repr__(self):
        txt = analyze_msg(self)
        return txt

    def __len__(self):
        return len(self.atoms)

    def init_atoms(self):
        import pandas as pd
        self.atoms = pd.DataFrame(columns=DEFAULT_LABELS)
        return None

    def set_atoms(self,newatoms):
        self.atoms = copy.deepcopy(newatoms)
        self.atoms.reset_index(drop=True,inplace=True)
        self._reset_atoms_dtypes()
        return None

    def _reset_atoms_dtypes(self):
        """Reset dtypes of atoms dataframe according to the default setting."""
        cols = self.atoms.columns.to_list()
        for c in cols:
            try:
                t = DEFAULT_DTYPES[c]
                self.atoms[c] = self.atoms[c].astype(t)
            except KeyError as ke:
                pass
        return None
        
    def get_aux_names(self):
        aux_names = list(self.atoms.columns)
        for l in DEFAULT_LABELS:
            aux_names.remove(l)
        if 'neighbors' in aux_names:
            aux_names.remove('neighbors')
        return aux_names

    def set_lattice_constant(self,alc):
        self.alc = alc
        return None

    def set_lattice(self, alc, a1, a2, a3):
        """
        Set lattice by specifying lattice constant and 3 vectors.

        Args:
          alc (float): lattice constant
          a1 (float vector of 3 components): a1 vector
          a2 (float vector of 3 components): a2 vector
          a3 (float vector of 3 components): a3 vector
        """
        self.alc = alc
        self.a1[:] = a1[:]
        self.a2[:] = a2[:]
        self.a3[:] = a3[:]
        return None

    def set_lattice_parameters(self,a,b,c,alpha,beta,gamma):
        """
        Set lattice by specifying lattice parameters, a,b,c,alpha,beta,gamma.
        See https://arxiv.org/pdf/1506.01455.pdf

        Args:
          a (float): lattice parameter *a*
          b (float): lattice parameter *b*
          c (float): lattice parameter *c*
          alpha (float): angle in degree
          beta (float): angle in degree
          gamma (float): angle in degree
        """
        self.alc = 1.0
        alpr = np.radians(alpha)
        betr = np.radians(beta)
        gmmr = np.radians(gamma)
        # val = (cos(alpha_r) * cos(beta_r) - cos(gamma_r))\
        #     / (sin(alpha_r) * sin(beta_r))
        # val = max(abs(val),1.0)
        # gamma_star = np.arccos(val)
        self.a1[:] = [float(a), 0.0, 0.0]
        self.a2[:] = [b*np.cos(gmmr),
                      b*np.sin(gmmr), 0.0]
        # self.a3[:] = [c*cos(beta_r),
        #               -c*sin(beta_r)*cos(gamma_star),
        #               c*sin(beta_r)*sin(gamma_star)]
        self.a3[:] = [c*np.cos(betr),
                      c*(np.cos(alpr) -np.cos(betr)*np.cos(gmmr))/np.sin(gmmr),
                      c*np.sqrt(np.sin(gmmr)**2 -np.cos(alpr)**2 -np.cos(betr)**2
                             +2.0*np.cos(alpr)*np.cos(betr)*np.cos(gmmr))/np.sin(gmmr)]
        return None
        
    def get_hmat(self):
        """
        Be careful about the definition of H-matrix here.
        It is transpose of the definition of that in POSCAR.
        Here hmat = [a1,a2,a3] so that
        pos = hmat * spos
        where pos and spos are real Cartesian position and scaled position.
        """
        hmat = np.zeros((3, 3), dtype=float)
        hmat[:, 0] = self.a1 * self.alc
        hmat[:, 1] = self.a2 * self.alc
        hmat[:, 2] = self.a3 * self.alc
        return hmat
    
    def set_hmat(self, hmat):
        self.alc = 1.0
        self.a1[:] = hmat[:, 0]
        self.a2[:] = hmat[:, 1]
        self.a3[:] = hmat[:, 2]

    def get_hmat_inv(self):
        hmat = self.get_hmat()
        return np.linalg.inv(hmat)

    def set_specorder(self,*specorder):
        #...Check if the number of species is relevant
        maxsid = self.atoms.sid.max()
        if len(specorder) < maxsid:
            txt = 'Number of species is not sufficient,' \
                  +' and must be greater than or equal {0:d}'.format(maxsid)
            raise ValueError(txt)
        #...Operation could be different case by case
        if self.specorder is None:
            self.specorder = copy.copy(specorder)
        # elif set(self.specorder) == set(specorder):  # Only re-ordering
        #     newsids = np.zeros(len(self.atoms),dtype=int)
        #     for i,sid in enumerate(self.atoms.sid):
        #         symbol = self.specorder[sid-1]
        #         sidnew = specorder.index(symbol)+1
        #         newsids[i] = sidnew
        #     self.atoms.sid = newsids
        #     self.specorder = specorder
        else:  # Re-define specorder even if specorder and sids are inconsistent...
            newsids = np.zeros(len(self.atoms),dtype=int)
            for i,sid in enumerate(self.atoms.sid):
                symbol = self.specorder[sid-1]
                try:
                    sidnew = specorder.index(symbol)+1
                except:
                    raise
                newsids[i] = sidnew
            self.atoms.sid = newsids
            self.specorder = specorder
        return None
                
    def get_lattice_vectors(self):
        return self.a1*self.alc, self.a2*self.alc, self.a3*self.alc

    def get_lattice_lengths(self):
        a = np.linalg.norm(self.a1*self.alc)
        b = np.linalg.norm(self.a2*self.alc)
        c = np.linalg.norm(self.a3*self.alc)
        return a,b,c

    def get_lattice_angles(self,unit='radius'):
        a = np.linalg.norm(self.a1)
        b = np.linalg.norm(self.a2)
        c = np.linalg.norm(self.a3)
        bc = np.cross(self.a2,self.a3)
        ac = np.cross(self.a1,self.a3)
        ab = np.cross(self.a1,self.a2)
        bc = np.linalg.norm(bc)
        ac = np.linalg.norm(ac)
        ab = np.linalg.norm(ab)
        # make it inside the range of arcsin
        ta = min(1.0-1.0e-10,bc/b/c)
        tb = min(1.0-1.0e-10,ac/a/c)
        tc = min(1.0-1.0e-10,ab/a/b)
        alpha = np.arcsin(ta)
        beta  = np.arcsin(tb)
        gamma = np.arcsin(tc)
        if unit[0] in ('d','D'):
            alpha = alpha/np.pi *180.0
            beta  = beta/np.pi *180.0
            gamma = gamma/np.pi *180.0
        return alpha,beta,gamma

    def get_reciprocal_vectors(self):
        a1,a2,a3 = self.get_lattice_vectors()
        v = np.dot(a1,np.cross(a2,a3))
        b1 = 2.0 *np.pi *np.cross(a2,a3)/v
        b2 = 2.0 *np.pi *np.cross(a3,a1)/v
        b3 = 2.0 *np.pi *np.cross(a1,a2)/v
        return b1,b2,b3

    def add_atoms(self,symbols,poss,vels=[],frcs=[]):
        """
        Add atoms of given symbols, positions, velocities and forces.
        Positions, velocities and forces are assumed to be scaled in lattice vectors.
        """
        import pandas as pd
        #...To remove future warning for pd.concat
        import warnings
        warnings.simplefilter(action='ignore',
                              category=FutureWarning)
        if not self.specorder:
            self.specorder = []
        if type(symbols) not in (list, np.ndarray):
            raise TypeError('symbols must be either list or numpy.ndarray.')
        sids = []
        for symbol in symbols:
            if symbol not in self.specorder:
                self.specorder.append(symbol)
            sid = self.specorder.index(symbol)+1
            sids.append(sid)
        
        newatoms = pd.DataFrame(columns=self.atoms.columns)
        if type(poss) == list:
            poss = np.array(poss)
        newatoms[['x','y','z']] = poss
        if len(vels) == len(poss):
            if type(vels) == list:
                vels = np.array(vels)
        else:
            vels = np.zeros(poss.shape)
        newatoms[['vx','vy','vz']] = vels
        if len(frcs) == len(poss):
            if type(frcs) == list:
                frcs = np.array(frcs)
        else:
            frcs = np.zeros(poss.shape)
        newatoms[['fx','fy','fz']] = frcs
        newatoms.sid = [ sid for sid in sids ]
        self.atoms = pd.concat([self.atoms, newatoms])
        self.atoms.reset_index(drop=True,inplace=True)
        self._reset_atoms_dtypes()
        self.assign_pbc()
        return None

    def remove_atoms(self,*indices):
        """
        Remove atoms of given INDICES from the atoms list.
        """
        try:
            self.atoms.drop(index=list(indices),inplace=True)
            self.atoms.reset_index(drop=True,inplace=True)
        except Exception as e:
            raise

        return None

    def tidy_specorder(self):
        """Tidy up the specorder by removing species of 0 atoms."""
        natoms = self.natm_per_species()
        to_remove = []
        for i,s in enumerate(self.specorder):
            if natoms[i] == 0:
                to_remove.append(s)
        if not len(to_remove) == 0:
            now_specorder = self.specorder
            now_symbols = self.get_symbols()
            new_specorder = copy.copy(now_specorder)
            for r in to_remove:
                new_specorder.remove(r)
            sids = [ new_specorder.index(s)+1 for s in now_symbols ]
            self.specorder = new_specorder
            self.atoms['sid'] = sids

        return None

    def species2sid(self, spc:str) -> int:
        """Convert species name to sid using specorder.
        """
        try:
            return self.specorder.index(spc) +1
        except:
            return -1
    
    def num_atoms(self,sid=0):
        if sid == 0:
            return len(self.atoms)
        else:
            return len(self.atoms[self.atoms.sid==sid])

    def num_species(self):
        """
        Return number of species in the system, counted not using self.specorder.
        """
        return self.atoms.sid.max()
    
    def natm_per_species(self):
        """
        Return number of atoms per species as a list of integer in the order of specorder.
        """
        nps= []
        max_nsp= 0
        try:
            if len(self.specorder) > 0:
                for i in range(len(self.specorder)):
                    nps.append(0)
                for sid in self.atoms.sid:
                    nps[sid-1] += 1
            else:
                max_nsp = self.atoms.sid.max()
                for i in range(max_nsp):
                    nps.append(0)
                for sid in self.atoms.sid:
                    nps[sid-1] += 1
        except Exception:
            nps.append(len(self.atoms))
        return nps

    def get_volume(self):
        return self.alc**3 *np.abs(np.dot(self.a1,np.cross(self.a2,self.a3)))

    def get_positions(self,real=False):
        """
        Get positions of atoms. If real=True, the positions are multiplied with h-mat.
        """
        if real:
            return self.get_real_positions()
        return self.get_scaled_positions()

    def get_velocities(self,real=False):
        """
        Get velocities of atoms. If real=True, they are multiplied with h-mat.
        """
        if real:
            return self.get_real_velocities()
        return self.get_scaled_velocities()
    
    def get_forces(self,real=False):
        """
        Get forces of atoms. If real=True, they are multiplied with h-mat.
        """
        if real:
            return self.get_real_forces()
        return self.get_scaled_forces()
    
    def get_real_positions(self):
        hmat = self.get_hmat()
        spos = self.get_scaled_positions()
        rpos = np.zeros((len(self.atoms),3))
        for ia in range(len(self.atoms)):
            pos = spos[ia]
            rpos[ia,:] = np.dot(hmat,pos)
        return rpos

    def set_real_positions(self,rposs):
        if len(rposs) != len(self.atoms):
            raise ValueError('Array size inconsistent.')
        hmati = self.get_hmat_inv()
        for i in range(len(self.atoms)):
            rpi = rposs[i]
            spi = np.dot(hmati,rpi)
            self.atoms.at[i,['x','y','z']] = spi
        return None

    def get_scaled_positions(self):
        return self.atoms[['x','y','z']].to_numpy(dtype=float)

    def set_scaled_positions(self,sposs):
        if len(sposs) != len(self.atoms):
            raise ValueError('Array size inconsistent.')
        if type(sposs) == list:
            sposs = np.array(sposs)
        self.atoms[['x','y','z']] = sposs
        return None

    def get_scaled_velocities(self):
        return self.atoms[['vx','vy','vz']].to_numpy(dtype=float)

    def get_real_velocities(self):
        hmat = self.get_hmat()
        vels = self.get_scaled_velocities()
        rvels = np.zeros(vels.shape)
        for ia in range(len(vels)):
            v = vels[ia]
            rvels[ia,:] = np.dot(hmat,v)
        return rvels

    def set_scaled_velocities(self,svels):
        if len(svels) != len(self.atoms):
            raise ValueError('Array size inconsistent.')
        if type(svels) == list:
            svels = np.array(svels)
        self.atoms[['vx','vy','vz']] = svels
        return None

    def set_real_velocities(self,vels):
        if len(vels) != len(self.atoms):
            raise ValueError('Array size inconsistent.')
        if type(vels) == list:
            vels = np.array(vels)
        svels = np.zeros(vels.shape)
        hmati = self.get_hmat_inv()
        for ia in range(len(self.atoms)):
            svels[ia,:] = np.dot(hmati,vels[ia,:])
        self.atoms[['vx','vy','vz']] = svels
        return None

    def get_scaled_forces(self):
        return self.atoms[['fx','fy','fz']].to_numpy(dtype=float)

    def set_scaled_forces(self,sfrcs):
        assert len(sfrcs) == len(self.atoms), 'Array size inconsistent.'
        if type(sfrcs) == list:
            sfrcs = np.array(sfrcs)
        self.atoms[['fx','fy','fz']] = sfrcs
        return None
    
    def get_real_forces(self):
        hmat = self.get_hmat()
        frcs = self.get_scaled_forces()
        rfrcs = np.zeros(frcs.shape)
        for ia in range(len(frcs)):
            f = frcs[ia]
            rfrcs[ia,:] = np.dot(hmat,f)
        return rfrcs

    def set_real_forces(self,frcs):
        assert len(frcs) == len(self.atoms), 'Array size inconsistent.'
        if type(frcs) == list:
            frcs = np.array(frcs)
        sfrcs = np.zeros(frcs.shape)
        hmati = self.get_hmat_inv()
        for ia in range(len(self.atoms)):
            sfrcs[ia,:] = np.dot(hmati, frcs[ia,:])
        self.atoms[['fx','fy','fz']] = sfrcs
        return None

    def get_tags(self):
        """
        Returns
        -------
        tags : numpy.array
            List of tags of atoms which is used in pmd.
        """
        tags = np.zeros(self.num_atoms())
        sids = self.atoms.sid
        for i in range(self.num_atoms()):
            tags[i] = 1.0*sids[i] +0.1 +1e-14*(i+1)
        return tags
        
    def get_symbols(self):
        """
        Returns
        -------
        symbols : list
              List of chemical symbols of all atoms in the system.
        """
        if not self.specorder:
            raise ValueError('specorder is not available.')
        symbols = []
        for sid in self.atoms.sid:
            symbols.append(self.specorder[sid-1])
        return symbols

    def get_symbol_of(self,ia):
        """
        Returns
        -------
        symbol : string
              The symbol of a specified atoms-ia.
        """
        sid = self.atoms.loc[ia,'sid']
        return self.specorder[sid-1]

    def set_symbols(self,symbols):
        """
        Set symbols of all atoms and append new symbols
        if they are not in the current specorder.
        """
        if len(symbols) != len(self.atoms):
            raise RuntimeError('len(symbols) != len(self.atoms)')
        for s in symbols:
            if s not in self.specorder:
                self.specorder.append(s)
        for i in range(len(self.atoms)):
            s = symbols[i]
            sid = self.specorder.index(s)+1
            self.atoms.at[i,'sid'] = sid
        return None

    def get_chemical_formula(self):
        """
        Returns chemical formula as a string based on the chemical symbols same as ASE.
        """
        symbols = self.get_symbols()
        
        uniq_symbols = []
        for s in symbols:
            if s not in uniq_symbols:
                uniq_symbols.append(s)
        formula = ''
        for s in uniq_symbols:
            n = symbols.count(s)
            formula += s + '{0:d}'.format(n)
        return formula

    def get_charges(self):
        # return self.charges
        try:
            chgs = np.array(self.atoms.chg.values)
        except Error:
            raise
        return chgs

    def set_charges(self,charges):
        # self.charges = charges
        if len(charges) != len(self.atoms):
            raise RuntimeError('len(charges) != len(self.atoms): ',len(charges),len(self.atoms))
        self.atoms['chg'] = charges
        # if len(self.charges) > 0:
        #     if self.specorder is not None \
        #        and len(self.charges) < len(self.specorder):
        #         lenc = len(self.charges)
        #         for i in range(len(self.specorder)-lenc):
        #             self.charges.append(0.0)
        return None

    def set_potential_energy(self,epot):
        self.epot = epot
        return None

    def set_kinetic_energy(self,ekin):
        self.ekin = ekin
        return None

    def set_stress_tensor(self,stnsr):
        if stnsr.shape != (3,3):
            raise TypeError('Stress tensor should be 3x3 array.')
        self.stnsr = copy.copy(stnsr)
        return None
    
    def get_potential_energy(self):
        return getattr(self,'epot',None)

    def get_kinetic_energy(self):
        return getattr(self,'ekin',None)

    def get_stress_tensor(self):
        return getattr(self,'stnsr',None)

    def get_stress(self):
        return np.array([self.stnsr[0,0],
                         self.stnsr[1,1],
                         self.stnsr[2,2],
                         self.stnsr[1,2],
                         self.stnsr[0,2],
                         self.stnsr[0,1]])
    
    def get_atom(self,idatm=-1):
        """
        Return a dictionary data of an atom specified by IDATM.
        Since this is rather slow, use get_atom_attr() when you want a specifiec attribute of an atom,
        """
        if idatm < 0 or idatm >= len(self.atoms):
            raise ValueError('idatm must be specified within, ',0,len(self.atoms)-1)
        return self.atoms.iloc[idatm].to_dict()

        
    def get_atom_attr(self,idatm=-1,attr_name=None):
        """
        Returns an atom attribute of given IDATM.
        Possible attributes can be:
          - pos
          - vel
          - frc
          - sid
          - other data such as chg, epot, ekin if exist
        """
        if attr_name is None:
            raise ValueError('attr_name must be specified.')
        if idatm < 0 or idatm >= len(self.atoms):
            raise ValueError('idatm must be specified within, ',0,len(self.atoms)-1)
        if attr_name not in self.atoms.columns:
            raise ValueError('attr_name must be in the list,',self.atoms.columns)

        return self.atoms.loc[idatm,attr_name]

    def get_density(self,):
        from nappy.elements import elements
        from nappy.units import amu_to_g, Ang_to_cm
        mass = 0.0
        vol = self.get_volume()
        nspcs = self.natm_per_species()
        for i,s in enumerate(self.specorder):
            mass += nspcs[i] *elements[s]['mass']
        return mass*amu_to_g/(vol*Ang_to_cm**3)

    def get_num_density(self, species=None):
        from nappy.elements import elements
        from nappy.units import amu_to_g, Ang_to_cm
        mass = 0.0
        vol = self.get_volume()
        nspcs = self.natm_per_species()
        if species:
            try:
                sid = self.specorder.index(species)
            except Exception:
                raise
            return float(nspcs[sid])/vol
        return float(len(self))/vol

    def set_voldata(self, voldata, nvoldiv=None, volorig=None):
        """
        Store volumetric data in the system.
        """
        try:
            if len(voldata.shape) == 1:
                if nvoldiv is None:
                    raise ValueError('If voldata is 1D array, '
                                     +'nvoldiv should be provided.')
                else:
                    self.nvoldiv = nvoldiv
                self.voldata = voldata.reshape(nvoldiv)
            elif len(voldata.shape) == 3:
                self.nvoldiv = voldata.shape
                self.voldata = voldata
        except Exception as e:
            raise

        if volorig == None:
            self.vorig = np.zeros(3)
        else:
            self.vorig = volorig
        
        return None
    
    def view(self,backend='3Dmol',**options):
        """
        Visualize the system using the given backend.
        """
        
        if '3Dmol' in backend:
            self._view_3Dmol(**options)
        else:
            raise ValueError('No such visualization backend.')
        return None

    def _view_3Dmol(self,**options):
        """
        Visualize the system using 3Dmol.js and py3Dmol.
        """
        from nappy.elements import elements
        try:
            import py3Dmol
        except Exception as e:
            raise

        if 'model' in options.keys():
            modopts = options['model']
        if 'volume' in options.keys():
            volopts = options['volume']
        
        pdb = nappy.io.get_PDB_txt(self,**options)
        v = py3Dmol.view()
        v.removeAllModels()
        v.setViewStyle({'style':'outline','color':'black','width':0.03})
        v.addModel(pdb,'pdb')
        colors = { s:elements[s]['CPKcolor'] for s in self.specorder }
        if 'modopts' in locals():
            v.setStyle(modopts)
        else:
            v.setStyle({'stick':{'radius':.1,'colorscheme':{'prop':'elem','map':colors}},
                        'sphere':{'radius':.5,'colorscheme':{'prop':'elem','map':colors}}})
        v.addUnitCell()
        if hasattr(self,'voldata'):
            cube = nappy.io.get_cube_txt(self)
            if 'volopts' in locals():
                v.addVolumetricData(cube, "cube", volopts)
            else:
                v.addVolumetricData(cube, "cube", {'isoval': 300, 'color': "red", 'opacity': 0.95})
        v.zoomTo()
        v.show()

        return None

    def get_distance(self,ia,ja):
        """
        Compute distance between atoms ia and ja taking the periodic boundary
        condition into account.
        """
        natm = len(self.atoms)
        if ia > natm:
            raise ValueError('ia > natms, ia,natms = ',ia,natm)
        if ja > natm:
            raise ValueError('ja > natms, ja,natms = ',ja,natm)
        xi = self.atoms.loc[ia,['x','y','z']].to_numpy(dtype=float)
        xj = self.atoms.loc[ja,['x','y','z']].to_numpy(dtype=float)
        xij = xj-xi -np.rint(xj-xi)
        hmat = self.get_hmat()
        rij = np.dot(hmat,xij)
        rij2 = rij[0]**2 +rij[1]**2 +rij[2]**2
        return np.sqrt(rij2)

    def get_angle(self,i,j,k):
        """
        Compute angle in degree between bonds i-j and i-k.
        """
        natm = len(self.atoms)
        if i > natm:
            raise ValueError('i > natms, i,natms = ',i,natm)
        if j > natm:
            raise ValueError('j > natms, j,natms = ',j,natm)
        if k > natm:
            raise ValueError('k > natms, k,natms = ',k,natm)
        xi = self.atoms.loc[i,['x','y','z']].to_numpy(dtype=float)
        xj = self.atoms.loc[j,['x','y','z']].to_numpy(dtype=float)
        xk = self.atoms.loc[k,['x','y','z']].to_numpy(dtype=float)
        xij = xj-xi -np.round(xj-xi)
        xik = xk-xi -np.round(xk-xi)
        hmat = self.get_hmat()
        rij = np.dot(hmat,xij)
        rik = np.dot(hmat,xik)
        dij = np.linalg.norm(rij)
        dik = np.linalg.norm(rik)
        angle = np.arccos(np.dot(rij,rik)/dij/dik) /np.pi *180.0
        return angle

    def make_pair_list(self,rcut=3.0,rcuts=None,nnmax=100):
        """
        Make a neighbor list.
        The neighbor list of each atom is stored as neighbors in the self.atoms column.

        INPUT:
            rcut: float
                Cutoff radius for all neighbors to be taken into account.
            rcuts: dict
                Dictionary of pairwise cutoff values. The format is like below.
                ```
                rcuts = {('Si','Si'):2.5, ('Si','O'):2.0,}
                ```
                If a pair is not given, the cutoff for the pair is the maximum of those.
            nnmax: int
                Max number of neighbors. [default: 100]
        """
        try:
            import nappy.pmd.pairlist as pl
            plst = pl.fmake_pairlist(self,rcut=rcut,nnmax=nnmax)
            lspr = []
            for ia in range(self.num_atoms()):
                lspr.append( [ plst[ia,1+j]-1 for j in range(plst[ia,0]) ] )
            self.atoms['neighbors'] = lspr
            return None
        except:
            pass
        rcs2 = np.zeros((len(self.specorder),len(self.specorder)),dtype=float)
        if rcuts is not None:
            for i,si in enumerate(self.specorder):
                for j,sj in enumerate(self.specorder):
                    if j < i: continue
                    if (si,sj) in rcuts:
                        rc2 = rcuts[(si,sj)]**2
                    elif (sj,si) in rcuts:
                        rc2 = rcuts[(sj,si)]**2
                    else:
                        rc2 = -1.0
                    rcs2[i,j] = rc2
                    rcs2[j,i] = rc2
            rcmax2 = max(rcs2.max(),rcut**2)
            for i in range(len(self.specorder)):
                for j in range(len(self.specorder)):
                    if rcs2[i,j] < 0.0:
                        rcs2[i,j] = rcmax2
            rcut = np.sqrt(rcmax2)
        else:
            rcs2[:,:] = rcut**2
        rc2= rcut**2

        h= self.get_hmat()
        hi= np.linalg.inv(h)
        lcx= int(1.0/np.sqrt(hi[0,0]**2 +hi[1,0]**2 +hi[2,0]**2)/rcut)
        lcy= int(1.0/np.sqrt(hi[0,1]**2 +hi[1,1]**2 +hi[2,1]**2)/rcut)
        lcz= int(1.0/np.sqrt(hi[0,2]**2 +hi[1,2]**2 +hi[2,2]**2)/rcut)
        if lcx == 0: lcx= 1
        if lcy == 0: lcy= 1
        if lcz == 0: lcz= 1
        lcyz= lcy*lcz
        lcxyz= lcx*lcy*lcz
        rcx= 1.0/lcx
        rcy= 1.0/lcy
        rcz= 1.0/lcz
        rcxi= 1.0/rcx
        rcyi= 1.0/rcy
        rczi= 1.0/rcz
        lscl= np.zeros((len(self.atoms),),dtype=int)
        lshd= np.zeros((lcxyz,),dtype=int)
        lscl[:]= -1
        lshd[:]= -1
        
        # Use numpy array instead of accessing pandas series when it will be heavily accessed.
        poss = self.get_scaled_positions()

        #...make a linked-cell list
        self.assign_pbc()
        for i in range(self.num_atoms()):
            pi = poss[i]
            #...assign a vector cell index
            mx = int(pi[0]*rcxi)
            my = int(pi[1]*rcyi)
            mz = int(pi[2]*rczi)
            mx = min(max(mx,0),lcx-1)
            my = min(max(my,0),lcy-1)
            mz = min(max(mz,0),lcz-1)
            m= mx*lcyz +my*lcz +mz
            lscl[i]= lshd[m]
            lshd[m]= i

        #...Determine possible max of num of neighbors
        nnmax = 0
        for icell in range(len(lshd)):
            inc = 0
            i = lshd[icell]
            if i == -1: continue
            while i >= 0:
                inc += 1
                i = lscl[i]
            nnmax = max(nnmax,inc)

        #...Initialize lspr
        # lspr = [ [] for i in range(self.num_atoms()) ]
        nplspr = np.zeros((self.num_atoms(),nnmax*27),dtype=int)
        nplspr[:,:] = -1
        nlspr = np.zeros(self.num_atoms(),dtype=int)
        # self.atoms['neighbors'] = emptylist
        sids = self.atoms.sid

        for ia in range(self.num_atoms()):
            pi = poss[ia]
            #...assign a vector cell index
            mx = int(pi[0]*rcxi)
            my = int(pi[1]*rcyi)
            mz = int(pi[2]*rczi)
            mx = min(max(mx,0),lcx-1)
            my = min(max(my,0),lcy-1)
            mz = min(max(mz,0),lcz-1)
            isp = sids[ia] -1
            for kuz in (-1,0,1):
                m1z = mz +kuz
                if m1z < 0: m1z += lcz
                if m1z >= lcz: m1z -= lcz
                for kuy in (-1,0,1):
                    m1y= my +kuy
                    if m1y < 0: m1y += lcy
                    if m1y >= lcy: m1y -= lcy
                    for kux in (-1,0,1):
                        m1x= mx +kux
                        if m1x < 0: m1x += lcx
                        if m1x >= lcx: m1x -= lcx
                        m1= m1x*lcyz +m1y*lcz +m1z
                        if lshd[m1] == -1: continue

                        ja = lshd[m1]
                        while ja >= 0:
                            if ja <= ia:
                                ja = lscl[ja]
                                continue
                            jsp = sids[ja] -1
                            pij = poss[ja] -pi
                            pij = pij - np.round(pij)
                            rij = np.dot(h,pij)
                            rij2 = rij[0]**2 +rij[1]**2 +rij[2]**2
                            if rij2 < rcs2[isp,jsp] and ja not in nplspr[ia,:]:
                                nplspr[ia,nlspr[ia]] = ja
                                nplspr[ja,nlspr[ja]] = ia
                                dij = np.sqrt(rij2)
                                nlspr[ia] += 1
                                nlspr[ja] += 1
                                
                            ja = lscl[ja]
        #...Finally add the lspr to atoms DataFrame
        lspr = []
        for ia in range(self.num_atoms()):
            lspr.append([ nplspr[ia,ja] for ja in range(nlspr[ia]) ])
        self.atoms['neighbors'] = lspr
        return None

    def remove_pair_list(self):
        try:
            del self.atoms['neighbors']
        except:
            raise
        return None

    def neighbors_of(self,ia,rcut=3.0):
        """
        Generator of the neighbors of a given atom-i.
        """
        if 'neighbors' not in self.atoms.columns:
            self.make_pair_list(rcut=rcut)
        lspri = self.atoms['neighbors'][ia]
        for jj in range(len(lspri)):
            yield lspri[jj]

    def assign_pbc(self):
        poss = self.get_scaled_positions()
        for i in range(len(self.atoms)):
            pi = poss[i]
            newpi = np.zeros(3)
            newpi[0] = pbc(pi[0])
            newpi[1] = pbc(pi[1])
            newpi[2] = pbc(pi[2])
            #self.atoms.at[i,['x','y','z']] = newpi
            self.atoms.loc[i,['x','y','z']] = newpi
        return None

    def get_expansion_num(self,length):
        """
        Compute expansion digits along a1, a2, and a3 from a given *length*.
        System size should be larger than the *length*.
        """
        h= np.zeros((3,3),dtype=float)
        h[0]= self.a1 *self.alc
        h[1]= self.a2 *self.alc
        h[2]= self.a3 *self.alc
        vol= np.abs(np.dot(h[0],np.cross(h[1],h[2])))
        l1= np.abs(vol/(np.linalg.norm(np.cross(h[1],h[2]))))
        l2= np.abs(vol/(np.linalg.norm(np.cross(h[2],h[0]))))
        l3= np.abs(vol/(np.linalg.norm(np.cross(h[0],h[1]))))
        # print('alc,a1,a2,a3=',self.alc,self.a1,self.a2,self.a3)
        # print('h=',h)
        # print('vol,l1,l2,l3=',vol,l1,l2,l3)
        n1= int(np.ceil(length/l1))
        n2= int(np.ceil(length/l2))
        n3= int(np.ceil(length/l3))
        return n1,n2,n3

    def shift_atoms(self,sx,sy,sz):
        """
        Shift all the atoms by (sx,sy,sz) in scaled unit.
        """
        def _shift(x,d):
            return x+d
        self.atoms['x'] = self.atoms['x'].apply(_shift,d=sx)
        self.atoms['y'] = self.atoms['y'].apply(_shift,d=sy)
        self.atoms['z'] = self.atoms['z'].apply(_shift,d=sz)
        return None

    def cycle_coord(self,ncycle=0):
        """
        Cyclic shift of coordinates, e.g. x-y-z to y-z-x.
        """
        ncycle = ncycle % 3
        if ncycle == 0:
            return None
        for icycle in range(ncycle):
            a1,a2,a3 = self.get_lattice_vectors()
            a1 = [ a1[1], a1[2], a1[0] ]
            a2 = [ a2[1], a2[2], a2[0] ]
            a3 = [ a3[1], a3[2], a3[0] ]
            self.set_lattice(1.0, a2, a3, a1)

        from_col = ['x','y','z']
        to_col = [ from_col[ (0+ncycle) %3],
                   from_col[ (1+ncycle) %3],
                   from_col[ (2+ncycle) %3],]
        newdf = np.DataFrame(columns=['x','y','z'])
        for i in range(3):
            newdf[to_col[i]] = self.atoms[from_col[i]]
        self.atoms[['x','y','z']] = newdf[['x','y','z']]
        return None
        
    def repeat(self,n1o,n2o,n3o,n1m=0,n2m=0,n3m=0):
        """
        Multiply the system by given n1o,n2o,n3o and replace the system 
        with multiplied one.
        """
        import pandas as pd
        #...Convert to int
        n1 = int(n1o)
        n2 = int(n2o)
        n3 = int(n3o)
        if n1 == 0: n1 = 1
        if n2 == 0: n2 = 1
        if n3 == 0: n3 = 1
        if n1 == n2 == n3 == 1:
            return None
        #...unit vectors to be repeated
        m1 = n1-n1m
        m2 = n2-n2m
        m3 = n3-n3m
        self.a1= self.a1*m1
        self.a2= self.a2*m2
        self.a3= self.a3*m3
        #n123= m1*m2*m3
        maxsid = self.atoms.sid.max()
        #natm0= self.num_atoms()
        # atoms0= copy.copy(self.atoms)
        newnatm = len(self.atoms) *n1*n2*n3
        newsids = [ 0 for i in range(newnatm) ]
        newposs = np.zeros((newnatm,3))
        newvels = np.zeros((newnatm,3))
        newfrcs = np.zeros((newnatm,3))
        colnames = list(self.atoms.columns)
        #...Labels except (sid,pos,vel,frc) are all auxiliary data
        auxnames = colnames.copy()
        auxnames.remove('sid')
        auxnames.remove('x')
        auxnames.remove('y')
        auxnames.remove('z')
        auxnames.remove('vx')
        auxnames.remove('vy')
        auxnames.remove('vz')
        auxnames.remove('fx')
        auxnames.remove('fy')
        auxnames.remove('fz')
        newauxs = {}
        for auxname in auxnames:
            newauxs[auxname] = []
        inc = 0
        poss = self.get_scaled_positions()
        vels = self.atoms[['vx','vy','vz']].to_numpy(dtype=float)
        frcs = self.atoms[['fx','fy','fz']].to_numpy(dtype=float)
        for i1 in range(n1m,n1):
            for i2 in range(n2m,n2):
                for i3 in range(n3m,n3):
                    for i0 in range(len(self.atoms)):
                        pi0 = poss[i0]
                        x= pi0[0]/m1 +1.0/m1*i1
                        y= pi0[1]/m2 +1.0/m2*i2
                        z= pi0[2]/m3 +1.0/m3*i3
                        vi0 = vels[i0]
                        vx= vi0[0]/m1 +1.0/m1*i1
                        vy= vi0[1]/m2 +1.0/m2*i2
                        vz= vi0[2]/m3 +1.0/m3*i3
                        newsids[inc] = self.atoms.sid[i0]
                        newposs[inc,:] = [x,y,z]
                        newvels[inc,:] = [vx,vy,vz]
                        newfrcs[inc,:] = frcs[i0,:]
                        for auxname in auxnames:
                            newauxs[auxname].append(self.atoms.loc[i0,auxname])
                        inc += 1
        #...Use DataFrame self.atoms
        self.atoms = pd.DataFrame(columns=colnames)
        self.atoms[['x','y','z']] = newposs
        self.atoms[['vx','vy','vz']] = newvels
        self.atoms[['fx','fy','fz']] = newfrcs
        # self.atoms['pos'] = newposs
        # self.atoms['vel'] = newvels
        # self.atoms['frc'] = newfrcs
        self.atoms['sid'] = newsids
        for auxname in auxnames:
            self.atoms[auxname] = newauxs[auxname]
        return None

    def divide(self,*ds):
        """
        Divide lattice vectors by (d1,d2,d3).
        Atoms whose scaled positions are greater than ds[#] are to be removed.
        """
        import pandas as pd
        if len(ds) != 3:
            raise ValueError('len(ds) != 3.')

        #...1st, count numbers to remove
        nrm = 0
        spos = self.get_scaled_positions()
        for ia in range(len(self.atoms)):
            pi = spos[ia]
            for l in range(3):
                if pi[l] >= ds[l]:
                    nrm += 1
                    break

        newnatm = len(self.atoms) -nrm
        self.a1 = self.a1 *ds[0]
        self.a2 = self.a2 *ds[1]
        self.a3 = self.a3 *ds[2]
        newsids = [ 0 for i in range(newnatm) ] 
        newposs = np.zeros((newnatm,3))
        newvels = np.zeros((newnatm,3))
        newfrcs = np.zeros((newnatm,3))
        # newposs = [ np.zeros(3) for i in range(newnatm) ] 
        # newvels = [ np.zeros(3) for i in range(newnatm) ] 
        # newfrcs = [ np.zeros(3) for i in range(newnatm) ] 
        colnames = list(self.atoms.columns)
        #...Labels except (sid,pos,vel,frc) are all auxiliary data
        auxnames = colnames.copy()
        auxnames.remove('sid')
        auxnames.remove('x')
        auxnames.remove('y')
        auxnames.remove('z')
        auxnames.remove('vx')
        auxnames.remove('vy')
        auxnames.remove('vz')
        auxnames.remove('fx')
        auxnames.remove('fy')
        auxnames.remove('fz')
        newauxs = {}
        for auxname in auxnames:
            newauxs[auxname] = []
        inc = 0
        poss = self.get_scaled_positions()
        vels = self.atoms[['vx','vy','vz']].to_numpy(dtype=float)
        frcs = self.atoms[['fx','fy','fz']].to_numpy(dtype=float)
        for ia in range(len(self.atoms)):
            pi = poss[ia]
            survive = True
            for l in range(3):
                if pi[l] >= ds[l]:
                    survive = False
                    break
            if survive:
                pinew = np.array(pi)
                for l in range(3):
                    if ds[l] < 0.9:
                        pinew[l] = pinew[l] /ds[l]
                newsids[inc] = self.atoms.sid[ia]
                newposs[inc,:] = pinew[:]
                newvels[inc,:] = vels[ia,:]
                newfrcs[inc,:] = frcs[ia,:]
                for auxname in auxnames:
                    newauxs[auxname].append(self.atoms.loc[ia,auxname])
                inc += 1
        #...Use DataFrame self.atoms
        self.atoms = pd.DataFrame(columns=colnames)
        self.atoms[['x','y','z']] = newposs
        self.atoms[['vz','vy','vz']] = newvels
        self.atoms[['fx','fy','fz']] = newfrcs
        self.atoms['sid'] = newsids
        for auxname in auxnames:
            self.atoms[auxname] = newauxs[auxname]
        return None

    def add_vacuum(self,axis,length,shift=0.0):
        """
        Add vacuum of the given length to the given axis of the current system.
        Shift value is added to each atom position of the given axis.

        Input
        -----
        axis : int (0, 1, or 2)
            Axis to which vacuum is added.
        length : float (in Ang.)
            Length of vacuum region.
        shift : float (in Ang.)
            Shift value that is add to atoms towards vacuum region.
        """
        lengths = self.get_lattice_lengths()
        ratios = np.array((1.0, 1.0, 1.0))
        ratios[axis] += length /lengths[axis]
            
        self.assign_pbc()
        #...Store reaal positions before extending the system
        rpos = self.get_real_positions()
        #...Extend the system cell; atoms is now positioned at the bottom side of extended axis
        self.a1 *= ratios[0]
        self.a2 *= ratios[1]
        self.a3 *= ratios[2]
        hmati = self.get_hmat_inv()
        #...Calc scaled positions in the new system
        spos = np.zeros((len(self.atoms),3),dtype=float)
        #...Get min and max of atoms in the new system cell
        amin = 1.0
        amax = 0.0
        for ia in range(len(self.atoms)):
            xi = rpos[ia]
            spos[ia] = np.dot(hmati,xi)
            amin = min(amin,spos[ia,axis])
            amax = max(amax,spos[ia,axis])
        #...Compute scaled shift value
        lengths2 = self.get_lattice_lengths()
        sshift = shift /lengths2[axis]
        for ia in range(len(self.atoms)):
            spos[ia,axis] += sshift
        self.atoms[['x','y','z']] = spos
        return None
        
    def to_ase_atoms(self):
        """
        Convert NAPSystem object to ASE atoms.
        Note that some information will be abandonned, 
        for example, forces cannot be stored in ASE atoms object (calculator object can have them, though).
        """
        try:
            from ase import Atoms
        except ImportError:
            raise ImportError('Cannot load ase Atoms.')

        cell = np.array([self.a1, self.a2, self.a3])
        cell *= self.alc
        a = np.linalg.norm(cell[0,:])
        b = np.linalg.norm(cell[1,:])
        c = np.linalg.norm(cell[2,:])
        spos = list(self.get_scaled_positions())
        symbols = [ self.specorder[sid-1] for sid in self.atoms.sid ]
        atoms = Atoms(symbols=symbols,
                      cell=cell,
                      scaled_positions=spos,
                      pbc=True)
        vels = [ [vel[0]*a, vel[1]*b, vel[2]*c] for vel in self.atoms[['vx','vy','vz']].values ]
        atoms.set_velocities(vels)

        return atoms

    def remove_overlapping_atoms(self,criterion=0.01):
        """
        Remove overlapping atoms in the system.
        Atoms with bigger indices in overlapping atoms are removed.
        Judgement of overlap is done by comparing the distance between
        atoms and *criterion*.
        """
        if 'neighbors' not in self.atoms.columns:
            self.make_pair_list(rcut=1.0)
        remove_ids = []
        h = np.zeros((3,3),dtype=float)
        h[:,0] = self.a1 *self.alc
        h[:,1] = self.a2 *self.alc
        h[:,2] = self.a3 *self.alc
        cr2 = criterion*criterion
        spos = self.get_scaled_positions()
        for i in range(len(self.atoms)):
            pi = spos[i]
            for jj in range(len(self.atoms.lspr[i])):
                j = self.atoms.lspr[i][jj]
                pj = spos[j]
                xij = pj-pi
                xij = xij -np.round(xij)
                rij = np.dot(h,xij)
                rij2= np.linalg.norm(rij)
                if rij2 < cr2:
                    remove_ids.append(max(i,j))
        # print remove_ids
        self.remove_atoms(*remove_ids)
        return
        

def analyze_msg(nsys):
    from nappy.elements import elements
    from nappy.units import amu_to_g, Ang_to_cm
    msg = ""
    a1 = nsys.a1 *nsys.alc
    a2 = nsys.a2 *nsys.alc
    a3 = nsys.a3 *nsys.alc
    a = np.linalg.norm(a1)
    b = np.linalg.norm(a2)
    c = np.linalg.norm(a3)
    vol = nsys.get_volume()
    alpha = np.arccos(np.dot(a2,a3)/b/c)/np.pi*180.0
    beta  = np.arccos(np.dot(a1,a3)/a/c)/np.pi*180.0
    gamma = np.arccos(np.dot(a1,a2)/a/b)/np.pi*180.0
    msg +=' a1 vector = [{0:10.3f}, {1:10.3f}, {2:10.3f}]\n'.format(a1[0],
                                                                    a1[1],
                                                                    a1[2])
    msg +=' a2 vector = [{0:10.3f}, {1:10.3f}, {2:10.3f}]\n'.format(a2[0],
                                                                    a2[1],
                                                                    a2[2])
    msg +=' a3 vector = [{0:10.3f}, {1:10.3f}, {2:10.3f}]\n'.format(a3[0],
                                                                    a3[1],
                                                                    a3[2])
    msg +=' a = {0:10.3f} A\n'.format(a)
    msg +=' b = {0:10.3f} A\n'.format(b)
    msg +=' c = {0:10.3f} A\n'.format(c)
    msg +=' alpha = {0:7.2f} deg.\n'.format(alpha)
    msg +=' beta  = {0:7.2f} deg.\n'.format(beta)
    msg +=' gamma = {0:7.2f} deg.\n'.format(gamma)
    msg +=' volume= {0:10.3f} A^3\n'.format(vol)
    msg +=' number of atoms   = {0:d}\n'.format(nsys.num_atoms())
    if nsys.specorder:
        msg +=' number of atoms per species:\n'
        nspcs = nsys.natm_per_species()
        mass = 0.0
        for i,s in enumerate(nsys.specorder):
            msg +='   {0:<2s}: {1:>4d}\n'.format(s,nspcs[i])
            mass += nspcs[i] *elements[s]['mass']
        msg +=' density = {0:5.3f} g/cm^3'.format(mass*amu_to_g
                                                  /(vol*Ang_to_cm**3))
        msg +=' = {0:7.5f} atom/Ang^3\n'.format(float(len(nsys))/vol)
    return msg
    
def analyze(nsys):
    msg = analyze_msg(nsys)
    print(msg,end='')
    return None

def main():
    import os,sys
    args = docopt(__doc__.format(os.path.basename(sys.argv[0])))

    infmt= args['--in-format']
    outfmt= args['--out-format']
    infname= args['INFILE']
    outfname= args['OUTFILE']
    scalefactor= args['--scale']
    shift= [ float(s) for s in args['--shift'].split(',') ]
    ncycle= int(args['--cycle-coord'])
    specorder= args['--specorder'].split(',')
    if specorder == 'None' or 'None' in specorder:
        specorder = []
    copies= [ float(i) for i in args['--periodic-copy'].split(',') ]
    charges= args['--charges']
    if charges == 'None':
        charges = []
    else:
        charges = [ float(c) for c in charges.split(',') ]

    nsyss = nappy.io.read(fname=infname,format=infmt,specorder=specorder)
    if type(nsyss) != list:
        nsyss = [nsyss,]

    if args['analyze']:
        analyze(nsyss[0])

    elif args['convert']:
        postfix_num = False
        if len(nsyss) > 1:
            postfix_num = True
            print(' Since the input file contains more than 1 system,'
                  +' files with numbers are to be written.')
        for i in range(len(nsyss)):
            if scalefactor != "None":
                nsyss[i].alc *= float(scalefactor)
            nsyss[i].shift_atoms(*shift)
            if ncycle > 0:
                nsyss[i].cycle_coord(ncycle)
        
            #...Periodic copy if needed
            copy_needed = False
            divide_needed = False
            for c in copies:
                if c > 1.5:
                    copy_needed = True
                elif c < 0.9:
                    divide_needed = True
            if copy_needed:
                nsyss[i].repeat(*copies)
            if divide_needed:
                nsyss[i].divide(*copies)
            ofname = outfname
            if postfix_num:
                ofname += f'_{i:d}'
            nappy.io.write(nsyss[i],fname=ofname,format=outfmt)
    else:
        raise NotImplementedError()

    return None

if __name__ == "__main__":
    main()
