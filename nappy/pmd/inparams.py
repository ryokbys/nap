"""
Module for handling in.params.XXX file and related information.
"""
import os
import copy
import datetime

class Params():
    """
    Parent class for all the params classes.
    """
    self.specorder = []
    self.params = []
    self.pairs = []
    self.triplets = []
    
    def __init__(self,specorder=[]):
        if type(specorder) != list:
            raise TypeError('specorder must be a list.')
        self.specorder = specorder
        return None


class Morse(Params):
    """
    Class that treats parameters and in.params.XXX for Morse potential.
    """
    def __init__(self,specorder=[]):
        super().__init__(specorder):
        return None

    def check(self):
        if len(self.pairs) == 0:
            raise ValueError('Pairs not given.')

        if len(self.pairs)*3 != len(self.params):
            raise ValueError('Nums of params and pairs are inconsistent.')

        return None

    def write(self,fname='in.params.Morse'):
        self.check()
        #...Set order of pairs and corresponding order of params.
        morse_prms = {}
        inc = -1
        for p in pairs:
            inc += 1
            d0 = self.params[inc]
            inc += 1
            alp = self.params[inc]
            inc += 1
            rmin = self.params[inc]
            morse_prms[p] = (d0,alp,rmin)
        
        with open(fname,'w') as f:
            f.write('# cspi, cspj,    D,      alpha,  rmin\n')
            for pair in pairs:
                d0,alp,rmin = morse_prms[tuple(pair)]
                f.write('  {0:3s}   {1:3s}'.format(*pair))
                f.write('  {0:7.4f}  {1:7.4f}  {2:7.4f}\n'.format(d0,alp,rmin))
        
        return None
            

    
class Angular(Params):
    """
    Class that treats parameters and in.params.XXX for Angular potential.
    """
    self.cutoff = -1.0   # negative for unset
    
    def __init__(self,specorder=[]):
        super().__init__(specorder):
        return None

    def check(self):
        """
        Check consistency between params and triplets.
        """
        if len(self.triplets) == 0:
            raise ValueError('Triplets not given.')

        if len(self.triplets)*3 != len(self.params):
            raise ValueError('Nums of params and triplets are inconsistent.')

        if self.cutoff < 0.0:
            raise ValueError('Cutoff is not set yet.')
        
        return None

    def write(self,fname='in.params.angular'):
        self.check()

        #...Set order of pairs and corresponding order of params.
        angular_prms = {}
        for t in triplets:
            inc += 1
            alp = self.params[inc]
            inc += 1
            bet = self.params[inc]
            inc += 1
            gmm = self.params[inc]
            angular_prms[tuple(t)] = (self.cutoff,alp,bet,gmm)
        
        with open(fname,'w') as f:
            f.write('# type,   cspi, cspj, cspk,  rc3,   alp,   bet,   gmm\n')
            for t in triplets:
                rc,alp,bet,gmm = angular_prms[tuple(t)]
                f.write(' angular1   {0:3s}   {1:3s}   {2:3s} '.format(*t))
                f.write(' {0:6.2f}  {1:7.3f} {2:7.3f} {3:7.3f}\n'.format(rc,alp,bet,gmm))
        
        return None
            
class Coulomb(Params):
    """
    Class that treats parameters and in.params.XXX for Coulomb potential.
    """
    self.__charge_types__ = ('fixed','fixed_bvs','variable','qeq')
    self.__term_types__ = ('full','short','long','direct','direct_cut','screened_cut')
    self.coul_params = {}
    
    def __init__(self,specorder=[]):
        super().__init__(specorder):
        return None

    def check(self):
        if len(self.pairs) == 0:
            raise ValueError('Pairs not given.')

        if len(self.pairs)*3 != len(self.params):
            raise ValueError('Nums of params and pairs are inconsistent.')

        return None

    def read(self,fname='in.params.Coulomb'):
        """
        Read in.params.Coulomb and store information in self.coul_params.
        """
        if not os.path.exists(fname):
            raise FileNotFoundError(fname)
        
        with open(fname,'r') as f:
            lines = f.readlines()
            pass

        mode = None
        fbvs = 0.0
        rads = {}
        vids = {}
        npqs = {}
        self.coul_params = {}

        for i,line in enumerate(lines):
            data = line.split()
            if len(data) == 0:
                mode = None
                continue
            if line[0] in ('#','!'):
                continue
            if data[0] == 'charges':
                if not data[1] in self.__charge_types__:
                    raise ValueError('Unknown charge type: '+data[1])
                mode = 'charges:'+data[1]
                self.coul_params[data[0]] = data[1]
                continue
            elif data[0] == 'fbvs':
                self.coul_params[data[0]] = float(data[1])
                mode = None
                continue
            elif data[0] == 'interactions':
                mode = data[0]
                self.coul_params[data[0]] = []
                continue
            elif data[0] == 'rad_screened_cut':
                csp = data[1]
                rad = float(data[2])
                rads[csp] = rad
                if csp not in self.coul_params.keys():
                    self.coul_params[csp] = {}
                self.coul_params[csp]['rad'] = rad
                mode = None
                continue
            elif data[0] == 'terms':
                mode = data[0]
                if not data[1] in self.__term_types__:
                    raise ValueError('Unknown term type: '+data[1])
                self.coul_params[data[0]] = data[1]
            elif 'charges' in mode:
                if '_bvs' in mode:
                    if len(data) != 4:
                        raise ValueError('format of {0:s} seems wrong.'.format(infname))
                    csp = data[0]
                    if csp not in self.coul_params.keys():
                        self.coul_params[csp] = {}
                    vid = float(data[1])
                    rad = float(data[2])
                    npq = int(data[3])
                    self.coul_params[csp]['vid'] = vid
                    self.coul_params[csp]['rad'] = rad
                    self.coul_params[csp]['npq'] = npq
                    # vids[csp] = vid
                    # rads[csp] = rad
                    # npqs[csp] = npq
                elif 'fixed' in mode:
                    if len(data) != 2:
                        raise ValueError('format of {0:s} seems wrong.'.format(infname))
                    csp = data[0]
                    if csp not in self.coul_params.keys():
                        self.coul_params[csp] = {}
                    vid = float(data[1])
                    self.coul_params[csp]['vid'] = vid
                    # vids[csp] = vid
            elif 'interactions' in mode:
                csp1 = data[0]
                csp2 = data[1]
                self.coul_params[mode].append((csp1,csp2))
            else:
                pass
        return None

    def write(self,fname='in.params.Coulomb'):
        self.check()

        with open(fname,'w') as f:
            f.write('# in.params.Coulomb generated by inparams.py at '
                    +datetime.now().strftime('%y-%m-%d')+'\n')
            f.write('#\n')
            f.write('terms     '+self.coul_params['terms']+'\n')
            charges = self.coul_params['charges']
            f.write('charges   '+charges+'\n')
            if charges == 'fixed_bvs':
                for s in self.specorder:
                    vid = self.coul_params[s]['vid']
                    rad = self.coul_params[s]['rad']
                    npq = self.coul_params[s]['npq']
                    f.write('  {0:3s}  {1:4.1f}  {2:7.3f}  {3:2d}\n'.format(s,vid,rad,npq))
            elif charges == 'fixed':
                for s in specorder:
                    vid = self.coul_params[s]['vid']
                    f.write('  {0:3s}  {1:7.3f}\n'.format(s,vid))
            elif charges in ('variable', 'qeq'):
                for s in specorder:
                    chi = self.coul_params[s]['chi']
                    jii = self.coul_params[s]['jii']
                    e0 = self.coul_params[s]['e0']
                    qlow = self.coul_params[s]['qlow']
                    qup = self.coul_params[s]['qup']
                    f.write('  {0:3s}  {1:12.4e}  {2:12.4e}  {3:7.3f}  {4:7.3f}  {5:7.3f}\n'.format(s,chi,jii,e0,qlow,qup))
            else:
                raise ValueError('Unknown charges.')

            f.write('\n')
            if '_bvs' in self.coul_params['charges']:
                f.write('fbvs  {0:7.3f}\n'.format(self.coul_params['fbvs']))
            f.write(\n)
            for i in range(len(specorder)):
                si = specorder[i]
                for j in range(i,len(specorder)):
                    sj = specorder[j]
                    if not same_pair_exists([si,sj],pairs):
                        f.write('   {0:3s}  {1:3s}\n'.format(si,sj))
                
