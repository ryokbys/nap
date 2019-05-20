#!/usr/bin/env python
"""
Make in.const.NN and in.params.NN to be used in NN potential.
It is a new version of `make_const-params_NN.py` to be fully replaced.
If you want to specify precisely different symmetry functions for each pairs,
you may had better modify config file and load it.

Usage:
  make_params_NN.py [options]
  make_params_NN.py load [options] CONFIG

Options:
  -h, --help  Show this help message and exit.
"""
from __future__ import print_function

import sys
from docopt import docopt
import random,math
import json

_descfname = 'in.params.desc'
_paramfname = 'in.params.NN2'
_logfname = 'log.make_params_NN'

#...initial range of parameters
_pmin = -20.0
_pmax =  20.0

_type_avail = ('Gaussian','cosine','polynomial','Morse','angular',
               'angular1','angular2','angular3','angular4','angular5')
_type_2body = ('Gaussian','cosine','polynomial','Morse',)
_type_3body = ('angular','angular1','angular2','angular3','angular4','angular5')
_type2num = {
    'Gaussian': 1,
    'cosine': 2,
    'polynomial': 3,
    'Morse': 4,
    'angular': 101,
    'angular1': 101,
    'angular2': 102,
    'angular3': 103,
    'angular4': 104,
    'angular5': 105,
}
_nparam_type = {
    'Gaussian': 2,
    'cosine': 1,
    'polynomial': 1,
    'Morse': 3,
    'angular': 2,
    'angular1': 2,
    'angular2': 3,
    'angular3': 2,
    'angular4': 2,
    'angular5': 3,
}

def _num2type(num):
    for k,v in _type2num.items():
        if v == num:
            return k
    return None

# _selected_triplet = [
#     (3, 1, 1), (3, 1, 2), (3, 1, 3), (3, 2, 2), (3, 2, 3), (3, 3, 3),
# ]


class Pair(object):
    def __init__(self,csp1,csp2):
        self.csp1 = csp1
        self.csp2 = csp2
        self.symfuncs = []
        self.sfparams = []

    def __repr__(self):
        return "Pair({0:s},{1:s})".format(self.csp1,self.csp2)

    def get_csps(self):
        return self.csp1,self.csp2

    def set_symfunc(self,symfunc):
        self.symfuncs.append(symfunc)

    def set_sfparam(self,sfparam):
        self.sfparams.append(sfparam)

class Triplet(object):
    def __init__(self,csp1,csp2,csp3):
        self.csp1 = csp1
        self.csp2 = csp2
        self.csp3 = csp3
        self.symfuncs = []
        self.sfparams = []

    def __repr__(self):
        return "Triplet({0:s},{1:s},{2:s})".format(self.csp1,self.csp2,self.csp3)

    def get_csps(self):
        return self.csp1,self.csp2,self.csp3

    def set_symfunc(self,symfunc):
        self.symfuncs.append(symfunc)

    def set_sfparam(self,sfparam):
        self.sfparams.append(sfparam)

def ncomb(n,m):
    '''
    Calculate nCm = n!/m!/(n-m)!.
    '''
    d = max(n-m,1)
    return math.factorial(n)/math.factorial(m)/ math.factorial(d)

def get_pairs(specorder=[]):
    if len(specorder) < 1:
        raise ValueError('Specorder has to be specified.')
    pairs= []
    for i in range(len(specorder)):
        cspi = specorder[i]
        for j in range(i,len(specorder)):
            cspj = specorder[j]
            pair = Pair(cspi,cspj)
            pairs.append(pair)
    return pairs

def get_triplets(specorder=[]):
    if len(specorder) < 1:
        raise ValueError('Specorder has to be specified.')
    pairs= []
    for i in range(len(specorder)):
        cspi = specorder[i]
        for j in range(len(specorder)):
            cspj = specorder[j]
            pairs.append([cspi,cspj])
    #return pairs
    triplets= []
    for i in range(len(specorder)):
        tmp= specorder[i]
        for pair in pairs:
            # if _selected_triplet and \
            #    not (tmp[0],pair[0],pair[1]) in _selected_triplet:
            #     continue
            triplet = Triplet(tmp,pair[0],pair[1])
            triplets.append(triplet)
            #triplets.append(tmp+pair)
    return triplets
    
def get_inputs():

    inputs = {}
    print('Cutoff radius for 2-body terms: ')
    inputs['rc2'] = float(sys.stdin.readline())

    print('Cutoff radius for 3-body terms: ')
    inputs['rc3'] = float(sys.stdin.readline())

    print('Number of species: ')
    inputs['nsp'] = int(sys.stdin.readline())

    print('Number of hidden layers (1 or 2): ')
    inputs['nl'] = int(sys.stdin.readline())
    if not inputs['nl'] in (1,2):
        raise ValueError(' nl is not 1 or 2, nl={0:d}'.format(inputs['nl']))

    nhl= [0,]
    for il in range(1,inputs['nl']+1):
        print('Number of nodes in a layer-{0:d}: '.format(il))
        nhl.append(int(sys.stdin.readline()))
    inputs['nhl'] = nhl
    print('nhl = ',nhl)

    #...Ask which symmetry functions are used
    print('Choose the symmetry functions to be used from the following:')
    for t in _type_avail:
        print(' {0:3d}:{1:s}'.format(_type2num[t],t))
    print('Space separated format, for example, 1 2 101')
    type_use = [ _num2type(int(i)) for i in sys.stdin.readline().split() ]
    inputs['type_use'] = type_use
    print('type_use =',type_use)
    print('Selected symmetry functions:')
    for t in type_use:
        print(' - '+t)
    
    #...check types
    if len(type_use) < 1:
        raise ValueError('No type specified...')
    for t in type_use:
        if not t in _type_avail:
            raise ValueError('No such available type: '+t)
    
    return inputs

def get_params(sfname,csp1,csp2,csp3=None):
    if sfname == 'Gaussian':
        rs,etas = get_gauss_params(csp1,csp2)
        return rs,etas
    elif sfname == 'cosine':
        rk = get_cosine_params(csp1,csp2)
        return rk
    elif sfname == 'polynomial':
        a1s = get_polynomial_params(csp1,csp2)
        return a1s
    elif sfname == 'Morse':
        ds,alps,rs = get_Morse_params(csp1,csp2)
        return ds,alps,rs
    elif sfname == 'angular':
        if not csp3:
            raise ValueError('csp3 is required for 3-body term.')
        ang = get_angular_params(csp1,csp2,csp3)
        return ang
    else:
        raise ValueError('No such symmetry function: '+sfname)

def get_gauss_params(csp1,csp2):
    print('\nDetermine Gaussian parameters for {0:s}-{1:s}:'.format(csp1,csp2))
    print('Minimum interaction distance [Ang]:')
    rmin = float(sys.stdin.readline())
    print('Maximum interaction distance [Ang]:')
    rmax = float(sys.stdin.readline())
    print('Number of Gaussian peak positions'
          +' bewteen {0:8.3f} - {1:8.3f}:'.format(rmin,rmax))
    nrs = int(sys.stdin.readline())
    if nrs == 1:
        rs = [0.0]
    elif nrs > 1:
        dr = (rmax-rmin) /(nrs-1)
        rs = [ rmin +dr*i for i in range(nrs) ]
    else:
        raise ValueError('Something wrong. nrs={0:d}'.format(nrs))
    print('Gaussian peaks are : ')
    for r in rs:
        print(' {0:8.3f}'.format(r),end='')
    print('')
    print('Gaussian width parameters [1/Ang]'
          +' (several values separated by white-space):')
    etas = [ float(eta) for eta in sys.stdin.readline().split() ]
    print('Number of Gaussian width parameters = {0:d}'.format(len(etas)))
    return rs, etas

def get_cosine_params(csp1,csp2):
    rk = []
    return rk

def get_polynomial_params(csp1,csp2):
    a1s = []
    return a1s

def get_Morse_params(csp1,csp2):
    ds = []
    alps = []
    rs = []
    return ds,alps,rs

def get_angular_params(csp1,csp2,csp3):
    print('\nDetermine angular parameters for {0:s}-{1:s}-{2:s}:'.format(csp1,csp2,csp3))
    print('Special angles [degree] (space-separation):')
    # angs = [ math.cos(float(a)/180*math.pi) for a in sys.stdin.readline().split() ]
    angs = [ float(a) for a in sys.stdin.readline().split() ]
    print('Number of anguler parameters = {0:d}'.format(len(angs)))
    return angs

def get_nsf2(pairs):
    nsf2 = 0
    if pairs is None:
        return nsf2
    for pair in pairs:
        nsf2 += len(pair.sfparams)
    return nsf2

def get_nsf3(triplets):
    nsf3 = 0
    if triplets is None:
        return nsf3
    for triplet in triplets:
        nsf3 += len(triplet.sfparams)
    return nsf3
    
def create_param_files(inputs,nsf2,pairs,nsf3,triplets):

    # print('inputs = ',inputs)
    nl = inputs['nl']
    nsp = inputs['nsp']
    nhl = inputs['nhl']
    rc2 = inputs['rc2']
    rc3 = inputs['rc3']
    
    with open(_descfname,'w') as f:
        nsf = nsf2 + nsf3
        nhl[0] = nsf2 +nsf3

        #f.write(' {0:5d} {1:5d}'.format(nl,nsp))
        f.write(' {0:5d} {1:5d}\n'.format(nsp,nsf))
        # print('nl = ',nl)
        # print('nhl= ',nhl)
        # for il in range(nl+1):
        #     f.write(' {0:5d}'.format(nhl[il]))
        # f.write('\n')

        if nsf2 > 0:
            for pair in pairs:
                cspi,cspj = pair.get_csps()
                for isf,sf in enumerate(pair.sfparams):
                    t = sf[0]
                    if t == 'Gaussian':
                        if len(sf) != 3:
                            raise RuntimeError('Num of Gaussian params is wrong.')
                        eta,rs = sf[1],sf[2]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s}'.format(cspi,cspj) \
                                +' {0:6.2f}'.format(rc2) \
                                +' {0:10.4f} {1:10.4f}\n'.format(eta,rs))
                    elif t == 'cosine':
                        if len(sf) != 2:
                            raise RuntimeError('Num of cosine params is wrong.')
                        r = sf[1]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s}'.format(cspi,cspj) \
                                +' {0:6.2f}'.format(rc2) \
                                +' {0:10.4f}\n'.format(r))
                    elif t == 'polynomial':
                        if len(sf) != 2:
                            raise RuntimeError('Num of polynomial params is wrong.')
                        a1 = sf[1]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s}'.format(cspi,cspj) \
                                +' {0:6.2f}'.format(rc2) \
                                +' {0:10.4f}\n'.format(a1))
                    elif t == 'Morse':
                        if len(sf) != 4:
                            raise RuntimeError('Num of Morse params is wrong.')
                        d,alp,r = sf[1],sf[2],sf[3]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s}'.format(cspi,cspj) \
                                +' {0:6.2f}'.format(rc2) \
                                +' {0:10.4f} {1:10.4f}'.format(d,alp) \
                                +' {0:10.4f}\n'.format(r))
        if nsf3 > 0:
            for triple in triplets:
                cspi,cspj,cspk = triple.get_csps()
                for isf,sf in enumerate(triple.sfparams):
                    t = sf[0]
                    if t in ('angular','angular1'):
                        a = sf[1]
                        ang = -math.cos(a/180*math.pi)
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s} {2:<3s}'.format(cspi,cspj,cspk) \
                                +' {0:6.2f}'.format(rc3) \
                                +' {0:10.4f}\n'.format(ang))
                    elif t in ('angular2'):
                        a = sf[1]
                        ang = -math.cos(a/180*math.pi)
                        zeta = sf[2]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s} {2:<3s}'.format(cspi,cspj,cspk) \
                                +' {0:6.2f}'.format(rc3) \
                                +' {0:10.4f}'.format(ang) \
                                +' {0:10.4f}\n'.format(zeta))
                    elif t in ('angular3','angular4'):
                        a = sf[1]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s} {2:<3s}'.format(cspi,cspj,cspk) \
                                +' {0:6.2f}'.format(rc3) \
                                +' {0:10.4f}\n'.format(a))
                    elif t in ('angular5'):
                        eta = sf[1]
                        rs = sf[2]
                        f.write(' {0:3d} '.format(_type2num[t]) \
                                +' {0:<3s} {1:<3s} {2:<3s}'.format(cspi,cspj,cspk) \
                                +' {0:6.2f}'.format(rc3) \
                                +' {0:10.4f} {1:10.4f}\n'.format(eta,rs))
        f.close()

    with open(_paramfname,'w') as g:
        if nl == 1:
            nc= nhl[0]*nhl[1] +nhl[1]
            g.write(' {0:4d} {1:6d} {2:4d}\n'.format(nl,nsf,nhl[1]))
        elif nl == 2:
            nc= nhl[0]*nhl[1] +nhl[1]*nhl[2] +nhl[2]
            g.write(' {0:4d} {1:4d} {2:3d} {3:3d}\n'.format(nl,nsf,nhl[1],nhl[2]))
        for ic in range(nc):
            g.write(' {0:10.6f}'.format(random.uniform(_pmin,_pmax)))
            g.write(' {0:10.4f} {1:10.4f}\n'.format(_pmin,_pmax))
        g.close()

    return

def save_config(fname,inputs,pairs,triplets):
    conf = {}
    conf['inputs'] = inputs
    conf['pairs'] = pairs
    conf['triplets'] = triplets
    with open(fname,'w') as f:
        f.write(json.dump(conf))

def load_config(fname):
    with open(fname,'r') as f:
        dic = json.load(f)
    inputs = dic['inputs']
    pairs = dic['pairs']
    triplets = dic['triplets']
    return inputs, pairs, triplets


#========================================================= main routine
if __name__ == "__main__":

    args = docopt(__doc__)
    load = args['load']

    if load:
        inconfname = args['CONFIG']
        print('Loading {0:s}...'.format(inconfname))
        inputs,pairs,triplets = load_config(inconfname)
    else:
        inputs = get_inputs()
    
        b2_exists = False
        b3_exists = False
        for t in inputs['type_use']:
            if t in _type_2body:
                b2_exists = True
            if t in _type_3body:
                b3_exists = True
    
        nsp = inputs['nsp']
        pairs = None
        if b2_exists:
            ncmb2= nsp+ ncomb(nsp,2)
            pairs = get_pairs(nsp)
            print('Number of pair combinations: {0:d}'.format(ncmb2))
            print('Pairs=',pairs)
            if len(pairs) != ncmb2:
                print('len(pairs),ncmb2=',len(pairs),ncmb2)
                raise ValueError('len(pairs) != ncmb2')
    
        triplets = None
        if b3_exists:
            if not b2_exists:
                ncmb2= nsp+ ncomb(nsp,2)/2
            ncmb3= ncmb2*nsp
            print('Number of triplets: {0:d}'.format(ncmb3))
            triplets= get_triplets(nsp)
            print('Triplets=',triplets)
            if len(triplets) != ncmb3:
                print('Since len(triplets) != ncmb3, set ncmb3 = len(triplets)')
                ncmb3 = len(triplets)

        if b2_exists:
            for pair in pairs:
                csp1,csp2 = pair.get_csps()
                for t in inputs['type_use']:
                    if t in _type_2body:
                        pair.set_symfunc(t)
                        sfparam = get_params(t,csp1,csp2)
                        pair.set_sfparam(sfparam)
            
        if b3_exists:
            for triplet in triplets:
                csp1,csp2,csp3 = triplet.get_csps()
                for t in inputs['type_use']:
                    if t in _type_3body:
                        triplet.set_symfunc(t)
                        sfparam = get_params(t,csp1,csp2,csp3)
                        triplet.set_sfparam(sfparam)
    nsf2 = get_nsf2(pairs)
    nsf3 = get_nsf3(triplets)

    print('Number of pairs, triplets = ',nsf2,nsf3)
    
    create_param_files(inputs,nsf2,pairs,nsf3,triplets)
    if not load:
        save_config('out.conf.make_params_NN',inputs,pairs,triplets)
    print('\nCheck '+_descfname+' and '+_paramfname)
