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
import yaml
import math

_constfname = 'in.const.NN'
_paramfname = 'in.params.NN'
_logfname = 'log.make_params_NN'

#...initial range of parameters
_pmin = -0.01
_pmax =  0.01

_type_avail = ('Gaussian','cosine','polynomial','Morse','angular')
_type_2body = ('Gaussian','cosine','polynomial','Morse',)
_type_3body = ('angular')
_type2num = {
    'Gaussian': 1,
    'cosine': 2,
    'polynomial': 3,
    'Morse': 4,
    'angular': 101,
}
_nparam_type = {
    'Gaussian': 2,
    'cosine': 1,
    'polynomial': 1,
    'Morse': 3,
    'angular': 2,
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
    def __init__(self,isp1,isp2):
        self.isp1 = isp1
        self.isp2 = isp2
        self.symfuncs = []
        self.sfparams = []

    def __repr__(self):
        return "Pair({0:d},{1:d})".format(self.isp1,self.isp2)

    def get_isps(self):
        return self.isp1,self.isp2

    def set_symfunc(self,symfunc):
        self.symfuncs.append(symfunc)

    def set_sfparam(self,sfparam):
        self.sfparams.append(sfparam)

class Triplet(object):
    def __init__(self,isp1,isp2,isp3):
        self.isp1 = isp1
        self.isp2 = isp2
        self.isp3 = isp3
        self.symfuncs = []
        self.sfparams = []

    def __repr__(self):
        return "Triplet({0:d},{1:d},{2:d})".format(self.isp1,self.isp2,self.isp3)

    def get_isps(self):
        return self.isp1,self.isp2,self.isp3

    def set_symfunc(self,symfunc):
        self.symfuncs.append(symfunc)

    def set_sfparam(self,sfparam):
        self.sfparams.append(sfparam)

def ncomb(n,m):
    '''
    Calculate nCm.
    '''
    return math.factorial(n)/math.factorial(m)

def get_pairs(nsp):
    pairs= []
    for i in range(1,nsp+1):
        for j in range(i,nsp+1):
            pair = Pair(i,j)
            pairs.append(pair)
    return pairs

def get_triplets(nsp):
    pairs= []
    for i in range(1,nsp+1):
        for j in range(i,nsp+1):
            pairs.append([i,j])
    #return pairs
    triplets= []
    for i in range(1,nsp+1):
        tmp= [i]
        for pair in pairs:
            # if _selected_triplet and \
            #    not (tmp[0],pair[0],pair[1]) in _selected_triplet:
            #     continue
            triplet = Triplet(tmp[0],pair[0],pair[1])
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

def get_params(sfname,isp1,isp2,isp3=None):
    if sfname == 'Gaussian':
        rs,etas = get_gauss_params(isp1,isp2)
        return rs,etas
    elif sfname == 'cosine':
        rk = get_cosine_params(isp1,isp2)
        return rk
    elif sfname == 'polynomial':
        a1s = get_polynomial_params(isp1,isp2)
        return a1s
    elif sfname == 'Morse':
        ds,alps,rs = get_Morse_params(isp1,isp2)
        return ds,alps,rs
    elif sfname == 'angular':
        if not isp3:
            raise ValueError('isp3 is required for 3-body term.')
        ang = get_angular_params(isp1,isp2,isp3)
        return ang
    else:
        raise ValueError('No such symmetry function: '+sfname)

def get_gauss_params(isp1,isp2):
    print('\nDetermine Gaussian parameters for {0:d}-{1:d}:'.format(isp1,isp2))
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
    print('Gaussian width parameters [1/Ang]'
          +' (several values separated by white-space):')
    etas = [ float(eta) for eta in sys.stdin.readline().split() ]
    print('Number of Gaussian width parameters = {0:d}'.format(len(etas)))
    return rs, etas

def get_cosine_params(isp1,isp2):
    rk = []
    return rk

def get_polynomial_params(isp1,isp2):
    a1s = []
    return a1s

def get_Morse_params(isp1,isp2):
    ds = []
    alps = []
    rs = []
    return ds,alps,rs

def get_angular_params(isp1,isp2,isp3):
    print('\nDetermine angular parameters for {0:d}-{1:d}-{2:d}:'.format(isp1,isp2,isp3))
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
        for it,t in enumerate(pair.symfuncs):
            n = 1
            if t == 'Gaussian':
                for sfps in pair.sfparams[it]:
                    if isinstance(sfps,list) or isinstance(sfps,tuple):
                        n *= len(sfps)
            elif t == 'polynomial' or 'cosine':
                n = len(pair.sfparams[it])

            nsf2 += n
    return nsf2

def get_nsf3(triplets):
    nsf3 = 0
    if triplets is None:
        return nsf3
    for triplet in triplets:
        for it,t in enumerate(triplet.symfuncs):
            n = len(triplet.sfparams[it])
            nsf3 += n
    return nsf3
    
def create_param_files(inputs,nsf2,pairs,nsf3,triplets):

    # print('inputs = ',inputs)
    nl = inputs['nl']
    nsp = inputs['nsp']
    nhl = inputs['nhl']

    with open(_constfname,'w') as f:
        nhl[0] = nsf2 +nsf3

        f.write(' {0:5d} {1:5d}'.format(nl,nsp))
        # print('nl = ',nl)
        # print('nhl= ',nhl)
        for il in range(nl+1):
            f.write(' {0:5d}'.format(nhl[il]))
        f.write('\n')

        if nsf2 > 0:
            for pair in pairs:
                ia,ja = pair.get_isps()
                for isf,sf in enumerate(pair.symfuncs):
                    if sf == 'Gaussian':
                        rs,etas = pair.sfparams[isf]
                        for eta in etas:
                            for r in rs:
                                f.write(' {0:3d}'.format(_type2num[sf]) \
                                        +' {0:3d} {1:3d}'.format(ia,ja) \
                                        +' {0:10.4f} {1:10.4f}\n'.format(eta,r))
                    elif sf == 'cosine':
                        rk = pair.sfparams[isf]
                        for r in rk:
                            f.write(' {0:3d}'.format(_type2num[sf]) \
                                    +' {0:3d} {1:3d}'.format(ia,ja) \
                                    +' {0:10.4f}\n'.format(r))
                    elif sf == 'polynomial':
                        a1s = pair.sfparams[isf]
                        for a1 in a1s:
                            f.write(' {0:3d}'.format(_type2num[sf]) \
                                    +' {0:3d} {1:3d}'.format(ia,ja) \
                                    +' {0:10.4f}\n'.format(a1))
                    elif sf == 'Morse':
                        ds,alps,rs = pair.sfparams[isf]
                        for d in ds:
                            for alp in alps:
                                for r in rs:
                                    f.write(' {0:3d}'.format(_type2num[sf]) \
                                            +' {0:3d} {1:3d}'.format(ia,ja) \
                                            +' {0:10.4f} {1:10.4f}'.format(d,alp) \
                                            +' {0:10.4f}\n'.format(r))
        if nsf3 > 0:
            for triple in triplets:
                ia,ja,ka = triple.get_isps()
                for isf,sf in enumerate(triple.symfuncs):
                    angs = [ -math.cos(a/180*math.pi) for a in triple.sfparams[isf] ]
                    for ang in angs:
                        f.write(' {0:3d}'.format(_type2num[sf]) \
                                +' {0:3d} {1:3d} {2:3d}'.format(ia,ja,ka) \
                                +' {0:10.4f}\n'.format(ang))
        f.close()

    rc2 = inputs['rc2']
    rc3 = inputs['rc3']
    with open(_paramfname,'w') as g:
        if nl == 1:
            nc= nhl[0]*nhl[1] +nhl[1]
        elif nl == 2:
            nc= nhl[0]*nhl[1] +nhl[1]*nhl[2] +nhl[2]
        g.write(' {0:6d} {1:10.3f} {2:10.3f}\n'.format(nc,rc2,rc3))
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
        f.write(yaml.dump(conf))

def load_config(fname):
    with open(fname,'r') as f:
        dic = yaml.load(f)
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
                isp1,isp2 = pair.get_isps()
                for t in inputs['type_use']:
                    if t in _type_2body:
                        pair.set_symfunc(t)
                        sfparam = get_params(t,isp1,isp2)
                        pair.set_sfparam(sfparam)
            
        if b3_exists:
            for triplet in triplets:
                isp1,isp2,isp3 = triplet.get_isps()
                for t in inputs['type_use']:
                    if t in _type_3body:
                        triplet.set_symfunc(t)
                        sfparam = get_params(t,isp1,isp2,isp3)
                        triplet.set_sfparam(sfparam)
    nsf2 = get_nsf2(pairs)
    nsf3 = get_nsf3(triplets)

    print('Number of pairs, triplets = ',nsf2,nsf3)
    
    create_param_files(inputs,nsf2,pairs,nsf3,triplets)
    if not load:
        save_config('out.conf.make_params_NN',inputs,pairs,triplets)
    print('\nCheck '+_constfname+' and '+_paramfname)
