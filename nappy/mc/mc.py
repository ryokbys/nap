#!/bin/local/env python
"""
Monte Carlo (MC) simulation of small precipitate 
in Al-Mg-Si system with using pmd.

Usage:
  mc.py [options]

Options:
  -h,--help  Show this help message and exit.
  -n NSTEPS  Number of MC steps. [default: 10]
  --nMg=NMG  Number of Mg in the system. [default: 1]
  --nSi=NSI  Number of Si in the system. [default: 1]
  --nVac=NVAC  Number of Vac in the system. [default: 1]
  --kinetic  Perform kinetic MC instead of Metropolis MC. [default: False]
  --temperature=TEMPERATURE
             Temperature for MC simulation. [default: 300.0]
  --size=SIZE
             Size of the system. [default: 3x3x3]
  --latt-const=LATTCONST
             Lattice constant of the system. [default: 4.0448]
  --ase-relax
             Use ASE optimization instead of relaxing within MD code. [default: False]
  --init-clustered
             Whether the initial structure clustered. [default: False]
  --restart
             Restart from json file. [default: False]
  --parallel=PARALLEL
             Number of parallel nodes along x,y,z used for pmd calc. [default: 1x1x1]
  --relax-steps=NRSTEPS
             Steps for relaxation of atomic positions. 
             If it is negative, the num of atoms applies. [default: -1]
  --movable=MOVABLE
             List of movable species. [default: Vac]
"""
from __future__ import print_function

import os,sys,copy
from docopt import docopt
import numpy as np
import random
import json

from ase.units import kB
from ase.io import read, write

HOME = os.environ['HOME']
_fsymbols = 'dat.symbols'
_ferg = 'dat.erg'

sys.path.append(HOME+'/src/nap/nappy/')
from interface.ase.pmdrun import PMD

sys.path.append(os.path.dirname(__file__))
from pairlist import make_pair_list

__author__ = "Ryo KOBAYASHI"
__version__ = "160507"


def make_fcc_Al(lattconst=4.04159,size=(3,3,3),):
    """
    Make fcc Al structure.
    """
    from ase.lattice.cubic import FaceCenteredCubic
    #...create Al fcc system
    atoms=FaceCenteredCubic(symbol='Al',
                            latticeconstant=lattconst,
                            size=size)
    return atoms

def make_symbols_random(sites,symbols0,nlspr,lspr,
                        nsolute0={}):
    """
    Make initial distribution of Mg, Si, and Vac clustered at the center.
    """
    symbols = copy.deepcopy(symbols0)
    nsolute = copy.deepcopy(nsolute0)

    #...make random array of symbols to be replaced
    arr = []
    for k,v in nsolute.items():
        for i in range(v):
            arr.append(k)
    random.shuffle(arr)

    idcs = [ i for i in range(len(symbols)) ]
    random.shuffle(idcs)

    for i in range(len(arr)-1,0-1,-1):
        idx = idcs.pop(i)
        symbols[idx] = arr.pop(i)
    return symbols


def make_symbols_clustered(sites,symbols0,nlspr,lspr,
                           nsolute0={}):
    """
    Make initial distribution of Mg, Si, and Vac clustered at the center.
    """
    symbols = copy.deepcopy(symbols0)
    nsolute = copy.deepcopy(nsolute0)

    #...replace 1st center atom
    cntr = (0.5, 0.5, 0.5)
    icntr = -1
    dmin = 1.e+30
    for i in range(len(sites)):
        d = (sites[i,0]-cntr[0])**2 \
            +(sites[i,1]-cntr[1])**2 \
            +(sites[i,2]-cntr[2])**2
        if d < dmin:
            dmin = d
            icntr = i

    #...make random array of symbols to be replaced
    arr = []
    for k,v in nsolute.items():
        for i in range(v):
            arr.append(k)
    random.shuffle(arr)

    symbols[icntr] = arr[0]
    arr.pop(0)

    # nsol = len(nsolute.keys())
    # while True:
    #     symbols[icnt
    #     irep = np.random.randint(0,nsol)
    #     if not nsolute.keys()[irep] == 0:
    #         symbols[icntr] = nsolute.keys()[irep]
    #         nsolute[symbols[icntr]] -= 1
    #         break

    while True:
        for i in range(len(symbols)):
            if symbols[i] == 'Al': continue
            for jj in range(nlspr[i]):
                j = lspr[i,jj]
                if symbols[j] == 'Al':
                    symbols[j] = arr.pop(0)
                    if len(arr) == 0:
                        print('initial symbols:')
                        print(symbols)
                        return symbols
                    # jrep = np.random.randint(0,nsol)
                    # if not nsolute.keys()[jrep] == 0:
                    #     symbols[j] = nsolute.keys()[jrep]
                    #     nsolute[symbols[j]] -= 1
                    # if count_nsol(nsolute) == 0:
                    #     print('initial symbols:')
                    #     print(symbols)
                    #     return symbols
    raise RuntimeError('Something wrong.')
    

def count_nsol(nsolute):
    n = 0
    for k,v in nsolute.items():
        n += v
    return n

def moved(symbols1,symbols0,species):
    """
    Return if the specified species moved
    """
    for i,s1 in enumerate(symbols1):
        if s1 == symbols0[i]: continue
        if s1 in species:
            return True
    return False

def atoms_from_symbols(atoms0,symbols):
    atoms = atoms0.copy()
    for i,s in enumerate(symbols):
        if s in ('Mg','Si'):
            atoms[i].symbol = s
    npopped = 0
    for i,s in enumerate(symbols):
        if s == 'Vac':
            atoms.pop(i-npopped)
            npopped += 1
    return atoms

def choose(symbols,solutes):
    """
    Choose one atom to be replaced or exchanged.
    """
    #...count num of candidates
    cnt = 0
    for s in symbols:
        if s in solutes:
            cnt += 1
    #...choose one randomly
    chosen = np.random.randint(0,cnt)
    cnt = 0
    for i,s in enumerate(symbols):
        if s in solutes:
            if cnt == chosen:
                return i
            cnt += 1
    raise RuntimeError('Something wrong...')


def get_prec_energy(atoms,symbols,chempots):
    erg = atoms.get_potential_energy()
    for s in symbols:
        erg -= chempots[s]
    return erg

def relax(atoms0,steps=100,pmd_relax=False):
    atoms = copy.deepcopy(atoms0)
    if pmd_relax:
        calc = atoms.get_calculator()
        calc.relax(atoms=atoms,nsteps=steps)
        relaxed_spos = calc.get_relaxed_scaled_positions()
        atoms.set_scaled_positions(relaxed_spos)
    else:
        from ase.optimize import QuasiNewton
        relax = QuasiNewton(atoms,logfile=None)
        relax.run(steps=steps)
    return atoms


def calc_chemical_potentials(atoms0,
                             solutes=(),
                             pmd_relax=False,
                             pmd_cmd='~/src/nap/pmd/pmd',
                             paranode=[1,1,1]):
    chempots = {}

    #...Al
    calc = PMD(label='pmd', command=pmd_cmd,
               force_type='NN', cutoff_radius=5.8,
               specorder=['Al','Mg','Si'],
               mass=[26.982,24.305,28.085],
               num_nodes_x=paranode[0],
               num_nodes_y=paranode[1],
               num_nodes_z=paranode[2])
    atoms0.set_calculator(calc)
    atoms = relax(atoms0,pmd_relax=pmd_relax,
                  steps=len(atoms0))
    chempots['Al'] = atoms.get_potential_energy()/len(atoms)

    symbols0 = atoms.get_chemical_symbols()
    for s in solutes:
        symbols = copy.deepcopy(symbols0)
        symbols[0] = s
        atoms = atoms_from_symbols(atoms0,symbols)
        atoms.set_calculator(calc)
        atoms = relax(atoms,pmd_relax=pmd_relax,
                      steps=len(atoms))
        chempots[s] = atoms.get_potential_energy() \
                      - chempots['Al']*symbols.count('Al')
    return chempots


def dumps_symbols(symbols):
    txt = ''
    for s in symbols:
        txt += s[0]
    txt += '\n'
    return txt

def loads_symbols(txt):
    """
    Input txt should include only array of characters of
    length of number of atoms.
    """
    symbols = []
    for t in txt:
        if t == 'A':
            s = 'Al'
        elif t == 'M':
            s = 'Mg'
        elif t == 'S':
            s = 'Si'
        elif t == 'V':
            s = 'Vac'
        else:
            raise ValueError('Unknow type: '+s)
        symbols.append(s)
    
    return symbols

def symbols_same(s0,s1):
    """
    Check if two symbols are the same.
    """
    for i in range(len(s0)):
        if s0[i] != s1[i]: return False
    return True

def symbols_exist(symbols,history):
    """
    Check if there is any identical symbols in the history.
    """
    for ih,symb in enumerate(history):
        if symb == symbols: return ih
    return -1

def diff_list(l1,l2):
    """
    Return the elements different from each other.
    """
    arr = []
    for i in range(len(l1)):
        if l1[i] != l2[i]:
            arr.append(l1[i])
            arr.append(l2[i])
            return arr
    return arr

def symbols_to_erg(atoms0,symbols,calc,chempots,
                   pmd_relax,nrstps):
    atoms = atoms_from_symbols(atoms0,symbols)
    atoms.set_calculator(calc)
    if nrstps == 0:
        erg = get_prec_energy(atoms,symbols,chempots)
    else:
        atoms = relax(atoms,pmd_relax=pmd_relax,
                      steps=nrstps)
        erg = get_prec_energy(atoms,symbols,chempots)
    return erg


def num_solutes(symbols,sol=None):
    n = 0
    for s in symbols:
        if sol is None:
            if s is not 'Al': n += 1
        else:
            if s is sol: n += 1
    return n

def check_num_solutes(symbols,prev_symbols):
    sols = ['Si','Mg','Vac']
    for s in sols:
        n1 = num_solutes(symbols,s)
        n0 = num_solutes(prev_symbols,s)
        if n1 != n0:
            print(' WARNING: num of '+s+' was changed from',n0,' to ',n1)

def run_lmc(atoms0,nsolute={},
            sites=[],
            symbols=[],
            nsteps=100,
            temperature=300.0,
            pmd_relax=False,
            pmd_cmd='~/src/nap/pmd/pmd',
            paranode=[1,1,1],
            relax_steps=-1,
            movable=[]):
    """
    Run lattice MC using pmd.

    Atoms0 object must be kept to original fcc Al system.
    To modify this use atoms_from_symbols() fucntion with
    modified symbols array.
    """
    #print('system: '+atoms.get_chemical_formula())

    #....chemical potential for species, may have to be determined by calculation?
    # chempots={'Al':-3.54703,
    #           'Mg':-1.42109,
    #           'Si':-4.18928}
    chempots = calc_chemical_potentials(atoms0,
                                        solutes=nsolute.keys(),
                                        pmd_relax=pmd_relax,
                                        pmd_cmd=pmd_cmd,
                                        paranode=paranode)
    print('chempots:')
    for k,v in chempots.items():
        print('  {0:s}:  {1:12.5f}'.format(k,v))

    if relax_steps < 0:
        nrstps = len(atoms0)
    else:
        nrstps = relax_steps
    num_adopted = 0
    num_keep = 0
    for i in range(len(symbols)):
        if symbols[i] in nsolute.keys():
            num_keep += 1
    print('total num of Mg,Si,Vac to be kept = ',num_keep)

    # modes = ('exchange','replace')
    modes = ('exchange',)
    sys.stdout.flush()

    calc = PMD(label='pmd', command=pmd_cmd,
               force_type='NN', cutoff_radius=5.8,
               specorder=['Al','Mg','Si'],
               mass=[26.982,24.305,28.085],
               num_nodes_x=paranode[0],
               num_nodes_y=paranode[1],
               num_nodes_z=paranode[2])
    atoms = atoms_from_symbols(atoms0,symbols)
    atoms.set_calculator(calc)
    if nrstps == 0:
        erg0 = get_prec_energy(atoms,symbols,chempots)
    else:
        atoms = relax(atoms,pmd_relax=pmd_relax,
                      steps=nrstps)
        erg0 = get_prec_energy(atoms,symbols,chempots)
    ergp = erg0

    outerg = open(_ferg,'w')
    outsmbl = open(_fsymbols,'w')

    txt = ' {0:8d} {1:15.7e} {2:s}'.format(0,erg0,atoms.get_chemical_formula())
    print('istp,erg='+txt)
    outerg.write(txt+'\n')
    outsmbl.write('{0:8d}  '.format(0)+dumps_symbols(symbols))
    write('POSCAR_{0:05d}'.format(num_adopted),images=atoms,
          format='vasp',direct=True,vasp5=True,sort=True)
    sys.stdout.flush()
    istp = 1
    adopted = False
    while istp <= nsteps:
        ichosen = choose(symbols,movable)
        #ichosen = choose(symbols,nsolute.keys())
        #ichosen = choose(symbols,('Vac',))
        ispcs = symbols[ichosen]
        prev_symbols = copy.deepcopy(symbols)
        irnd = np.random.randint(0,len(modes))
        mode = modes[irnd]
        # print(' irnd,mode = ',irnd,mode)
        if mode == 'exchange':
            nn = nlspr[ichosen]
            #...check if there are atoms with different symbols in lspr
            ncandidate = 0
            for inn in range(nn):
                jc = lspr[ichosen,inn]
                js = symbols[jc]
                if js != ispcs:
                    ncandidate += 1
            if ncandidate == 0: continue
            #...choose one of the neighbors to be exchaned
            icandidate = np.random.randint(0,ncandidate)
            ncandidate = 0
            for inn in range(nn):
                jc = lspr[ichosen,inn]
                js = symbols[jc]
                if js != ispcs:
                    if icandidate == ncandidate:
                        jchosen = jc
                        break
                    ncandidate += 1
            jspcs = symbols[jchosen]
            #...exchange symbols
            symbols[ichosen] = jspcs
            symbols[jchosen] = ispcs
            #atoms.set_chemical_symbols(symbols)
            # print('{0:s} '.format(mode)
            #       +'{0:5d}({1:s})<=>{2:5d}({3:s})'.format(ichosen,ispcs,
            #                                               jchosen,jspcs))

        elif mode == 'replace':
            ncandidate = 0
            for js in replace_spcs:
                if js == ispcs: continue
                ncandidate += 1
            icandidate = np.random.randint(0,ncandidate)
            ncandidate = 0
            for js in replace_spcs:
                if js == ispcs: continue
                if ncandidate == icandidate:
                    jspcs = js
                    break
            #...replace symbol of i
            symbols[ichosen] = jspcs
            # print('{0:s} '.format(mode)
            #       +'{0:5d} {1:s}==>{2:s} '.format(ichosen,ispcs,jspcs))

        else:
            raise RuntimeError('There is no mode like '+mode)
        atoms = atoms_from_symbols(atoms0,symbols)
        #atoms.set_chemical_symbols(symbols)
        # print(atoms.get_chemical_formula())
        atoms.set_calculator(calc)
        if nrstps == 0:
            erg = get_prec_energy(atoms,symbols,chempots)
        else:
            atoms = relax(atoms,pmd_relax=pmd_relax,
                          steps=nrstps)
            erg = get_prec_energy(atoms,symbols,chempots)
        prob = np.exp(-(erg-ergp)/(temperature*kB))
        rand = np.random.random()
        if prob > rand:
            adopted = True
            num_adopted += 1
            txt_adopted = '  adopted: N={0:6d}, '.format(num_adopted) \
                          +'DE={0:12.4e}, P={1:11.3e}'.format((erg-ergp),prob)
            ergp = erg
            #...output new structure
            outsmbl.write('{0:8d}  '.format(num_adopted)+dumps_symbols(symbols))
            if moved(symbols,prev_symbols,['Mg','Si']):
                write('POSCAR_{0:05d}'.format(num_adopted),images=atoms,
                      format='vasp',direct=True,vasp5=True,sort=True)
        else:
            #...revert to previous state
            #atoms.set_chemical_symbols(prev_symbols)
            symbols[:] = prev_symbols[:]
            atoms = atoms_from_symbols(atoms0,symbols)
            # print('not adopted:',erg-ergp,prob,rand)
            erg = ergp
        # print(atoms.get_chemical_formula())
        txt = ' {0:8d} {1:15.7e} {2:s}'.format(istp,erg,atoms.get_chemical_formula())
        #txt = '{0:8d} {1:15.7e}'.format(istp,erg)
        if adopted:
            txt += txt_adopted
        adopted = False
        print('istp,erg='+txt)
        outerg.write(txt+'\n')
        istp += 1
        sys.stdout.flush()

    outsmbl.write('{0:8d}  '.format(num_adopted)+dumps_symbols(symbols))
    write('POSCAR_{0:05d}'.format(num_adopted),images=atoms,
          format='vasp',direct=True,vasp5=True,sort=True)
    outerg.close()
    print('num_adopted / nsteps, rate = {0:5d} /{1:5d}, '.format(num_adopted,nsteps)
          +'{0:7.1f}%'.format(float(num_adopted)/nsteps*100))
    return symbols


def run_kmc(atoms0,nsolute={},
            sites=[],
            symbols=[],
            nsteps=100,
            temperature=300.0,
            pmd_relax=False,
            pmd_cmd='~/src/nap/pmd/pmd',
            paranode=[1,1,1],
            relax_steps=-1,
            movable=[]):
    """
    Run kinetic MC using pmd.

    Atoms0 object must be kept to original fcc Al system.
    To modify this use atoms_from_symbols() fucntion with
    modified symbols array.

    Assuming there is only one vacancy that moves.
    """
    chempots = calc_chemical_potentials(atoms0,
                                        solutes=nsolute.keys(),
                                        pmd_relax=pmd_relax,
                                        pmd_cmd=pmd_cmd,
                                        paranode=paranode)
    print('chempots:')
    for k,v in chempots.items():
        print('  {0:s}:  {1:12.5f}'.format(k,v))

    #...define frequency prefactors in sec^{-1}
    #...phonon frequency of Al at BZ boundary is 4~6 THz
    freq = {'Al':5.e+12,
            'Mg':5.e+12,
            'Si':5.e+12}
    print('freq:')
    for k,v in freq.items():
        print('  {0:s}:  {1:12.5f}'.format(k,v))
    
    #...define average migration barrier in eV
    demig = {'Al':0.569,
             'Mg':0.450,
             'Si':0.479}
    print('demig:')
    for k,v in demig.items():
        print('  {0:s}:  {1:12.5f}'.format(k,v))

    if relax_steps < 0:
        nrstps = len(atoms0)
    else:
        nrstps = relax_steps

    nnmax = max(nlspr)
    #...history arrays for enhancing MC speed
    #...by avoid the energy calcs of the identical systems
    hist_symbols = []
    hist_erg = []
    hist_connects = []

    num_keep = 0
    for i in range(len(symbols)):
        if symbols[i] in nsolute.keys():
            num_keep += 1
    print('total num of Mg,Si,Vac to be kept = ',num_keep)

    sys.stdout.flush()

    calc = PMD(label='pmd', command=pmd_cmd,
               force_type='NN', cutoff_radius=5.8,
               specorder=['Al','Mg','Si'],
               mass=[26.982,24.305,28.085],
               num_nodes_x=paranode[0],
               num_nodes_y=paranode[1],
               num_nodes_z=paranode[2])
    erg0 = symbols_to_erg(atoms0,symbols,calc,chempots,
                          pmd_relax,nrstps)
    ergp = erg0
    hist_symbols.append(symbols)
    hist_erg.append(erg0)

    outerg = open(_ferg,'w')
    outsmbl = open(_fsymbols,'w')
    tclck = 0.0
    txt = ' {0:8d} {1:11.3e} {2:15.7e}'.format(0,tclck,erg0)
    print('istp,tclck,erg='+txt)
    outerg.write(txt+'\n')
    outsmbl.write('{0:8d}  '.format(0)+dumps_symbols(symbols))
    atoms = atoms_from_symbols(atoms0,symbols)
    write('POSCAR_{0:05d}'.format(0),images=atoms,
          format='vasp',direct=True,vasp5=True,sort=True)
    sys.stdout.flush()
    istp = 1
    adopted = False
    prev_symbols = copy.copy(symbols)
    while istp <= nsteps:
        check_num_solutes(symbols,prev_symbols)
        prev_symbols = copy.copy(symbols)
        ichosen = choose(symbols,movable)
        #ichosen = choose(symbols,nsolute.keys())
        #ichosen = choose(symbols,('Vac',))
        ispcs = symbols[ichosen]

        nn = nlspr[ichosen]
        #...check if there are atoms with different symbols in lspr
        ncandidate = 0
        for inn in range(nn):
            jc = lspr[ichosen,inn]
            js = symbols[jc]
            if js != ispcs:
                ncandidate += 1
        if ncandidate == 0: continue
        #...evaluate every migration event
        tmpsymbs = []
        tmpergs = []
        probs = []
        ncalc = 0
        tmpspcs = []
        for inn in range(nn):
            jchosen = lspr[ichosen,inn]
            jspcs = symbols[jchosen]
            if jspcs == ispcs: continue
            tmpsymb = copy.copy(symbols)
            #...exchange positions or migrate vacancy
            tmpsymb[ichosen] = jspcs
            tmpsymb[jchosen] = ispcs
            #...if you want to avoid some of energy calcs,
            #...you should some code here
            ihistory = symbols_exist(tmpsymb,hist_symbols)
            if ihistory == -1: # no same symbols
                erg = symbols_to_erg(atoms0,tmpsymb,calc,
                                     chempots,pmd_relax,nrstps)
                hist_symbols.append(tmpsymb)
                hist_erg.append(erg)
                ncalc += 1
            else: # there is one same symbols
                erg = hist_erg[ihistory]
            tmpergs.append(erg)
            tmpsymbs.append(tmpsymb)
            tmpspcs.append(jspcs)
            #...since ispcs should be vacancy,
            #...jspcs is the migrating atom
            de = demig[jspcs] +(erg-ergp)/2
            prob = freq[jspcs] *np.exp(-de/(temperature*kB))
            probs.append(prob)
        ptot = sum(probs)
        rand = np.random.random()*ptot
        ievent = len(probs)-1
        ptmp = 0.0
        # print('ptot,rand=',ptot,rand)
        for i,pi in enumerate(probs):
            ptmp += pi
            # print('i,pi,ptmp,rand=',i,pi,ptmp,rand)
            if rand < ptmp:
                ievent = i
                break
        #...event was chosen
        ergp = tmpergs[ievent]
        symbols = copy.copy(tmpsymbs[ievent])
        #...progress clock
        rand = np.random.random()
        dt = -np.log(rand)/ptot
        tclck += dt
        #...output new structure
        outsmbl.write('{0:8d}  '.format(istp)+dumps_symbols(symbols))
        if tmpspcs[ievent] != 'Al':
            atoms = atoms_from_symbols(atoms0,symbols)
            write('POSCAR_{0:05d}'.format(istp),images=atoms,
                  format='vasp',direct=True,vasp5=True,sort=True)
        #...stdout
        txt = ' {0:8d} {1:11.3e} {2:15.7e}'.format(istp,tclck,ergp)+\
              ', ncalc={0:3d}'.format(ncalc) +\
              ', ievent={0:2d}'.format(ievent) +\
              ', {0:s}'.format(tmpspcs[ievent])
        print('istp,tclck,erg='+txt)
        outerg.write(txt+'\n')
        istp += 1
        sys.stdout.flush()
        outsmbl.flush()
        outerg.flush()

    outsmbl.write('{0:8d}  '.format(istp)+dumps_symbols(symbols))
    write('POSCAR_{0:05d}'.format(istp),images=atoms,
          format='vasp',direct=True,vasp5=True,sort=True)
    outerg.close()
    print('num of history = ',len(hist_symbols))
    return symbols


if __name__ == "__main__":

    args = docopt(__doc__,version=__version__)
    print(args)
    nsteps = int(args['-n'])
    temperature = float(args['--temperature'])
    ase_relax = args['--ase-relax']
    pmd_relax = not ase_relax
    kinetic = args['--kinetic']
    lattconst = float(args['--latt-const'])
    size = [ int(x) for x in args['--size'].split('x')]
    nMg = int(args['--nMg'])
    nSi = int(args['--nSi'])
    nVac= int(args['--nVac'])
    init_clustered = args['--init-clustered']
    restart = args['--restart']
    paranode = [ int(x) for x in args['--parallel'].split('x') ]
    npara = paranode[0]*paranode[1]*paranode[2]
    nrstps = int(args['--relax-steps'])
    movable = args['--movable'].split(',')
    if npara < 1:
        raise ValueError('npara < 1')

    atoms = make_fcc_Al(lattconst=lattconst,
                        size=size)
    sites = atoms.get_scaled_positions()
    symbols = atoms.get_chemical_symbols()

    #...make permanent neighbor list
    nlspr,lspr = make_pair_list(atoms,rcut=3.0)

    if pmd_relax:
        print('optimization by pmd')
    else:
        print('optimization by ase')

    nsolute = {'Mg':nMg,
               'Si':nSi,
               'Vac':nVac}
    print('nslute=',nsolute)
    if restart:
        with open(_fsymbols,'r') as f:
            lines = f.readlines()
            symbols = loads_symbols(lines[-1].split()[1])
        if len(symbols) != len(sites):
            raise ValueError('len(symbols) != len(sites)')
        for s in nsolute.keys():
            if symbols.count(s) != nsolute[s]:
                raise RuntimeError("symbols.count(s) != nsolute[s], s=",s)
    elif init_clustered:
        symbols = make_symbols_clustered(sites,symbols,
                                         nlspr,lspr,
                                         nsolute0=nsolute)
    else:
        symbols = make_symbols_random(sites,symbols,
                                      nlspr,lspr,
                                      nsolute0=nsolute)
    sys.stdout.flush()

    if npara > 1:
        pmd_cmd = 'mpirun -np {0:d} '.format(npara) \
                  +HOME+'/src/nap/pmd/pmd'
    else:
        pmd_cmd = HOME+'/src/nap/pmd/pmd'
    print('pmd_cmd = ',pmd_cmd)
    if not kinetic:
        print('>>> Metropolis MC')
        symbols = run_lmc(atoms,nsolute=nsolute,
                          sites=sites,
                          symbols=symbols,
                          nsteps=nsteps,
                          temperature=temperature,
                          pmd_relax=pmd_relax,
                          pmd_cmd=pmd_cmd,
                          paranode=paranode,
                          relax_steps=nrstps,
                          movable=movable)
    elif kinetic:
        print('>>> kinetic MC')
        symbols = run_kmc(atoms,nsolute=nsolute,
                          sites=sites,
                          symbols=symbols,
                          nsteps=nsteps,
                          temperature=temperature,
                          pmd_relax=pmd_relax,
                          pmd_cmd=pmd_cmd,
                          paranode=paranode,
                          relax_steps=nrstps,
                          movable=movable)
    
