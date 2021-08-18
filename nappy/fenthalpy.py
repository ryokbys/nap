#!/usr/bin/env python
"""
Compute formation enthalpy from given structures using pmd.
The 1st file in the arguments is the product and the following files are the reactants.
If --ergs option is specified, pmd will not be performed and use the energies instead.

Usage:
  fenthalpy.py [options] FILES [FILES...]

Options:
  -h, --help   Show this message and exit.
  --dry        Dry run. Calculate only coefficients of reactants.
  --ergs ERGS  Energies of product and reactants in the order corresponding to given files if available. Comma separated.
               If provided, not to perform MD relaxation. [default: None]
  --nstp NSTP  Number of steps for relaxation MD. [default: 1000]
  --dt DT      Time interval for relaxation MD. [default: -2.0]
  --out4fp     Write out to a file the fp.py data format. 
  --outfname OUTFILE
               Output file name for out4fp. [default: data.pmd.fenth]
"""
import os,sys
from datetime import datetime
from docopt import docopt
import numpy as np

import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "210809"

def get_unit_comp(nsys):
    """
    Get an unit formula from arbitrary composition.

    Input:
        nsys: NAPSystem

    Output:
        gcd: integer
            Greatest common divider.
        unum_sp: list of integer
            List of number of species that is the unit composition.
    """
    num_species = nsys.natm_per_species()
    gcd = np.gcd.reduce(num_species)
    unum_sp = [ n/gcd for n in num_species ]
    return gcd,unum_sp

def get_reactant_coeffs(reactants,product):
    """
    Get coefficients for reactant structures by solving linear algebra Ax=b.
    This is not available for the case the number of reactants is greater than the number of species.

    Input:
        reactants:  list of NAPSystems
        product:  NAPSystem

    Output:
        coeffs:  numpy array of coeffcients
    """
    specorder = product.specorder
    b_vec = np.array([ float(bi) for bi in product.natm_per_species()])
    #...Construct A_mat
    nsp = len(specorder)
    nreact = len(reactants)
    A_mat = np.zeros((nsp,nreact))
    for ir,reactant in enumerate(reactants):
        specorder_ir = reactant.specorder
        natm_per_sp = reactant.natm_per_species()
        for isp,sp in enumerate(specorder):
            if not sp in specorder_ir:
                A_mat[isp,ir] = 0.0
            else:
                sp_in_spir = specorder_ir.index(sp)
                A_mat[isp,ir] = float(natm_per_sp[sp_in_spir])
    print(' A_mat = ',A_mat)
    print(' b_vec = ',b_vec)
    #...Since nreact could be nsp, A_mat may not have inverse, 
    #...so solve minimization of |Ax-b|^2 to obtain x^* vector, x^*=(A^T*A)*A^T*b.
    AA = np.dot(A_mat.T,A_mat)
    AAinv = np.linalg.inv(AA)
    x = np.dot(AAinv,np.dot(A_mat.T,b_vec))
    #...For check
    Ax = np.dot(A_mat,x)
    if len(Ax) != len(b_vec):
        raise ValueError('len(Ax) != len(b_vec)')
    wrong = False
    for i in range(len(b_vec)):
        if abs(Ax[i] -b_vec[i]) > 0.1:
            wrong = True
    if wrong:
        print('Result maybe wrong: i,Ax[i],b[i].')
        for i in range(len(b_len)):
            print(' {0:2d}  {1:5.1f} {2:5.1f}'.format(i,Ax[i],b_vec[i]))
    return x

def calc_formation_enthalpy(ergs_react,erg_prod,coeffs):
    """
    Calculate the formation enthalpy using energies and coefficients of reactants,
    and energy of product.
    """
    if len(ergs_react) != len(coeffs):
        raise ValueError('len(ergs_react) != len(coeffs)')
    dH = erg_prod
    for i in range(len(ergs_react)):
        ei = ergs_react[i]
        ai = coeffs[i]
        dH -= ai*ei
    dH = -dH
    return dH

def get_relaxed_energy(nsys,nstp=1000,dt=-2.0):
    """
    Perform pmd of relaxation and return the final potential energy.
    """
    print(' Relaxing the system: ',nsys.get_chemical_formula())
    pmd = nappy.pmd.PMD(nsys)
    pmd.load_inpmd()
    pmd.set_params(stress_control='vc-Berendsen', pressure_target=0.0,
                   stress_target=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]],
                   stress_relax_time=50.0, print_level=0)
    pmd.run(nstp=nstp,dt=dt,ifdmp=1,dmp=0.99)
    return pmd.result['epot']

def write_fenth_out4fp(fname,dH):
    """
    Write out fenthalpy data in fp.py general data format.

    Parameters:
    -----------
    fname : string
         Name of the output file.
    dH : float
         Formation enthalpy per unit formula.
    """
    with open(fname,'w') as f:
        cmd = ' '.join(s for s in sys.argv)
        f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        f.write('#  {0:s}\n'.format(cmd))
        #...Num of data, weight for the data
        f.write('  {0:6d}  {1:7.3}\n'.format(1, 1.0))
        f.write('  {0:8.3f}\n'.format(dH))
    return None

def main(args):
    files = args['FILES']
    if len(files) < 2:
        raise ValueError('User should provide at least one reactant and one product.')
    dry = args['--dry']
    nstp = int(args['--nstp'])
    dt = float(args['--dt'])
    out4fp = args['--out4fp']
    ergs = args['--ergs']
    if ergs != 'None':
        ergs = [ float(x) for x in args['--ergs'].split(',') ]
        if len(ergs) != len(files):
            raise ValueError('Number of files and ergs are not inconsistent.')
    
    print(' pwd = ',os.getcwd())
    
        
    product = nappy.io.read(files[0])
    reactants = [ nappy.io.read(f) for f in files[1:] ]
    print(' Product: ',product.get_chemical_formula())
    print(' Reactants: ',end='')
    for r in reactants:
        print(r.get_chemical_formula()+', ',end='')
    print('')

    #...Compute coefficients of reactants
    coeffs = get_reactant_coeffs(reactants,product)
    print(' Coefficients: ',)
    for i,r in enumerate(reactants):
        print('   {0:5.2f} for {1:s}'.format(coeffs[i],r.get_chemical_formula()))

    if dry:
        return None

    
    if type(ergs) is list:  # Energies are provided
        erg_prod = ergs[0]
        ergs_react = [ x for x in ergs[1:] ]
    else:  #...Compute relaxation and get potential energies of given structures.
        erg_prod = get_relaxed_energy(product,nstp=nstp,dt=dt)
        ergs_react = [ get_relaxed_energy(r,nstp=nstp,dt=dt) for r in reactants ]
    print(' E of product = {0:.3f}'.format(erg_prod))
    print(' Es of reactnats:')
    for i in range(len(ergs_react)):
        print('   i, E_i = {0:.3f}'.format(ergs_react[i]))

    #...Get formation enthalpy
    dH = calc_formation_enthalpy(ergs_react,erg_prod,coeffs)
    gcd = np.gcd.reduce(product.natm_per_species())
    print(' Formation enthalpy per unit formula:')
    print('   dH = -1*[ {0:.2f} '.format(erg_prod),end='')
    for i,r in enumerate(reactants):
        print('-{0:.2f}*({1:.2f}) '.format(coeffs[i],ergs_react[i]),end='')
    print(']/{0:d}'.format(gcd))
    print('      = {0:.2f} (eV/u.f.)'.format(dH/gcd))

    if out4fp:
        outfname = args['--outfname']
        write_fenth_out4fp(outfname,dH/gcd)
        print(' Wrote {0:s}'.format(outfname))

    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    main(args)
    
