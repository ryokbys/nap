#!/usr/bin/env python
"""
Compute formation enthalpy from given structures.
If --erg-xxx option is not specified, pmd will be performed to get energies.

Usage:
  fenthalpy.py [options]

Options:
  -h, --help   Show this message and exit.
  --dry        Dry run. Calculate only coefficients of reactants.
  --product PROD
               Only one atom config file of the product. [default: None]
  --reactants REACT
               Atom config files of reactants. Comma separated. [default: None]
  --erg-prod EPROD
               Energy per atom of the product if available. If provided, not to perform MD relaxation. [default: None]
  --ergs-react ERGS
               Energies per atom of reactants if available in the order corresponding to given files. Comma separated.
               If provided, not to perform MD relaxation. [default: None]
  --nstp NSTP  Number of steps for relaxation MD. [default: 1000]
  --dt DT      Time interval for relaxation MD. [default: -2.0]
  --damp DAMP  Damping coefficient. [default: 0.9]
  --out4fp     Write out to a file in the fp.py data format. 
  -o,--outfname OUTFILE
               Output file name for out4fp. [default: data.pmd.fenth]
  --print-level IPRINT
               Print level in pmd. [default: 0]
"""
import os,sys
from datetime import datetime
from docopt import docopt
import numpy as np

import nappy

__author__ = "RYO KOBAYASHI"
__version__ = "rev210809"

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
        if abs(Ax[i] -b_vec[i]) > 0.01:
            wrong = True
    if wrong:
        print(' WARNING: Exact solution was not obtained.')
        print(' Result maybe wrong: i,Ax[i],b[i].')
        for i in range(len(b_vec)):
            print('   {0:2d}  {1:5.1f} {2:5.1f}'.format(i,Ax[i],b_vec[i]))
    else:
        print(' Ax=b is satisfied, which means the exact number relationship between LHS and RHS is found.')
        
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

def get_pmd_done(nsys,nstp=1000,dt=-2.0,print_level=0,damp=0.9):
    """
    Perform pmd of relaxation and return the pmd object.
    """
    print(' Relaxing the system: ',nsys.get_chemical_formula())
    pmd = nappy.pmd.PMD(nsys)
    pmd.load_inpmd()
    pmd.set_params(stress_control='vc-Berendsen', pressure_target=0.0,
                   stress_target=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]],
                   stress_relax_time=50.0, print_level=print_level)
    pmd.run(nstp=nstp,dt=dt,ifdmp=1,dmp=damp)
    nsys_fin = pmd.get_system()
    nappy.io.write(nsys_fin,fname="pmdfin_{0:s}".format(nsys_fin.get_chemical_formula()))
    return pmd

def write_fenth_out4fp(fname,dH,vol):
    """
    Write out formation enthalpy and volume per atom in fp.py general data format.

    Parameters:
    -----------
    fname : string
         Name of the output file.
    dH : float
         Formation enthalpy per atom.
    vol : float
         Volume per atom.
    """
    with open(fname,'w') as f:
        cmd = ' '.join(s for s in sys.argv)
        f.write('# Output at {0:s} from,\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        f.write('#  {0:s}\n'.format(cmd))
        #...Num of data, weight for the data
        f.write('  {0:6d}  {1:7.3f}\n'.format(2, 1.0))
        f.write('  {0:8.3f}  {1:8.3f}\n'.format(dH,vol))
    return None

def main(args):
    product = args['--product']
    if product == 'None':
        raise ValueError('A product must be given via --product option.')
    reactants = [ x for x in args['--reactants'].split(',') ]
    if reactants[0] == 'None':
        raise ValueError('At least one reactants must be given via --reactants option.')
    dry = args['--dry']
    nstp = int(args['--nstp'])
    dt = float(args['--dt'])
    damp = float(args['--damp'])
    out4fp = args['--out4fp']
    erg_prod = args['--erg-prod']
    ergs_react = args['--ergs-react']
    iprint = int(args['--print-level'])
    if erg_prod != 'None':
        erg_prod = float(erg_prod)
    if ergs_react != 'None':
        ergs_react = [ float(x) for x in ergs_react.split(',') ]
        if len(ergs_react) != len(reactants):
            raise ValueError('Number of files and --ergs-react are not consistent with --reactants.')
    
    print(' Working directory: ',os.getcwd())
        
    product = nappy.io.read(product)
    reactants = [ nappy.io.read(f) for f in reactants ]
    print(' Product: ',product.get_chemical_formula())
    print(' Reactants: ',end='')
    for r in reactants:
        print(r.get_chemical_formula()+', ',end='')
    print('')
    
    #...Compute coefficients of reactants
    coeffs = get_reactant_coeffs(reactants,product)
    print(' Coefficients, x_vec: ',)
    for i,r in enumerate(reactants):
        print('   {0:<12s} = {1:>5.2f}'.format(r.get_chemical_formula(),coeffs[i]))
    sys.stdout.flush()

    if dry:
        return None

    if erg_prod == 'None':  # Compute relaxation and get potential energies of given structures.
        pmd_prod = get_pmd_done(product,nstp=nstp,dt=dt,
                                print_level=iprint,damp=damp)
        erg_prod = pmd_prod.result['epot']
        product = pmd_prod.get_system()
    else: # Energy per atom is given
        erg_prod *= len(product)
    if ergs_react == 'None':
        pmds_react = [ get_pmd_done(r,nstp=nstp,dt=dt,print_level=iprint)
                       for r in reactants ]
        ergs_react = [ p.result['epot'] for p in pmds_react ]
        reactants = [ p.get_system() for p in pmds_react ]
    else:
        ergs_react = [ e*len(r) for e,r in zip(ergs_react,reactants) ]
    print(' E of product, {0:s} = {1:.3f}'.format(product.get_chemical_formula(),erg_prod))
    print(' Es of reactants:')
    for i in range(len(ergs_react)):
        r = reactants[i]
        print('   {0:<12s} = {1:>8.3f}'.format(r.get_chemical_formula(),ergs_react[i]))

    #...Get formation enthalpy
    dH = calc_formation_enthalpy(ergs_react,erg_prod,coeffs)
    gcd = np.gcd.reduce(product.natm_per_species())
    print(' Formation enthalpy per formula unit:')
    print('   dH = -1*[ {0:.2f} '.format(erg_prod),end='')
    for i,r in enumerate(reactants):
        print('-{0:.2f}*({1:.2f}) '.format(coeffs[i],ergs_react[i]),end='')
    print(']/{0:d}'.format(gcd))
    print('      = -1*[ {0:.2f} '.format(erg_prod),end='')
    for i,r in enumerate(reactants):
        print('-({0:.2f}) '.format(coeffs[i]*ergs_react[i]),end='')
    print(']/{0:d}'.format(gcd))
    print('      = {0:.2f} (eV/f.u.)'.format(dH/gcd))
    print(' Formation enthalpy per atom:')
    print('      = {0:.3f} (eV/atom)'.format(dH/len(product)))
    print(' at volume per atom:')
    vol = product.get_volume()/len(product)
    print('      = {0:.3f} (Ang^3/atom)'.format(vol))

    if out4fp:
        outfname = args['--outfname']
        write_fenth_out4fp(outfname,dH/len(product),vol)
        print(' Wrote {0:s}'.format(outfname))

    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    main(args)
    
