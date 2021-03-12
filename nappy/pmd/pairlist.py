import nappy.pmd.mods as pmods

def fmake_pairlist(nsys,rcut=3.0,iprint=1,l1st=True,nnmax=100):
    """
    Make pairlist of given napsys by calling Fortran module pairlist.
    """
    natm = nsys.num_atoms()
    pos = nsys.get_scaled_positions()
    hmat = nsys.get_hmat()
    hmati = nsys.get_hmat_inv()
    tags = nsys.get_tags()

    rcut = rcut
    nnmax = nnmax
    iprint = iprint
    l1st = l1st

    lspr,d2lspr = pmods.pairlist.mk_lspr_sngl(natm,nnmax,tags,pos.T,
                                              rcut,hmat,hmati,
                                              iprint,l1st)

    return lspr.T, d2lspr.T
