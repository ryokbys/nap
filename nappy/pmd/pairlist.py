def fmake_pairlist(nsys,rcut=3.0,iprint=1,l1st=True,nnmax=100):
    """
    Make pairlist of given napsys by calling Fortran module pairlist.
    """
    import numpy as np
    try:
        import nappy.pmd.mods as pmods
    except:
        raise
    natm = nsys.num_atoms()
    pos = nsys.get_scaled_positions()
    hmat = nsys.get_hmat()
    hmati = nsys.get_hmat_inv()
    tags = nsys.get_tags()

    rcut = rcut
    nnmax = nnmax
    iprint = iprint
    l1st = l1st

    #...The hmat array goes to mk_lspr_sngl with the same (i,j)-element,
    #...and it is already consist of (a,b,c) lattice vectors,
    #...so we dont need to transpose it.
    lspr,d2lspr = pmods.pairlist.mk_lspr_sngl(natm,nnmax,tags,pos.T,
                                              rcut,hmat,hmati,
                                              iprint,l1st)
    #...lspr, d2lspr come back as (nnmax+1,natm)-shape arrays
    return lspr.T, d2lspr.T
