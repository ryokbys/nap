def fmake_pairlist(nsys,rcut=3.0,iprint=1,l1st=True,nnmax=10):
    """
    Make pairlist of given napsys by calling Fortran module pairlist.
    """
    import numpy as np
    try:
        # import nappy.pmd.mods as pmods
        import nappy.pmd.pmd_wrapper as pw
    except:
        raise
    natm = nsys.num_atoms()
    pos = nsys.get_scaled_positions()
    hmat = nsys.get_hmat()
    hmati = nsys.get_hmat_inv()
    tags = nsys.get_tags()

    rcut = rcut
    iprint = iprint
    l1st = l1st

    #...Estimate nnmax
    vol = nsys.get_volume()
    rho = max(float(natm)/vol, 0.2)
    nnmax = max(int(1.2 *rho *4*np.pi *rcut**3 /3),nnmax)  # with margin 20%
    # print('nnmax=',nnmax)
    
    #...The hmat array goes to mk_lspr_sngl with the same (i,j)-element,
    #...and it is already consist of (a,b,c) lattice vectors,
    #...so we dont need to transpose it.
    # lspr = pmods.pairlist.mk_lspr_sngl(natm,nnmax,tags,pos.T,
    #                                    rcut,hmat,hmati,
    #                                    iprint,l1st)
    lspr = pw.wrap_lspr_sngl(pos.T,tags,hmat,hmati,rcut,
                             iprint,l1st,nnmax)
    #...lspr comes back as (nnmax+1,natm)-shape arrays
    return lspr.T
