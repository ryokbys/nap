#!/usr/bin/env python
"""
ZBL potential functions.
"""
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "260208"

azbl = (0.1818, 0.5099, 0.2802, 0.02817)
bzbl = (3.2, 0.9423, 0.4029, 0.2016)


def dzbl(r, zi, zj):
    k = 14.4  # 1/(4*pi*eps0) in eV*Ang/e^2 unit
    rs = 0.4683766 / ( zi**0.23 + zj**0.23 )
    return (k*zi*zj) * (-1.0/r**2 *phi(r/rs) +1.0/r/rs *dphi(r/rs))


def dddzbl(r, zi, zj):
    k = 14.4  # 1/(4*pi*eps0) in eV*Ang/e^2 unit
    rs = 0.4683766 / ( zi**0.23 + zj**0.23 )
    return (k*zi*zj) *( -6.0/r**4 *phi(r/rs) +2.0/r**3/rs*dphi(r/rs) \
                       +2.0/r**3/rs*dphi(r/rs) -1.0/r**2/rs**2*ddphi(r/rs) \
                       +2.0/r**3/rs*dphi(r/rs) -1.0/r**2/rs**2*ddphi(r/rs) \
                       -1.0/r**2/rs**2*ddphi(r/rs) +1.0/r/rs**3*dddphi(r/rs) )


def phi(x):
    #return 0.1818*np.exp(-3.2*x) +0.5099*np.exp(-0.9423*x) \
    #    +0.2802*np.exp(-0.4029*x) +0.02817*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += azbl[i]*np.exp(-bzbl[i]*x)
    return s


def dphi(x):
    #return -0.58176*np.exp(-3.2*x) -0.48047877*np.exp(-0.9423*x) \
    #    -0.11289258*np.exp(-0.4029*x) -0.005679072*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += -azbl[i] * bzbl[i] * np.exp(-bzbl[i]*x)
    return s


def ddphi(x):
    #return 1.861632*np.exp(-3.2*x) +0.452755144971*np.exp(-0.9423*x) \
    #    +0.045484420482*np.exp(-0.4029*x) +0.0011449009152*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += azbl[i] * bzbl[i]**2 * np.exp(-bzbl[i]*x)
    return s


def dddphi(x):
    #return -5.9572224*np.exp(-3.2*x) -0.426631173106*np.exp(-0.9423*x) \
    #    -0.0183256730122*np.exp(-0.4029*x) -0.000230812024504*np.exp(-0.2016*x)
    s = 0.0
    for i in range(4):
        s += -azbl[i] * bzbl[i]**3 * np.exp(-bzbl[i]*x)
    return s


def correct2b(uf3l_prms, config):
    from nappy.elements import get_number_from_symbol
    pairs = config['pairs']
    pair_rins = config['pair_rins']

    uf32b = uf3l_prms['2B']

    for pair in uf32b.keys():
        si, sj = pair
        ipair = get_comb_index(pair, pairs)
        zi = get_number_from_symbol(si)
        zj = get_number_from_symbol(sj)
        coefs = uf32b[pair]['coefs']
        knots = uf32b[pair]['knots']
        prin = pair_rins[ipair]
        if prin > max(knots):
            print(f' Since pair_rin > max(knots) for {si}-{sj}, '
                  +f' correct pari_rin to {max(knots):.2f}.')
            prin = 0.999 * max(knots)
        new_coefs = correct_coefs_wZBL(prin,
                                       knots, coefs,
                                       zi=zi, zj=zj)
        uf3l_prms['2B'][pair]['coefs'] = new_coefs

    return uf3l_prms


