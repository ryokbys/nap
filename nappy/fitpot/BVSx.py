#!/usr/bin/env python
"""
BVSx potential functions.
"""
import os
import numpy as np

__author__ = "RYO KOBAYASHI"
__version__ = "260226"


def read_params_Morse(infname='in.params.Morse'):
    """
    元素ペアのパラメータファイルを読み込み、辞書型で返す関数。
    
    Returns:
        dict: {(cspi, cspj): {'D': float, 'alpha': float, 'rmin': float}, ...}
    """
    params_dict = {}
    
    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    with open(infname, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
                
            # 空行や「#」で始まるコメント行をスキップ
            if not line or line.startswith('#'):
                continue
                
            # スペース区切りで分割
            parts = line.split()
                
            # カラム数が足りない場合はスキップ
            if len(parts) < 5:
                continue
                
            cspi = parts[0]
            cspj = parts[1]
            try:
                # 各数値をfloatに変換
                d_val = float(parts[2])
                alpha_val = float(parts[3])
                rmin_val = float(parts[4])
                        
                # 元素ペアをタプルとしてキーにする
                params_dict[(cspi, cspj)] = {
                    'D': d_val,
                    'alpha': alpha_val,
                    'rmin': rmin_val
                }
            except ValueError:
                # 数値変換できない行（ヘッダー文字列など）は無視
                continue

    return params_dict


def read_params_Coulomb(infname='in.params.Coulomb'):
    """
    Read params in in.params.Coulomb and return them in dictionary format.
    """
    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    keywords = ('terms',
                'interactions',
                'fbvs',
                'charges')
    params = {}
    terms = None
    charges = None
    with open(infname,'r') as f:
        lines = f.readlines()
        mode = None
        for line in lines:
            if line[0] in ('#', '!'):
                continue
            data = line.split()
            if len(data) == 0:
                continue
            kwdline = False
            if data[0] in keywords:
                kwdline = True
                mode = data[0]
            if kwdline:
                if mode == 'terms':
                    terms = data[1]
                    params['terms'] = data[1]
                elif mode == 'fbvs':
                    params['fbvs'] = float(data[1])
                elif mode == 'charges':
                    charges = data[1]
                    params[mode] = {}
                else:
                    params[mode] = []
                continue
            if mode == 'charges' and charges == 'fixed_bvs':
                name = data[0]
                vid = float(data[1])
                rad = float(data[2])
                npq = int(data[3])
                params[mode][name] = (vid,rad,npq)
                # params[mode].append((name,vid,rad,npq))
            elif mode == 'interactions':
                spi = data[0]
                spj = data[1]
                params[mode].append((spi,spj))
            elif mode in ('terms','fbvs'):
                pass
            else:
                raise ValueError('Something wrong with mode: ',mode)
    return params


def read_params_angular(infname='in.params.angular'):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    with open(infname,'r') as f:
        lines = f.readlines()

    angular_prms = {}
    for l in lines:
        if l[0] in ('!','#'):
            continue
        data = l.split()
        if len(data) < 8:
            continue
        if not data[0].isdigit() and not data[1].isdigit() \
           and not data[2].isdigit():
            atype = data[0]
            cspi = data[1]
            cspj = data[2]
            cspk = data[3]
            pdict = {
                'rc3': float(data[4]),
                'alp': float(data[5]),
                'bet': float(data[6]),
                'gmm': float(data[7])
            }
            angular_prms[(cspi,cspj,cspk)] = pdict
        else:
            continue
    return angular_prms


def calc_Morse(D0ij, alpij, rminij, rmin=0.1, rmax=6.0):
    """
    Morse potential curve.
    """
    rs = np.linspace(rmin,rmax,200)
    pot = np.zeros_like(rs)
    potc = D0ij*((np.exp(alpij*(rminij-rmax))-1.0)**2 -1.0)
    pot = D0ij*((np.exp(alpij*(rminij-rs))-1.0)**2 -1.0) -potc
    
    return rs, pot


def calc_Coulomb(terms, qi, qj, fbvs, radi, radj,
                rmin=0.1, rmax=6.0, dielec=1.0):
    """
    Coulomb potential curve of fixed_bvs.
    """
    from scipy.special import erfc
    k = 14.3998554737
    rhoij = fbvs * (radi + radj)
    sqpi = np.sqrt(np.pi)
    terfcc = erfc(rmax/rhoij)
    vrc = k * qi * qj/ rmax * terfcc
    dvdrc = -k / rmax * (terfcc/rmax
                         +2.0/rhoij*sqpi*np.exp(-(rmax/rhoij)**2))
    
    rs = np.linspace(rmin,rmax,200)
    pot = np.zeros_like(rs)
    terfc = erfc(rs / rhoij)
    pot = 0.5 * ( k * qi * qj / rs * terfc - vrc -dvdrc*(rs-rmax)) /dielec
    
    return rs, pot



def calc_angular(rc3, alp, bet, gmm,):
    """
    Angular potential curve.
    """
    ts = np.linspace(0.0, np.pi, 200)
    pot = np.zeros_like(ts)

    csn = np.cos(ts)
    tcsn = csn - gmm
    pot = alp * tcsn**2
    return ts, pot



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


