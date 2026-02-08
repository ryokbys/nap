#!/usr/bin/env python
import os, sys
import numpy as np

__description__="""
Utility functions for UF3, UF3L potential.
"""

__author__ = "RYO KOBAYASHI"
__version__ = "260206"


def read_params_uf3(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    uf3_prms = {'1B':{},
                '2B':{},
                '3B':{}}
    mode = 'none'
    body = 'none'
    with open(infname,'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line[0] == '#':
                if len(line) > 3:
                    mode = 'read'
                    continue
                else:
                    mode = 'none'
                    body = 'none'
                    continue
            data = line.split()
            if data[0] in ('1B','2B','3B'):
                body = data[0]
                if body == '1B':
                    spi = data[1]
                    erg = float(data[2])
                    uf3_prms[body][spi] = erg
                    continue
                elif body == '2B':
                    spi = data[1]
                    spj = data[2]
                    uf3_prms[body][(spi,spj)] = {}
                    uf3_prms[body][(spi,spj)]['nlead'] = int(data[3])
                    uf3_prms[body][(spi,spj)]['ntrail'] = int(data[4])
                    uf3_prms[body][(spi,spj)]['spacing'] = data[5]
                    d = f.readline().split()
                    uf3_prms[body][(spi,spj)]['rc2b'] = float(d[0])
                    uf3_prms[body][(spi,spj)]['nknot'] = int(d[1])
                    d = f.readline().split()
                    uf3_prms[body][(spi,spj)]['knots'] = \
                        np.array([ float(x) for x in d])
                    d = f.readline().split()
                    uf3_prms[body][(spi,spj)]['ncoef'] = int(d[0])
                    d = f.readline().split()
                    uf3_prms[body][(spi,spj)]['coefs'] = \
                        np.array([ float(x) for x in d])
                    continue
                elif body == '3B':
                    spi = data[1]
                    spj = data[2]
                    spk = data[3]
                    d3b = {}
                    d3b['nlead'] = int(data[4])
                    d3b['ntrail'] = int(data[5])
                    d3b['spacing'] = data[6]
                    d = f.readline().split()
                    rcjk,rcij,rcik = (float(d[0]), float(d[1]), float(d[2]))
                    nkjk,nkij,nkik = (int(d[3]), int(d[4]), int(d[5]))
                    d3b['rcij'] = rcij
                    d3b['rcik'] = rcik
                    d3b['rcjk'] = rcjk
                    d3b['nkij'] = nkij
                    d3b['nkik'] = nkik
                    d3b['nkjk'] = nkjk
                    d = f.readline().split()
                    d3b['knotsjk'] = \
                        np.array([ float(x) for x in d])
                    d = f.readline().split()
                    d3b['knotsik'] = \
                        np.array([ float(x) for x in d])
                    d = f.readline().split()
                    d3b['knotsij'] = \
                        np.array([ float(x) for x in d])
                    d = f.readline().split()
                    ncij,ncik,ncjk = (int(d[0]), int(d[1]), int(d[2]))
                    d3b['ncij'] = ncij
                    d3b['ncik'] = ncik
                    d3b['ncjk'] = ncjk
                    d3b['coefs'] = np.zeros((ncij,ncik, ncjk))
                    nline3b = ncij*ncik
                    for icij in range(ncij):
                        for icik in range(ncik):
                            d = f.readline().split()
                            d3b['coefs'][icij,icik,:] = \
                                [ float(x) for x in d]
                    uf3_prms[body][(spi,spj,spk)] = d3b

    return uf3_prms


def read_params_uf3l(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    uf3l_prms = {'1B':{},
                 '2B':{},
                 '3B':{}}
    mode = 'none'
    body = 'none'
    with open(infname,'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line[0] == '#':
                if len(line) > 3:
                    mode = 'read'
                    continue
                else:
                    mode = 'none'
                    body = 'none'
                    continue
            data = line.split()
            if data[0] in ('1B','2B','3B'):
                body = data[0]
                if body == '1B':
                    spi = data[1]
                    erg = float(data[2])
                    uf3l_prms[body][spi] = erg
                    continue
                elif body == '2B':
                    spi = data[1]
                    spj = data[2]
                    uf3l_prms[body][(spi,spj)] = {}
                    uf3l_prms[body][(spi,spj)]['nlead'] = int(data[3])
                    uf3l_prms[body][(spi,spj)]['ntrail'] = int(data[4])
                    uf3l_prms[body][(spi,spj)]['spacing'] = data[5]
                    d = f.readline().split()
                    uf3l_prms[body][(spi,spj)]['rc2b'] = float(d[0])
                    uf3l_prms[body][(spi,spj)]['nknot'] = int(d[1])
                    d = f.readline().split()
                    uf3l_prms[body][(spi,spj)]['knots'] = \
                        np.array([ float(x) for x in d])
                    d = f.readline().split()
                    uf3l_prms[body][(spi,spj)]['ncoef'] = int(d[0])
                    d = f.readline().split()
                    uf3l_prms[body][(spi,spj)]['coefs'] = \
                        np.array([ float(x) for x in d])
                    continue
                elif body == '3B':
                    spi = data[1]
                    spj = data[2]
                    spk = data[3]
                    d3b = {}
                    d3b['nlead'] = int(data[4])
                    d3b['ntrail'] = int(data[5])
                    d3b['spacing'] = data[6]
                    d = f.readline().split()
                    rcij, rcik, nknot, gmj, gmk = \
                        (float(d[0]),  float(d[1]), int(d[2]),
                         float(d[3]), float(d[4]))
                    d3b['rcij'] = rcij
                    d3b['rcik'] = rcik
                    d3b['nknot'] = nknot
                    d3b['gmj'] = gmj
                    d3b['gmk'] = gmk
                    d = f.readline().split()
                    d3b['knots'] = np.array([ float(x) for x in d])
                    d = f.readline().split()
                    ncoef = int(d[0])
                    d3b['ncoef'] = ncoef
                    d3b['coefs'] = np.zeros(ncoef)
                    d = f.readline().split()
                    d3b['coefs'] = [ float(x) for x in d ]
                    uf3l_prms[body][(spi,spj,spk)] = d3b

    return uf3l_prms


def read_params_uf3d(infname):

    if not os.path.exists(infname):
        raise FileNotFoundError(infname)

    uf3d_prms = {'1B':[],
                 '2B':[],
                 '3B':[]}
    mode = 'none'
    body = 'none'
    with open(infname,'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line[0] == '#':
                if len(line) > 3:
                    mode = 'read'
                    continue
                else:
                    mode = 'none'
                    body = 'none'
                    continue
            data = line.split()
            if data[0] in ('1B','2B','3B'):
                body = data[0]
                if body == '1B':
                    spi = data[1]
                    erg = float(data[2])
                    uf3d_prms[body].append({'species': spi,
                                            'erg': erg})
                    continue
                elif body == '2B':
                    spi = data[1]
                    spj = data[2]
                    d2b = {'pair': (spi,spj)}
                    d2b['nlead'] = int(data[3])
                    d2b['ntrail'] = int(data[4])
                    d2b['spacing'] = data[5]
                    d = f.readline().split()
                    d2b['rc2b'] = float(d[0])
                    d2b['nknot'] = int(d[1])
                    d = f.readline().split()
                    d2b['knots'] = np.array([ float(x) for x in d])
                    d = f.readline().split()
                    d2b['ncoef'] = int(d[0])
                    d = f.readline().split()
                    d2b['coefs'] = np.array([ float(x) for x in d])
                    uf3d_prms[body].append(d2b)
                    continue
                elif body == '3B':
                    spi = data[1]
                    spj = data[2]
                    spk = data[3]
                    d3b = {'trio': (spi,spj,spk)}
                    d3b['nlead'] = int(data[4])
                    d3b['ntrail'] = int(data[5])
                    d3b['spacing'] = data[6]
                    d = f.readline().split()
                    rcij, rcik, nknij, nknik, nkncs = \
                        (float(d[0]),  float(d[1]), int(d[2]),
                         int(d[3]), int(d[4]))
                    d3b['rcij'] = rcij
                    d3b['rcik'] = rcik
                    d3b['nknij'] = nknij
                    d3b['nknik'] = nknik
                    d3b['nkncs'] = nkncs
                    #...knots
                    d = f.readline().split()
                    d3b['knij'] = np.array([ float(x) for x in d])
                    d = f.readline().split()
                    d3b['knik'] = np.array([ float(x) for x in d])
                    d = f.readline().split()
                    d3b['kncs'] = np.array([ float(x) for x in d])
                    d = f.readline().split()
                    ncfij, ncfik, ncfcs = ( int(d[0]), int(d[1]), int(d[2]) )
                    d3b['ncfij'] = ncfij
                    d3b['ncfik'] = ncfik
                    d3b['ncfcs'] = ncfcs
                    d = f.readline().split()
                    d3b['cfij'] = np.array([ float(x) for x in d ])
                    d = f.readline().split()
                    d3b['cfik'] = np.array([ float(x) for x in d ])
                    d = f.readline().split()
                    d3b['cfcs'] = np.array([ float(x) for x in d ])
                    uf3d_prms[body].append(d3b)

    return uf3d_prms


def write_params_uf3(uf3prms,
                     outfname='in.params.uf3',
                     author=None,
                     overwrite=False):
    from datetime import datetime
    if os.path.exists(outfname) and not overwrite:
        raise Exception(f'{outfname} already exists.')

    f = open(outfname, 'w')

    data1B = uf3prms.get('1B',None)
    data2B = uf3prms.get('2B',None)
    data3B = uf3prms.get('3B',None)

    today = datetime.now().strftime('%Y-%m-%d')
    if author is None:
        author = __author__

    entry_comment = f'#UF3 POT UNITS: metal DATE: {today} AUTHOR: {author} CITATION:\n'
    if data1B is not None:
        spcs = data1B.keys()
        for spi in spcs:
            erg = data1B[spi]
            f.write(entry_comment)
            f.write(f'1B  {spi}  {erg:0.4f}\n')
            f.write('#\n')
    if data2B is not None:
        pairs = data2B.keys()
        for pair in pairs:
            spi,spj = pair
            dp = data2B[pair]
            nlead = dp['nlead']
            ntrail= dp['ntrail']
            spacing = dp['spacing']
            rc2b = dp['rc2b']
            nknot = dp['nknot']
            knots = dp['knots']
            ncoef = dp['ncoef']
            coefs = dp['coefs']
            f.write(entry_comment)
            f.write(f'2B  {spi}  {spj}  {nlead}  {ntrail}  {spacing}\n')
            f.write(f'{rc2b:0.2f}  {nknot}\n')
            for i in range(nknot):
                f.write(f'{knots[i]:0.4f} ')
            f.write('\n')
            f.write(f'{ncoef}\n')
            for i in range(ncoef):
                f.write(f'{coefs[i]:11.4e} ')
            f.write('\n')
            f.write('#\n')
    if data3B is not None:
        trios = data3B.keys()
        for trio in trios:
            spi,spj,spk = trio
            d3b = data3B[trio]
            nlead = d3b['nlead']
            ntrail= d3b['ntrail']
            spacing = d3b['spacing']
            rcij = d3b['rcij']
            rcik = d3b['rcik']
            rcjk = d3b['rcjk']
            nkij = d3b['nkij']
            nkik = d3b['nkik']
            nkjk = d3b['nkjk']
            knotsjk = d3b['knotsjk']
            knotsik = d3b['knotsik']
            knotsij = d3b['knotsij']
            ncij = d3b['ncij']
            ncik = d3b['ncik']
            ncjk = d3b['ncjk']
            coefs = d3b['coefs']
            f.write(entry_comment)
            f.write(f'3B  {spi}  {spj}  {spk}  {nlead}  {ntrail}  {spacing}\n')
            f.write(f'{rcjk:0.2f}  {rcij:0.2f}  {rcik:0.2f}  {nkjk}  {nkij}  {nkik}\n')
            for i in range(nkjk):
                f.write(f'{knotsjk[i]:0.4f} ')
            f.write('\n')
            for i in range(nkik):
                f.write(f'{knotsik[i]:0.4f} ')
            f.write('\n')
            for i in range(nkij):
                f.write(f'{knotsij[i]:0.4f} ')
            f.write('\n')
            f.write(f'{ncij}  {ncik}  {ncjk}\n')
            for icij in range(ncij):
                for icik in range(ncik):
                    for icjk in range(ncjk):
                        f.write(f'{coefs[icij,icik,icjk]:11.4e} ')
                    f.write('\n')
            f.write('#\n')
    f.close()
    return None


def write_params_uf3l(uf3lprms,
                      outfname='in.params.uf3l',
                      author=None,
                      overwrite=False):
    from datetime import datetime
    if os.path.exists(outfname) and not overwrite:
        raise Exception(f'{outfname} already exists.')

    f = open(outfname, 'w')

    data1B = uf3lprms.get('1B', None)
    data2B = uf3lprms.get('2B', None)
    data3B = uf3lprms.get('3B', None)

    today = datetime.now().strftime('%Y-%m-%d')
    if author is None:
        author = __author__

    entry_comment = f'#UF3 POT UNITS: metal DATE: {today} AUTHOR: {author} CITATION:\n'
    if data1B is not None:
        spcs = data1B.keys()
        for spi in spcs:
            erg = data1B[spi]
            f.write(entry_comment)
            f.write(f'1B  {spi}  {erg:0.4f}\n')
            f.write('#\n')
    if data2B is not None:
        pairs = data2B.keys()
        for pair in pairs:
            spi, spj = pair
            dp = data2B[pair]
            nlead = dp['nlead']
            ntrail = dp['ntrail']
            spacing = dp['spacing']
            rc2b = dp['rc2b']
            nknot = dp['nknot']
            knots = dp['knots']
            ncoef = dp['ncoef']
            coefs = dp['coefs']
            f.write(entry_comment)
            f.write(f'2B  {spi}  {spj}  {nlead}  {ntrail}  {spacing}\n')
            f.write(f'{rc2b:0.4f}  {nknot}\n')
            for i in range(nknot):
                f.write(f'{knots[i]:0.4f} ')
            f.write('\n')
            f.write(f'{ncoef}\n')
            for i in range(ncoef):
                f.write(f'{coefs[i]:11.4e} ')
            f.write('\n')
            f.write('#\n')
    if data3B is not None:
        trios = data3B.keys()
        for trio in trios:
            spi, spj, spk = trio
            d3b = data3B[trio]
            nlead = d3b['nlead']
            ntrail = d3b['ntrail']
            spacing = d3b['spacing']
            rcij = d3b['rcij']
            rcik = d3b['rcik']
            nknot = d3b['nknot']
            knots = d3b['knots']
            ncoef = d3b['ncoef']
            coefs = d3b['coefs']
            gmj = d3b['gmj']
            gmk = d3b['gmk']
            f.write(entry_comment)
            f.write(f'3B  {spi}  {spj}  {spk}  {nlead}  {ntrail}  {spacing}\n')
            f.write(f'{rcij:0.3f}  {rcik:0.3f}  {nknot}  {gmj:0.3f}  {gmk:0.3f}\n')
            for i in range(nknot):
                f.write(f'{knots[i]:0.4f} ')
            f.write('\n')
            f.write(f'{ncoef}\n')
            for ic in range(ncoef):
                f.write(f'{coefs[ic]:11.4e} ')
            f.write('\n')
            f.write('#\n')
    f.close()
    return None


def write_params_uf3d(uf3dprms,
                      outfname='in.params.uf3d',
                      author=None,
                      overwrite=False):
    from datetime import datetime
    if os.path.exists(outfname) and not overwrite:
        raise Exception(f'{outfname} already exists.')

    f = open(outfname, 'w')

    data1B = uf3dprms.get('1B', None)
    data2B = uf3dprms.get('2B', None)
    data3B = uf3dprms.get('3B', None)

    today = datetime.now().strftime('%Y-%m-%d')
    if author is None:
        author = __author__

    entry_comment = f'#UF3 POT DATE: {today} AUTHOR: {author} CITATION:\n'
    if data1B is not None:
        for d1 in data1B:
            spi = d1['species']
            erg = d1['erg']
            f.write(entry_comment)
            f.write(f'1B  {spi}  {erg:0.4f}\n')
            f.write('#\n')
    if data2B is not None:
        for dp in data2B:
            spi, spj = dp['pair']
            nlead = dp['nlead']
            ntrail = dp['ntrail']
            spacing = dp['spacing']
            rc2b = dp['rc2b']
            nknot = dp['nknot']
            knots = dp['knots']
            ncoef = dp['ncoef']
            coefs = dp['coefs']
            f.write(entry_comment)
            f.write(f'2B  {spi}  {spj}  {nlead}  {ntrail}  {spacing}\n')
            f.write(f'{rc2b:0.4f}  {nknot}\n')
            for i in range(nknot):
                f.write(f'{knots[i]:0.4f} ')
            f.write('\n')
            f.write(f'{ncoef}\n')
            for i in range(ncoef):
                f.write(f'{coefs[i]:11.4e} ')
            f.write('\n')
            f.write('#\n')
    if data3B is not None:
        for d3b in data3B:
            spi, spj, spk = d3b['trio']
            nlead = d3b['nlead']
            ntrail = d3b['ntrail']
            spacing = d3b['spacing']
            rcij = d3b['rcij']
            rcik = d3b['rcik']
            #...ij
            nknij = d3b['nknij']
            knij = d3b['knij']
            ncfij = d3b['ncfij']
            cfij = d3b['cfij']
            #...ik
            nknik = d3b['nknik']
            knik = d3b['knik']
            ncfik = d3b['ncfik']
            cfik = d3b['cfik']
            #...cos
            nkncs = d3b['nkncs']
            kncs = d3b['kncs']
            ncfcs = d3b['ncfcs']
            cfcs = d3b['cfcs']
            f.write(entry_comment)
            f.write(f'3B  {spi}  {spj}  {spk}  {nlead}  {ntrail}  {spacing}\n')
            f.write(f'{rcij:0.3f}  {rcik:0.3f}  {nknij}  {nknik}  {nkncs}\n')
            #...knots info
            for i in range(nknij):
                f.write(f'{knij[i]:0.4f} ')
            f.write('\n')
            for i in range(nknik):
                f.write(f'{knik[i]:0.4f} ')
            f.write('\n')
            for i in range(nkncs):
                f.write(f'{kncs[i]:0.4f} ')
            f.write('\n')
            #...coefs info
            f.write(f'  {ncfij}  {ncfik}  {ncfcs}\n')
            for ic in range(ncfij):
                f.write(f'{cfij[ic]:11.4e} ')
            f.write('\n')
            for ic in range(ncfik):
                f.write(f'{cfik[ic]:11.4e} ')
            f.write('\n')
            for ic in range(ncfcs):
                f.write(f'{cfcs[ic]:11.4e} ')
            f.write('\n')
            f.write('#\n')
    f.close()
    return None


def to_list(variable):
    if isinstance(variable, np.ndarray):
        return variable.tolist()
    elif isinstance(variable, list):
        return variable
    else:
        raise TypeError('The variable is not list nor ndarray.')


def get_element_list(uf3prms):
    elist = []
    for k,d in uf3prms.items:
        if k == '1B':
            for s in d.keys():
                if e not in elist:
                    elist.append(s)
        elif k == '2B':
            for pair in d.keys():
                si, sj = pair
                if si not in elist:
                    elist.append(si)
                if sj not in elist:
                    elist.append(sj)
        elif k == '3B':
            for trio in d.keys():
                si, sj, sk = trio
                if si not in elist:
                    elist.append(si)
                if sj not in elist:
                    elist.append(sj)
                if sk not in elist:
                    elist.append(sk)
    return elist


def handle_combine(args):
    prms1 = read_params_uf3(args.file1)
    prms2 = read_params_uf3(args.file2)
    prms = merge_dicts(prms1, prms2)
    write_params_uf3(prms, outfname=args.output,
                     author=args.author)
    print(f' wrote {args.output}')
    return None


def merge_dicts(dict1, dict2):
    """
    3段階の入れ子構造を持つ辞書を結合する。
    - 第1段階目のキーは上書きしない。
    - 第2段階目のキーが重複している場合、以下のデータを上書き。

    Args:
        dict1 (dict): 1つ目の辞書（優先）
        dict2 (dict): 2つ目の辞書

    Returns:
        dict: 結合された辞書
    """
    result = dict1.copy()  # dict1をベースにコピー

    for key1, value1 in dict2.items():
        if key1 in result:  # 第1段階目のキーが重複
            if isinstance(result[key1], dict) and isinstance(value1, dict):
                # 両方が辞書の場合、第2段階をマージ
                result[key1] = merge_dicts(result[key1], value1)
            else:
                # 片方が辞書でない場合は上書きしない（辞書構造が異なる場合の保護）
                continue
        else:
            # 第1段階目のキーが重複しない場合、dict2の値を追加
            result[key1] = value1

    return result


def knot_index(r, knots):
    n = 0
    for i, t in enumerate(knots):
        if r > t:
            n = i
        else:
            return n
    return n


def b_spl(r, ts):
    """
    Non-recursive implementation of B-spline.
    Assuming maximum order (d) = 3.

    Args:
      r --- position to be evaluated.
      ts --- knot positions.

    Return:
      nr --- index in knots list.
      b --- B array.
    """
    teps = 1.0e-8
    nr = knot_index(r, ts)

    b = np.zeros((5,4))  # corresponds to btmp(-3:+1,0:3) in fortran
    b[3,0] = 1.0
    for id in range(1,4):  # 1,2,3
        for l in range(-id, 1):  # -id, id+1, ..., 0
            n = nr + l
            tn0 = ts[n]
            tn1 = ts[n+id]
            dt1 = tn1 - tn0
            tmp1 = 0.0
            if abs(dt1) > teps:
                tmp1 = (r-tn0)/dt1 * b[l+3, id-1]
            tn2 = ts[n+1]
            tn3 = ts[n+id+1]
            dt2 = tn3 - tn2
            tmp2 = 0.0
            if abs(dt2) > teps:
                tmp2 = (tn3-r)/dt2 * b[l+3+1, id-1]
            b[l+3, id] = tmp1 + tmp2
    return nr, b[:, 3]


def get_bspl_curve(rs, knots, coefs):
    es = np.zeros(len(rs))
    for i, r in enumerate(rs):
        nr, b = b_spl(r, knots)
        tmp = 0.0
        for l in range(-3, 1):
            n = nr + l
            if n < 0 or n >= len(coefs):
                continue
            c = coefs[n]
            tmp += c * b[l+3]
        es[i] = tmp
    return es



def main():
    import argparse
    # Create arg parser
    parser = argparse.ArgumentParser(description=__description__)

    subparsers = parser.add_subparsers(title="sub-command", dest="command")
    combine = subparsers.add_parser('combine', help="Combine two in.params.uf3 files.")
    combine.add_argument("file1", help="1st file")
    combine.add_argument("file2", help="2nd file")
    combine.add_argument("-o", "--output",
                         type=str,
                         help="output file name.",
                         default="in.params.uf3.combine")
    combine.add_argument("-a", "--author",
                         type=str,
                         help="author name",
                         default="")
    combine.set_defaults(func=handle_combine)


    # Analyze argument
    args = parser.parse_args()

    if args.command is None:
        args.print_help()
    else:
        args.func(args)

if __name__ == "__main__":

    main()
