#!/usr/bin/env python
import os,sys
import numpy as np

__description__="""
Utility functions for UF3, UF3L potential.
"""

__author__ = "RYO KOBAYASHI"
__version__ = "250331"

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
                    epot = float(data[2])
                    uf3_prms[body][spi] = epot
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
                    epot = float(data[2])
                    uf3l_prms[body][spi] = epot
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
                    rc, nknot, betj, betk = (float(d[0]), int(d[1]), float(d[2]), float(d[3]))
                    d3b['rc'] = rc
                    d3b['nknot'] = nknot
                    d3b['betj'] = betj
                    d3b['betk'] = betk
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
            epot = data1B[spi]
            f.write(entry_comment)
            f.write(f'1B  {spi}  {epot:0.4f}\n')
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
