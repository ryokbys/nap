#!/usr/bin/env python
"""
Output a markdown table of the summary of fp results.

Usage:
  fp2pmarkdown.py [options] <out.fp>

Options:
  -h, --help  Show this message and exit.
"""
import os

from docopt import docopt

__author__ = "Ryo KOBAYASHI"
__version__ = "221126"

def read_outfp(fname='out.fp'):
    with open(fname,'r') as f:
        lines = f.readlines()
    bestiid = -1
    for line in reversed(lines):
        if 'best_iid,best_loss' in line:
            bestiid = int(line.split()[3])
            break
    if bestiid < 0:
        raise ValueError('best_iid does not exist !!!')
    targets = []
    weights = []
    losses = []
    mode = ''
    for line in lines:
        if 'weights:' in line:
            mode = 'weights'
            continue
        if '# iid,losses' in line:
            data = line.split()
            targets = [ t for t in data[3:] ]
            continue
        elif 'iid,losses' in line and ' {0:d} '.format(bestiid) in line:
            if len(targets) < 1:
                raise ValueError('len(targets) < 1 !!!')
            data = line.split()
            losses = [ float(l) for l in data[2:] ]
            break
        if mode == 'weights':
            data = line.split()
            if len(data) != 2:
                mode = ''
                continue
            weights.append(float(data[1]))

    return bestiid,targets,weights,losses

def main():
    import os,sys
    args = docopt(__doc__)
    outfp = args['<out.fp>']

    bestiid,targets,weights,losses = read_outfp(outfp)

    #...Markdown table
    maxlen = 0
    wxls = []
    for i,t in enumerate(targets):
        maxlen = max(maxlen,len(t))
        l = losses[i]
        if t != 'total':
            w = weights[i]
            wxls.append(w*l)
        else:
            wxls.append(l)
    txt = '\n'
    txt += f'Best iid = {bestiid:d}\n\n'
    txt += f'|  Target  |  Weight |  Loss  | Weight x Loss |\n'
    txt +=  '|----------|---------|--------|---------------|\n'
    for i in range(len(targets)):
        t = targets[i]
        if t != 'total':
            w = weights[i]
            l = losses[i]
            wxl = wxls[i]
            txt += f'|  {t:s}  |  {w:5.3f}  |  {l:.4f}  |  {wxl:.4f}  |\n'
        else:
            l = losses[i]
            txt += f'|  {t:s}  |  -----  |  -----  |  {l:.4f}  |\n'
    print(txt)
    
    return None

if __name__ == "__main__":

    main()
