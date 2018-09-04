#!/usr/bin/env python
"""
Create `in.params.NN2` from `in.params.desc`.

Usage:
  desc2nn.py [options]

Options:
  -h, --help  Show this message and exit.
  -i IFNAME   Input file name. [default: in.params.desc]
  -o OFNAME   Output file name. [default: in.params.nn2]
  -l, --num-layer=NUM_LAYER
              Number of hidden layers. [default: 1]
  -n, --num-nodes=NUM_NODES
              Number of nodes in hidden layers. [default: 10]
  -m MAX      Max value of the weights that are randomly determined. [default: 1.0e-2]
"""
from __future__ import print_function

from docopt import docopt
import random
import nappy.nn.desc as desc

__author__ = "RYO KOBAYASHI"
__version__ = "180807"


def create_nn2_input(nsf,num_layer,num_nodes,wmax,fname='in.params.nn2'):
    Ni = [nsf]
    for i in range(len(num_nodes)):
        Ni.append(num_nodes[i])
    Ni.append(1)
    nwgts = 0
    for i in range(1,len(Ni)):
        nwgts += Ni[i-1]*Ni[i]
    print(' Number of weights = {0:d}'.format(nwgts))
    with open(fname,'w') as f:
        f.write(' {0:4d} {1:8d}'.format(num_layer,nsf))
        for i in range(len(num_nodes)):
            f.write(' {0:5d}'.format(num_nodes[i]))
        f.write('\n')
        for i in range(nwgts):
            w = (random.random() -0.5)*2.0 *min(1.0e-2,wmax)
            f.write('  {0:12.4e}  {1:12.4e}  {2:12.4e}\n'.format(w,-wmax,wmax))
    print(' Wrote {0:s}'.format(fname))
    return None

if __name__ == "__main__":

    args = docopt(__doc__)
    infname = args['-i']
    outfname = args['-o']
    num_layer = int(args['--num-layer'])
    num_nodes = [ int(x) for x in args['--num-nodes'].split() ]
    wmax = float(args['-m'])

    nsp,nsf,descs,r_inner = desc.read_desc(infname)

    create_nn2_input(nsf,num_layer,num_nodes,wmax,outfname)
    
    
