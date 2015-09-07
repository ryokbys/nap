#!/usr/bin/env python
"""
Make clusters of samples according to the distance defined using
radial and angular distribution functions.

Usage:
    cluster_samples.py [options] DIR [DIR...]

Options
    -h, --help  Show this help message and exit.
    -n, --num-clusters=N
                Number of clusters to be obtained. [default: 1]
    -l, --load=YAML
                Load `Y` from a given YAML file. [default: None]
"""

import os,sys
import numpy as np
import scipy.cluster.hierarchy as sch
import yaml
from docopt import docopt
sys.path.append(os.path.dirname(__file__)+'/..')
import pmdsys
import rdf
import adf


_rdfname= 'out.rdf'
_adfname= 'out.adf'
_tiny= 0.0001

def read_df(infname):
    #...Count line number
    with open(infname,'r') as f:
        nlines= len(f.readlines())

    df= np.zeros((nlines,2),dtype=np.float)
    with open(infname,'r') as f:
        for i,line in enumerate(f.readlines()):
            buff= line.split()
            df[i,0]= float(buff[0])
            df[i,1]= float(buff[1])
    return df

def get_distance(rdf1,rdf2,adf1,adf2):
    dist= 0.0
    for i in range(len(rdf1)):
        dist += np.abs(rdf1[i,1]-rdf2[i,1])
    for i in range(len(adf1)):
        dist += np.abs(adf1[i,1]-adf2[i,1])
    return dist

def get_dist_matrix(rdfs,adfs):
    dim= len(rdfs)
    D= np.zeros((dim,dim),dtype=np.float)
    for i in range(dim):
        for j in range(dim):
            D[i,j]= get_distance(rdfs[i],rdfs[j],adfs[i],adfs[j])
            D[j,i]= D[i,j]
    return D


if __name__ == '__main__':

    args= docopt(__doc__)

    dirs= args['DIR']
    nclst= int(args['--num-clusters'])
    yafname= args['--load']

    if yafname == 'None':
    
        #...Read out.adf and out.rdf of all the directories.
        print ' Reading rdf and adf files...'
        rdfs= []
        adfs= []
        for dir in dirs:
            rdfs.append(read_df(dir+'/'+_rdfname))
            adfs.append(read_df(dir+'/'+_adfname))
        
        #...Compute distances bewteen all pairs of samples and make a distance matrix.
        print ' Computing a distance matrix...'
        D= get_dist_matrix(rdfs,adfs)
        # print D
    
        #...Perform cluster analysis.
        print ' Performing cluster analysis...'
        Y= sch.linkage(D, method='centroid')
        #...Save Y
        with open('cluster.yaml','w') as f:
            yaml.dump(Y.tolist(),f)
        # print Y

    else:
        print ' loading {0}...'.format(yafname)
        with open(yafname,'r') as f:
            loaded= yaml.load(Y)
        Y= np.array(loaded)
    
    #...Get the list of cluster info
    threshold= Y[-nclst,2] +_tiny
    print ' threshold =',threshold
    lsclst= sch.fcluster(Y, threshold, criterion='distance')
    
    #...Write results
    maxclst= max(lsclst)
    chcklst= np.zeros((maxclst+1,),dtype=np.int)
    selected= []
    # print 'lsclst=',lsclst
    # print 'maxclst=',maxclst
    # print 'chcklst=',chcklst
    for i in range(len(lsclst)):
        if chcklst[lsclst[i]] == 0:
            selected.append(dirs[i])
            chcklst[lsclst[i]] += 1
    with open('out.selected','w') as f:
        for dir in selected:
            f.write(' {0}\n'.format(dir))
    print 'done.'
