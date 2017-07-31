#!/bin/env python
"""
Analyze NN potential with drawing NN structure graph.

Usage:
  analyze.py [options]
  analyze.py draw [options]

Options:
  -h,--help  Show this message and exit.
  -w         Show weight values. [default: False]
  -t THRESHOLD
             Threshold value multiplied to max edge for omitting criterion of edge. [default: 0.01]
"""

import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from docopt import docopt

_cnstfname= 'in.const.NN'
_combfname= 'in.comb.NN'
_paramfname= 'in.params.NN'

def comb(n,m):
    return math.factorial(n)/math.factorial(m)

def read_NN_config():
    fcnst= open(_cnstfname,'r')
    buff= fcnst.readline().split()
    nl= int(buff[0])
    nhl= np.zeros((nl+1),dtype=int)
    nsp= int(buff[1])
    nhl[0]=int(buff[2])
    nhl[1]=int(buff[3])
    if nl == 2:
        nhl[2]= int(buff[4])
    print 'num of species=',nsp
    print 'num of layers =',nl
    print 'num of neurons-{0:d} ='.format(0),nhl[0]
    print 'num of neurons-{0:d} ='.format(1),nhl[1]
    if nl == 2:
        print 'num of neurons-{0:d} ='.format(2),nhl[2]
    n2=0
    n3=0
    ngauss= 0
    ncos= 0
    npoly= 0
    nangle= 0
    for line in fcnst.readlines():
        buff= line.split()
        itype= int(buff[0])
        if itype <= 100:
            n2 += 1
            itype2= itype % 100
            if itype2 == 1:
                ngauss += 1
            elif itype2 == 2:
                ncos += 1
            elif itype2 == 3:
                npoly += 1
        elif itype <= 200:
            n3 += 1
            nangle += 1
    fcnst.close()
    print 'read in.const.NN'
    print 'num of 2body terms=',n2
    print 'num of 3body terms=',n3
    
    ncmb2= nsp +comb(nsp,2)
    ncmb3= ncmb2*nsp
    print 'num of 2body pairs   =',ncmb2
    print 'num of 3body triplets=',ncmb3
    nhl[0]= n2*ncmb2 +n3*ncmb3
    print 'num of 2body and 3body inputs =',n2*ncmb2, n3*ncmb3
    if nl == 1:
        print 'num of input neurons =',nhl[0]*nhl[1] +nhl[1]
    elif nl == 2:
        print 'num of input neurons =',nhl[0]*nhl[1] +nhl[1]*nhl[2] +nhl[2]

    fcmb= open(_combfname,'r')
    cmb2= np.zeros((ncmb2,2),dtype=int)
    cmb3= np.zeros((ncmb3,3),dtype=int)
    print 'pairs:'
    for i2 in range(ncmb2):
        buff= fcmb.readline().split()
        cmb2[i2,0]= int(buff[0])
        cmb2[i2,1]= int(buff[1])
        print '  ',i2,': {0:1d}-{1:1d}'.format(cmb2[i2,0],cmb2[i2,1])
    print 'triplets:'
    for i3 in range(ncmb3):
        buff= fcmb.readline().split()
        cmb3[i3,0]= int(buff[0])
        cmb3[i3,1]= int(buff[1])
        cmb3[i3,2]= int(buff[2])
        print '  ',i3,':', \
            ' {0:1d}-{1:1d}-{2:1d}'.format(cmb3[i3,0],cmb3[i3,1],cmb3[i3,2])
    fcmb.close()
    return nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nangle
    
def read_NN_params(nl,nhl):
    #.....read in.params.NN
    fparam= open(_paramfname,'r')
    buff= fparam.readline().split()
    if nl == 1:
        wgt11= np.zeros((nhl[0],nhl[1]))
        wgt12= np.zeros(nhl[1])
        for ihl0 in range(nhl[0]):
            for ihl1 in range(nhl[1]):
                buff= fparam.readline().split()
                wgt11[ihl0,ihl1]= float(buff[0])
        for ihl1 in range(nhl[1]):
            buff= fparam.readline().split()
            wgt12[ihl1]= float(buff[0])
    elif nl == 2:
        wgt21= np.zeros((nhl[0],nhl[1]))
        wgt22= np.zeros((nhl[1],nhl[2]))
        wgt23= np.zeros(nhl[2])
        for ihl0 in range(nhl[0]):
            for ihl1 in range(nhl[1]):
                buff= fparam.readline().split()
                wgt21[ihl0,ihl1]= float(buff[0])
        for ihl1 in range(nhl[1]):
            for ihl2 in range(nhl[2]):
                buff= fparam.readline().split()
                wgt22[ihl1,ihl2]= float(buff[0])
        for ihl2 in range(nhl[2]):
            buff= fparam.readline().split()
            wgt23[ihl2]= float(buff[0])
        
    
    fparam.close()
    print 'read in.params.NN'
    #for ihl0 in range(n2+n3):
    #    print ihl0,': ',wgt11[ihl0,0:nhl[1]+1]
    # for ihl1 in range(nhl[1]):
    #     print ihl1,': ',wgt12[ihl1]OB
    if nl == 1:
        return wgt11,wgt12
    elif nl == 2:
        return wgt21,wgt22,wgt23


def analyze(nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nagnle,
            wgt11=None,wgt12=None,wgt21=None,wgt22=None,wgt23=None):
    """
    Analyze the NN structure.
    """
    pass


def draw(nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nagnle,
         wgt11=None,wgt12=None,wgt21=None,wgt22=None,wgt23=None):
    g= nx.Graph()
    pos= {}
    dy=-1.0
    for ihl0 in range(nhl[0]):
        g.add_node('0-{0:03d}'.format(ihl0))
        pos['0-{0:03d}'.format(ihl0)]= [0,ihl0*dy]
    dy= -float(nhl[0]-1)/(nhl[1]-1)
    for ihl1 in range(nhl[1]):
        g.add_node('1-{0:03d}'.format(ihl1))
        pos['1-{0:03d}'.format(ihl1)]= [1,ihl1*dy]
    if nl == 1:
        dy= -float(nhl[0])/2
        g.add_node('2')
        pos['2']= [2,dy]
    elif nl == 2:
        dy= -float(nhl[0]-1)/(nhl[2]-1)
        for ihl2 in range(nhl[2]):
            g.add_node('2-{0:03d}'.format(ihl2))
            pos['2-{0:03d}'.format(ihl2)]= [2,ihl2*dy]
        dy= -float(nhl[0])/2
        g.add_node('3')
        pos['3']= [3,dy]

    n= 0
    nlabel= {}
    for key in pos:
        # print key,pos[key]
        if key[0] != '0':
            nlabel[key]= ''
        else:
            ineuron= int(key[2:5])
            # print 'ineuron=',ineuron
            if ineuron < n2*len(cmb2):
                pair= ineuron / n2
                isf2= ineuron % n2
                nlabel[key]= '{0:1d}-'.format(cmb2[pair,0]) \
                             +'{0:1d}:'.format(cmb2[pair,1]) \
                             +' {0:02d}'.format(isf2)
            else:
                ine= ineuron -n2*len(cmb2)
                triplet= ine / n3
                isf3   = ine % n3
                # print ' n3,triplet,isf3=',n3,triplet,isf3
                nlabel[key]= '{0:1d}-'.format(cmb3[triplet,0]) \
                             +'{0:1d}-'.format(cmb3[triplet,1]) \
                             +'{0:1d}:'.format(cmb3[triplet,2]) \
                             +' {0:02d}'.format(isf3)
        # print n,nlabel[n]
        n += 1
    # exit()
    
    maxedge=0.0
    if nl == 1:
        for ihl0 in range(nhl[0]):
            for ihl1 in range(nhl[1]):
                maxedge= max(maxedge,np.abs(wgt11[ihl0,ihl1]))
        for ihl1 in range(nhl[1]):
            maxedge= max(maxedge,np.abs(wgt12[ihl1]))
    elif nl == 2:
        for ihl0 in range(nhl[0]):
            for ihl1 in range(nhl[1]):
                maxedge= max(maxedge,np.abs(wgt21[ihl0,ihl1]))
        for ihl1 in range(nhl[1]):
            for ihl2 in range(nhl[2]):
                maxedge= max(maxedge,np.abs(wgt22[ihl1,ihl2]))
        for ihl2 in range(nhl[2]):
            maxedge= max(maxedge,np.abs(wgt23[ihl2]))
    print 'max of edge value= ',maxedge
        
    colors= []
    elabels= {}
    ic= 0
    if nl == 1:
        for ihl0 in range(nhl[0]):
            for ihl1 in range(nhl[1]):
                if np.abs(wgt11[ihl0,ihl1]) > threshold*maxedge:
                    val= wgt11[ihl0,ihl1]
                    g.add_edge('0-{0:03d}'.format(ihl0),'1-{0:03d}'.format(ihl1))
                    elabels[('0-{0:03d}'.format(ihl0),'1-{0:03d}'.format(ihl1))]='{0:7.4f}'.format(val)
        for ihl1 in range(nhl[1]):
            if np.abs(wgt12[ihl1]) > threshold*maxedge:
                val= wgt12[ihl1]
                g.add_edge('1-{0:03d}'.format(ihl1), '2')
                elabels[('1-{0:03d}'.format(ihl1),'2')]= '{0:7.4f}'.format(val)
    elif nl == 2:
        for ihl0 in range(nhl[0]):
            for ihl1 in range(nhl[1]):
                if np.abs(wgt21[ihl0,ihl1]) > threshold*maxedge:
                    val= wgt21[ihl0,ihl1]
                    g.add_edge('0-{0:03d}'.format(ihl0),'1-{0:03d}'.format(ihl1))
                    elabels[('0-{0:03d}'.format(ihl0),'1-{0:03d}'.format(ihl1))]='{0:7.4f}'.format(val)
        for ihl1 in range(nhl[1]):
            for ihl2 in range(nhl[2]):
                if np.abs(wgt22[ihl1,ihl2]) > threshold*maxedge:
                    val= wgt22[ihl1,ihl2]
                    g.add_edge('1-{0:03d}'.format(ihl1), '2-{0:03d}'.format(ihl2))
                    elabels[('1-{0:03d}'.format(ihl1),'2-{0:03d}'.format(ihl2))]= '{0:7.4f}'.format(val)
        for ihl2 in range(nhl[2]):
            if np.abs(wgt23[ihl2]) > threshold*maxedge:
                val= wgt23[ihl2]
                g.add_edge('2-{0:03d}'.format(ihl2), '3')
                elabels[('2-{0:03d}'.format(ihl2),'3')]= '{0:7.4f}'.format(val)

    for e in g.edges():
        e1= e[0]
        e2= e[1]
        for l in elabels.keys():
            if e1 in l and e2 in l:
                colors.append(np.sqrt(np.abs(float(elabels[l]))))

    # print 'len(edges)=',len(g.edges())
    # print g.edges()
    # print 'len(colors)=',len(colors)
    # print colors
    #exit()
            
    #nx.draw_networkx_nodes(g,pos,node_size=30,node_color='b',node_shape='o')
    #nx.draw_networkx_edges(g,pos)
    nodes= nx.draw_networkx_nodes(g,pos,node_size=30,node_color='b')
    edges= nx.draw_networkx_edges(g,pos,edge_color=colors,edge_cmap=plt.get_cmap('jet'))
    if flag_weight:
        nx.draw_networkx_edge_labels(g,pos,alpha=1.0,edge_labels=elabels,label_pos=0.5)
    for key in pos:
        pos[key][0] -= 0.2
    nx.draw_networkx_labels(g,pos,nlabel,font_size=8)
    
    plt.colorbar(edges)
    plt.tick_params(axis='x',bottom='off',top='off',labelbottom='off')
    plt.tick_params(axis='y',bottom='off',top='off',labelleft='off')
    plt.show()


if __name__ == '__main__':

    args= docopt(__doc__)

    flag_weight= args['-w']
    threshold= args['-t']

    nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nagnle= read_NN_config()
    if nl == 1:
        wgt11,wgt12= read_NN_params(nl,nhl)
    elif nl == 2:
        wgt21,wgt22,wgt23= read_NN_params(nl,nhl)

    if args['draw']:
        if nl == 1:
            draw(nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nagnle,
                 wgt11=wgt11,wgt12=wgt12)
        elif nl == 2:
            draw(nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nagnle,
                 wgt21=wgt21,wgt22=wgt22,wgt23=wgt23)
    else:
        analyze(nl,nsp,nhl,n2,n3,cmb2,cmb3,ngauss,ncos,npoly,nagnle)
