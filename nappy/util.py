#!/usr/bin/env python
"""
Utility functions for nappy.

Usage:
  util.py [options]

Options:
  -h, --help  Show this message and exit.
"""
import os,sys
from docopt import docopt

__author__ = "RYO KOBAYASHI"
__version__ = ""

def decode_tag(tag):
    """
    Decode tag used in pmd file.
    """
    sid = int(tag)
    ifmv = int((tag -sid)*10)
    num = int(((tag-sid)*10 -ifmv)*1e+14)
    return sid,ifmv,num

def get_tag(sid,ifmv,num):
    tag = sid +ifmv*0.1 +num*1e-14
    return tag

def cartesian_to_scaled(hi,xc,yc,zc):
    """
    Convert an atomic position in Cartesian coordinate
    to scaled position using inversed h-matrix.
    Inversed h-matrix has to be given.
    """
    x = 0.0
    y = 0.0
    z = 0.0
    x += hi[0,0]*xc +hi[0,1]*yc +hi[0,2]*zc
    y += hi[1,0]*xc +hi[1,1]*yc +hi[1,2]*zc
    z += hi[2,0]*xc +hi[2,1]*yc +hi[2,2]*zc
    return x,y,z


def scaled_to_cartesian(h,xs,ys,zs):
    """
    Convert a scaled positions to Cartesian coordinate.
    H-matrix has to be given.
    """
    xc = 0.0
    yc = 0.0
    zc = 0.0
    xc += h[0,0]*xs +h[0,1]*ys +h[0,2]*zs
    yc += h[1,0]*xs +h[1,1]*ys +h[1,2]*zs
    zc += h[2,0]*xs +h[2,1]*ys +h[2,2]*zs
    return xc,yc,zc


def get_axis_and_angle(v,u):
    """
    Get rotation axis and angle between given two vectors v and u.
    """
    import numpy as np
    lv = np.linalg.norm(v)
    lu = np.linalg.norm(u)
    cs = np.dot(v,u)/lv/lu
    vxu = np.cross(v,u)
    axis = vxu/lv/lu
    sn = np.linalg.norm(vxu)/lv/lu
    ang = np.arccos(cs)
    return axis, ang


def rotate(vector,axis,ang):
    """
    Rotate the given *vector* around the *axis* by *ang*.
    *axis* should be normalized vector.
    """
    import numpy as np
    rmat = np.zeros((3,3),dtype=float)
    nx,ny,nz = axis[:]
    rmat[0,:] = [ 0., -nz, ny]
    rmat[1,:] = [ nz, 0., -nx]
    rmat[2,:] = [-ny, nx, 0.]
    mmat = np.zeros((3,3),dtype=float)
    imat = np.identity(3)
    rmat2 = np.dot(rmat,rmat)
    mmat[:,:] = imat[:,:] +np.sin(ang)*rmat[:,:] \
                +(1.0 -np.cos(ang))*rmat2[:,:]
    return np.dot(mmat,vector)

def pbc(x):
    if x < 0.:
        return x -int(x) +1.0
    elif x >= 1.0:
        return x -int(x)
    else:
        return x
    
def parse_option(line):
    """
    Parse an optional variable from a comment line.
    The comment line should start with '#'
    and the option name should follow it with follwoing ':' like,
    # optname:  optvalue
    """
    data = line.split()
    if data[0] != '#':
        return None
    if data[1][-1] != ':':
        return None
    optname = data[1].strip(':')
    optval = data[2]
    #...Return optname and optval as a dictionary.
    return {optname:optval}

def gen_header(argv):
    """
    Generate header that contains command line itself and datetime string.
    ARGV is usually 'sys.argv'.
    """
    from datetime import datetime
    tnow = datetime.now()
    dstr = tnow.strftime('%Y-%m-%d %H:%M:%S')
    header = '# Generated at {0:s} \n'.format(dstr)
    if len(argv) > 10:
        argv = argv[:10]
        argv.append('...')
    cmd = ' '.join(argv)
    header += '#           by {0:s}\n'.format(cmd)
    return header

if __name__ == "__main__":

    args = docopt(__doc__)
