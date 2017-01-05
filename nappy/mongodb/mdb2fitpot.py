#!/usr/bin/env python
"""
Extract data from MongoDB and export them for the use of fitpot.

Usage:
    mdb2fitpot.py [options] CONFIG_FILE

Options:
    -h, --help   Show this help message and exit.
    -q, --query=<queries>
                 Query operations passed to MondoDB in json format. [default: {}]
    -d, --dir=<dir>
                 Directory in which output files are to be saved. [default: samples]
    -o, --offset=<offnset>
                 Offset number of 5-digit directory name. [default: 0]
"""

__author__    = "Ryo KOBAYASHI"
__email__     = "ryo.kbys@gmail.com"
__copyright__ = "Copyright 2015, Ryo KOBAYASHI"
__license__   = "MIT"
__version__   = "0.1"


import os,sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')

from docopt import docopt
from pymongo import MongoClient
import json
import yaml

from napsys import NAPSystem
from atom import Atom

def read_db_config(fname='config.json'):
    """
    Read the DB configuration from json format file.
    """
    with open(fname,'r') as f:
        conf= json.load(f)
    return conf

def load_DB(conf):
    """
    Load DB variables according to the config read from a file.
    """
    client= MongoClient(conf['host'],conf['port'])
    db= client[conf['database']]
    db.authenticate(conf['user'],conf['password'])
    col= db[conf['collection']]
    return client,db,col

def doc_to_pos(doc,conf):
    """
    Make a pos file, which has pmd format, from a document in MongoDB.
    """
    psys= NAPSystem()
    matrix=doc['calculations'][-1]['output']['crystal']['lattice']['matrix']
    a1= matrix[0]
    a2= matrix[1]
    a3= matrix[2]
    psys.set_lattice(1.0,a1,a2,a3)

    species_ids=conf['species_ids']

    sites= doc['calculations'][-1]['output']['crystal']['sites']
    for site in sites:
        ra= site['abc']
        ai= Atom()
        ai.set_pos(ra[0],ra[1],ra[2])
        ai.set_sid(species_ids[site['species'][0]['element']])
        psys.add_atom(ai)
    return psys

def unicode_to_ascii(dic):
    """
    Recursively replace unicode string of value in a given dictionary
    to ascii.
    """
    for k,v in dic.items():
        if type(v) == unicode:
            dic[k]= v.encode('ascii')
        elif type(v) == dict:
            v= unicode_to_ascii(v)
            dic[k]= v
    return dic

def main(conf,dirname='sample',query='{}',offset=0):
    """
    Read data from MongoDB and export them for the use of fitpot.
    Output files are:
        - pos
        - erg.ref
        - frc.ref
    in ##### directories.
    """
    client,db,col= load_DB(conf)
    os.system("mkdir -p {0}".format(dirname))
    nsmpl= col.count()
    docs= col.find(query)
    i= offset
    for doc in docs:
        i += 1
        savedir= "{0}/{1:05d}".format(dirname,i)
        os.system('mkdir -p '+savedir)
        pos= doc_to_pos(doc,conf)
        pos.write_pmd(savedir+"/pos")
        forces= doc['calculations'][-1]['output']['ionic_steps'][-1]['forces']
        #...Reference energy output
        with open(savedir+'/erg.ref','w') as f:
            f.write(' {0}\n'.format(doc['output']['final_energy']))
        #...Reference forces output
        with open(savedir+'/frc.ref','w') as g:
            g.write(' {0}\n'.format(doc['nsites']))
            for force in forces:
                g.write(' {0:12.7f} {1:12.7f} {2:12.7f}\n'.format(force[0],force[1],force[2]))
        #...DB information output in YAML format
        with open(savedir+'/db_info.yaml','w') as f:
            dbinfo= {}
            dbinfo['db']= db.name
            dbinfo['collection']= col.name
            dbinfo['task_id']= doc['task_id']
            dbinfo['dir_name']= doc['dir_name']
            dbinfo['energy']= doc['output']['final_energy']
            dbinfo['energy_per_atom']= doc['output']['final_energy_per_atom']
            dbinfo['forces']= forces
            unicode_to_ascii(dbinfo)
            yaml.dump(dbinfo,f)
        print '.',


if __name__ == '__main__':

    args= docopt(__doc__)
    confname= args['CONFIG_FILE']
    query= args['--query']
    dirname= args['--dir']
    offset= int(args['--offset'])
    
    query= json.loads(query)
    # print query

    conf= read_db_config(confname)
    main(conf,dirname,query,offset)

