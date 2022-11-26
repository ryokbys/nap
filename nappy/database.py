#!/usr/bin/env python
"""
This script provides functions to use a nappy database with the file, nappydb.yaml.
Currently, the DB will not function with nested dictionaries.

Usage:
  database.py store [options]
  database.py get [options]
  database.py del [options]
  database.py show

Options:
  -h, --help   Show this message and exit.
  --key KEY    Key. [defualt: None]
  --value VAL  Value. [default: None]
"""
from docopt import docopt
import yaml

__author__ = "RYO KOBAYASHI"
__version__ = "221122"

_dbname = "nappydb.yaml"

class NappyDB(object):
    """This class treats data that can be loaded or stored in nappy modules.
    """

    def __init__(self):
        """Load {0:s} in the working directory or create a new db object.
        """.format(_dbname)

        try:
            self.load()
        except:
            self.db = {}
        return None

    def __repr__(self):
        txt = yaml.dump(self.db, sort_keys=False, default_flow_style=False)
        return txt

    def load(self):
        with open(_dbname,'r') as f:
            self.db = yaml.safe_load(f)
        return None

    def dump(self):
        """This overwrites {0:s} file.
        """.format(_dbname)
        with open(_dbname,'w') as f:
            yaml.dump(self.db,f)
        return None

    def has(self,key,db=None):
        """Look for the given key in the DB not only 1st layer of the DB,
        but also for nested dictionaries as well.
        """
        ans = False
        if db == None:
            db = self.db
        if type(db) is dict:
            for k,v in db.items():
                if k == key:
                    return True
                if type(v) is dict:
                    ans = self.has(key,db=v)
                    if ans:
                        return ans
        elif type(db) is list:
            for d in db:
                ans = self.has(key,db=d)
                if ans:
                    return ans
        return ans

    def __iadd__(self,db):
        """Add the given DB to the current DB.
        If some key-value pairs in the given DB exist in the current DB,
        they will be overwritten.
        """
        for k,v in db.items():
            self.db[k] = v
        return self

    def __getitem__(self,key):
        return self.db.get(key)

    def __setitem__(self,key,val):
        self.db[key] = val
        return None

    def __delitem__(self,key):
        del self.db[key]
        return None


def store(key,val):
    ndb = NappyDB()
    item = {key: val}
    ndb += item
    ndb.dump()
    return None

def get(key):
    ndb = NappyDB()
    if not ndb.has(key):
        raise ValueError('There is no such key stored in the NappyDB.')
    return ndb[key]

def main():
    args = docopt(__doc__)
    key = args['--key']
    if key == 'None':
        key = None
    
    if args['get']:
        val = get(key)
        if 0.1 < abs(val) < 1000.0:
            print(f' {key:s} = {val:8.3f}')
        else:
            print(f' {key:s} = {val:12.4e}')
    elif args['store']:
        val = float(args['--value'])  # assuming all the values are float
        store(key,val)
    elif args['show']:
        ndb = NappyDB()
        print(ndb)
    elif args['del']:
        ndb = NappyDB()
        if not ndb.has(key):
            raise ValueError('There is no such key stored in the NappyDB.')
        else:
            val = ndb[key]
            del ndb[key]
            print(f' Deleted {key}: {val}')
            ndb.dump()
    return None

if __name__ == "__main__":
    main()
