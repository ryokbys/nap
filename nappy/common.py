import re

def get_key(v):
    prefix, index = re.match(r'([a-zA-Z]+)_(\d+)', v).groups()
    return prefix, -int(index)

