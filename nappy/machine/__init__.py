__version__ = "170130"
__authors__ = "Ryo KOBAYASHI"
__email__   = "kobayashi.ryo@nitech.ac.jp"

class Machine():
    """
    Parent class of any other machine classes.
    """
    QUEUES = {
        'batch':  {'num_nodes': 12, 'default_sec': 86400, 'limit_sec':  604800},
        'default':{'num_nodes': 12, 'default_sec': 86400, 'limit_sec':  604800},
    }
    SCHEDULER = 'pbs'
    
    def __init__(self):
        pass
        

class FUJITSU_FX100(Machine):
    """
    Class for the specific super computer at Nagoya University, Fujitsu FX100.
    """

    # Resource groups and their nums of nodes and time limits
    QUEUES = {
        'fx-debug': {'num_nodes': 32, 'default_sec': 3600, 'limit_sec':  3600},
        'fx-small': {'num_nodes': 16, 'default_sec':86400, 'limit_sec':604800},
        'fx-middle':{'num_nodes': 96, 'default_sec':86400, 'limit_sec':259200},
        'fx-large': {'num_nodes':192, 'default_sec':86400, 'limit_sec':259200},
        'fx-xlarge':{'num_nodes':864, 'default_sec':86400, 'limit_sec': 86400},
    }

    SCHEDULER = 'fujitsu'


class MIKE(Machine):
    """
    Class for the specific machine in Ogata Lab. at NITech.ac.jp
    """

    SCHEDULER = 'pbs'

    #...Currently it is hard to specify nodes to use, so this queues are
    #...like the greatest commond divider
    QUEUES = {
        'batch':  {'num_nodes': 12, 'default_sec': 86400, 'limit_sec':  604800},
        'default':{'num_nodes': 12, 'default_sec': 86400, 'limit_sec':  604800},
    }
