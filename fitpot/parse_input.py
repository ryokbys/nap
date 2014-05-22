u"""Parse input file, 'in.fitpot.'

*This code may be somewhat stupid, because it is written only for reading
'in.fitpot' file.*
"""

def read_input(fname='in.fitpot'):
    u""" Read in.fitpot file.
    """
    dict={}
    f= open(fname,'r')
    for line in f.readlines():
        data= line.split()
        # skip if the line is empty or comment line
        if len(data)==0 or \
                line[0]=='#' or \
                line[0]=='!' or \
                line[0]=='%':
            continue
        dict[data[0]]=data[1]
    f.close()
    # convert type
    if 'num_samples' in dict:
        dict['num_samples']= int(dict['num_samples'])
    if 'num_iteration' in dict:
        dict['num_iteration']= int(dict['num_iteration'])
    if 'eps' in dict:
        dict['eps']= float(dict['eps'])
    if 'xtol' in dict:
        dict['xtol']= float(dict['xtol'])
    if 'gtol' in dict:
        dict['gtol']= float(dict['gtol'])
    if 'ftol' in dict:
        dict['ftol']= float(dict['ftol'])
    # ga parameters
    if 'ga_num_individuals' in dict:
        dict['ga_num_individuals']= int(dict['ga_num_individuals'])
    if 'ga_num_bit' in dict:
        dict['ga_num_bit']= int(dict['ga_num_bit'])
    if 'ga_temperature' in dict:
        dict['ga_temperature']= float(dict['ga_temperature'])
    if 'ga_murate' in dict:
        dict['ga_murate']= float(dict['ga_murate'])
    return dict


if __name__ == '__main__':
    dict=read_input()
    print dict
