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
    return dict


if __name__ == '__main__':
    dict=read_input()
    print dict
