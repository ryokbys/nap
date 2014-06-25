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
        else:
            dataset=[]
            for i in range(1,len(data)):
                dataset.append(data[i])
            dict[data[0]]= dataset
    f.close()
    # convert type
    if 'num_samples' in dict:
        dict['num_samples'][0]= int(dict['num_samples'][0])
    if 'num_iteration' in dict:
        dict['num_iteration'][0]= int(dict['num_iteration'][0])
    if 'eps' in dict:
        dict['eps'][0]= float(dict['eps'][0])
    if 'xtol' in dict:
        dict['xtol'][0]= float(dict['xtol'][0])
    if 'gtol' in dict:
        dict['gtol'][0]= float(dict['gtol'][0])
    if 'ftol' in dict:
        dict['ftol'][0]= float(dict['ftol'][0])
    if 'atom_energy' in dict:
        dict['atom_energy'][0]= int(dict['atom_energy'][0])
        dict['atom_energy'][1]= float(dict['atom_energy'][0])
    
    # ga parameters
    if 'ga_num_individuals' in dict:
        dict['ga_num_individuals'][0]= int(dict['ga_num_individuals'][0])
    if 'ga_num_bit' in dict:
        dict['ga_num_bit'][0]= int(dict['ga_num_bit'][0])
    if 'ga_temperature' in dict:
        dict['ga_temperature'][0]= float(dict['ga_temperature'][0])
    if 'ga_murate' in dict:
        dict['ga_murate'][0]= float(dict['ga_murate'][0])
    return dict


if __name__ == '__main__':
    dict=read_input()
    print dict
