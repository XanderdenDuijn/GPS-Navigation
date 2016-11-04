import math
import numpy

obs_input = open('89090590.11o','r')
obs_file = obs_input.readlines()
obs_input.close()

##aprox_xyz = [line for line in obs_file if 'APPROX POSITION XYZ' in line]
##print aprox_xyz

header = []
cnt = 0

#reading the header
for line in obs_file:
    cnt += 1
    if 'APPROX POSITION XYZ' in line:
        header.append(line)
    elif 'ANTENNA: DELTA H/E/N' in line:
        header.append(line)
    elif '# / TYPES OF OBSERV' in line:
        header.append(line)
    elif 'TIME OF FIRST OBS' in line:
        header.append(line)
    elif 'END OF HEADER' in line:
        break   

#print header
    
#read all observation info
for i in range(cnt,len(obs_file)):
    epoch = 1
    line = obs_file[i]
    next_line = obs_file[i+1]
    # we assume Galileo sats are used
    # what if only Glonass sats are used??
    nr_sats = int(line[30:line.find('G')])
    print nr_sats
    # get epoch
    epoch = line[0:26]
    print epoch
    flag = line[28:29]
    print flag
    if nr_sats > 12:
        sat_id = line[32:-1]+next_line[32:-1]
        print sat_id
        # read obervations for every epoch
        for i in range(i+2,(i+2+nr_sats)):
            print line[:-1]

    # in case nr sats <= 12
    break
    
    

