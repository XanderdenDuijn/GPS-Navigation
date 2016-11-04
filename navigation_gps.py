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

for i in range(cnt,len(obs_file)):
    # we assume Galileo sats are used
    # what if only Glonass sats are used??
    nr_sats = obs_file[i][30:obs_file[i].find('G')]
    if nr_sats > 12:
        print 'more than one line is used'
    break
    
    

