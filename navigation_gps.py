import math
import numpy

####################################
### READING THE OBSERVATION FILE ###
####################################

obs_input = open('89090590.11o','r')
obs_file = obs_input.readlines()
obs_input.close()

##########################
### READING THE HEADER ###
##########################

header = []
i = 0 # i can be seen as the line number
for line in obs_file:
    i += 1
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

################################
### READING THE OBSERVATIONS ###
################################
    
data = {}
data['observations'] = []
while i < len(obs_file): #!!!MAKE SURE YOU DON'T PRINT THE DATA!!!
    print 'START'
    dic = {}
    cur_line = obs_file[i]
    next_line = obs_file[i+1]

    epoch_nr = 1

    #read epoch info
    epoch = cur_line[0:26].split() # a list with epoch data
    flag = cur_line[28:29]

    # find first letter in line
    for char in cur_line:
        if char in 'GR':
            char_idx = cur_line.find(char)
            break

    nr_sats = int(cur_line[30:char_idx])
    if nr_sats > 12:
        sat_id = cur_line[char_idx:-1]+next_line[char_idx:-1]
        i += 2
    else:
        sat_id = cur_line[char_idx:]
        i += 1

    sat_names = [] # this list will contain the satellite names
    
    #iterate over satellite identifaction with steps of 3 characters
    for idx in range(0,len(sat_id),3):
        sat_names.append(sat_id[idx:idx+3])

    ### FILL THE DICTIONARY
    dic['epoch'] = epoch_nr
    dic['year'] = epoch[0]
    dic['month'] = epoch[1]
    dic['day'] = epoch[2]
    dic['hour'] = epoch[3]
    dic['min'] = epoch[4]
    dic['flag'] = flag
    dic['satellite_info'] = {}
    for sat_name in sat_names:
        dic['satellite_info'][sat_name] = []

    #read epoch observations
    for cnt,i in enumerate(range(i,i+nr_sats)):
        dic['satellite_info'][sat_names[cnt]].append(obs_file[i][:-1])
    
    epoch_nr += 1
    i += 1
    data['observations'].append(dic)
    print 'NEXT EPOCH'
    
    
print len(data['observations'])

