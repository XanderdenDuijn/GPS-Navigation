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

#p2 = 
obs_nr = int(header[2][:6].strip()) #number of different observation types
# will is use the next line if obs_nr > 8 or 9???
obs_types = []
# iterate over line with steps of 6 chars
# 9 iterations = nr of possible observation types in 1 line
for col in range(6,59,6):
    obs_types.append(header[2][col:col+6].strip())

################################
### READING THE OBSERVATIONS ###
################################
    
data = {}
data['observations'] = []
while i < 91: #!!!MAKE SURE YOU DON'T PRINT THE DATA!!!
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
        if char in 'GRSE':
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
        dic['satellite_info'][sat_name] = {}

    cur_line = obs_file[i+1]
        
    #read epoch observations
    for cnt,i in enumerate(range(i,i+nr_sats)):
        cur_line = obs_file[i]
        # get each observation and add it to the dict
        for idx,col in enumerate(range(0,len(cur_line[:-1]),16)):          
            dic['satellite_info'][sat_names[cnt]][obs_types[idx]] = cur_line[col:col+16].strip()

    data['observations'].append(dic)
    i += 1 # go to nect line
    epoch_nr += 1
    print 'NEXT EPOCH'
    
    
for sat in data['observations'][0]['satellite_info']:
    print data['observations'][0]['satellite_info'][sat]['P2']

