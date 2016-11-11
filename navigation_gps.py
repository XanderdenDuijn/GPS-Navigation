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
for i,line in enumerate(obs_file):
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

obs_nr = int(header[2][:6].strip()) #number of different observation types
# will is use the next line if obs_nr > 8 or 9???
obs_types = []
# iterate over line with steps of 6 chars
# 9 iterations = nr of possible observation types in 1 line
for idx in range(6,59,6):
    obs_types.append(header[2][idx:idx+6].strip())

################################
### READING THE OBSERVATIONS ###
################################
    
i += 1 # go to next line   
data = {}
data['observations'] = []
epoch_nr = 1
while i < len(obs_file):
    #print 'START EPOCH'
    
    cur_line = obs_file[i]
    next_line = obs_file[i+1]

    # read epoch info
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

    sat_names = [] # list with satellite names
    
    #iterate over satellite identifaction with steps of 3 characters
    for idx in range(0,len(sat_id),3):
        sat_names.append(sat_id[idx:idx+3])

    dic = {}
    ### FILL THE DICTIONARY
    dic['epoch'] = epoch_nr
    dic['year'] = epoch[0]
    dic['month'] = epoch[1]
    dic['day'] = epoch[2]
    dic['hour'] = epoch[3]
    dic['min'] = epoch[4]
    dic['sec'] = epoch[5]
    dic['flag'] = flag
    dic['satellite_info'] = {}
    for sat_name in sat_names:
        dic['satellite_info'][sat_name] = {}
   
    cur_line = obs_file[i]
    
    ### REACH OBSERVATIONS

    for obs_idx,j in enumerate(range(i,i+nr_sats)):
        cur_line = obs_file[j][:-1]
        for obstype_idx, obs_val in enumerate(range(0,len(cur_line[:-1]),16)):
            dic['satellite_info'][sat_names[obs_idx]][obs_types[obstype_idx]] = cur_line[obs_val:obs_val+16].strip()   

    data['observations'].append(dic) #append the observation data to the dict
    epoch_nr += 1
    i += nr_sats
    #print 'NEXT EPOCH'
    
# TESTING
##print data['observations'][0]['satellite_info']['G02']['P2']
##for sat in data['observations'][0]['satellite_info']:
##    print sat,data['observations'][0]['satellite_info'][sat]['P2']

####################################
### READING THE NAVIGATION FILE ###
####################################

nav_input = open('brdc0590.11n','r')
nav_file = nav_input.readlines()
nav_input.close()

nav_header = []
i = 0 # i can be seen as the line number
for i,line in enumerate(nav_file):
    if 'ION ALPHA' in line:
        nav_header.append(line)
    elif 'ION BETA' in line:
        nav_header.append(line)
    elif 'LEAP SECONDS' in line:
        nav_header.append(line)
    elif 'END OF HEADER' in line:
        break
#print nav_header

# reading first line of navigation file is not possible to check whether it is a GPS navigation file

i +=1 # got to next line

def GDtoJD(d,m,y,hh,mm,ss):
    #representing time as a real number
    hour = hh+mm/60.0+ss/3600.0

    if m <= 2:
        y = y-1
        m = m+12
    JD = int(365.25*y)+int(30.6001*(m+1))+d+hour/24.0+1720981.5
    return JD


for epoch in data['observations']:
    obs_time = GDtoJD(int(epoch['day']), int(epoch['month']), int(epoch['year']), int(epoch['hour']), int(epoch['min']), float(epoch['sec']))
    #print time
    #print obs_time
    #print int(epoch['day']), int(epoch['month']), int(epoch['year']), int(epoch['hour']), int(epoch['min']), float(epoch['sec'])
    sat_lst = []
    cnt = 0
    # get navigation satellite and epoch data
    # jumping over 8 lines
    for i in range(i, len(nav_file),8):
        line = nav_file[i]

        #split up the line to get epoch data
        y = int(line[3:6])
        m = int(line[6:9])
        d = int(line[9:11])
        hh = int(line[11:14])
        mm = int(line[14:17])
        ss = float(line[17:22])

        #convert Gregorian date to Julian Date
        nav_time = GDtoJD(d,m,y,hh,mm,ss)
        #calculate time difference
        delta_time = abs(obs_time-nav_time)

        #one hour in Julian Date
        hour = 1/24.0

        #find appropriate navigation epoch
        #if the time difference is less than one hour...
        if delta_time < hour:
            nav_sat = int(line[:2])
            print 'NAV SAT', nav_sat
            #iterate over the satellite names of the obervation epoch
            for k in epoch['satellite_info']:
                # it must be a galileo satellite (G)
                # the number of the Galileo satellite and navigation satellite must be equal
                # there are cases that satellite x can appear more than once (for ..59.44 and ..00.00)
                # if that is the case we use the ..59.44 (more accurate??)
                if k[0:1] == 'G' and int(k[1:]) == nav_sat and nav_sat not in sat_lst: 
                    sat_lst.append(nav_sat)
                    print sat_lst
                    epoch['satellite_info'][k]['params'] = {}

                    params_lst = []
                    for j in range(i,i+8):
                        line = nav_file[j]
                        for idx in range(3,len(line)-1,19):
                            if 'D' in line[idx:idx+19]:
                                params_lst.append(line[idx:idx+19].replace('D','e'))
                            else:
                                params_lst.append(line[idx:idx+19])
                    # remove first item in list. We only want the orbital parameters
                    params_lst.pop(0)
                    for idx,char in enumerate('abcdefghijklmnopqrstuvwxyzABCDE'):
                        
                        epoch['satellite_info'][k]['params'][char] = params_lst[idx]
                    
                else:
                    continue


        #print i
    break
    #print delta_time


print data['observations'][0]['satellite_info']['G02']

#for k in data['observations'][0]['satellite_info']['G02']:
    #print k
