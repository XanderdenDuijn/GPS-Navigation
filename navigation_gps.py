import math
import numpy
from conversions import *

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
obs_types = []
i = 0 # i = line number
for i,line in enumerate(obs_file):
    if 'APPROX POSITION XYZ' in line:
        header.append(line)
    elif 'ANTENNA: DELTA H/E/N' in line:
        header.append(line)
    elif '# / TYPES OF OBSERV' in line:
        # iterate over line with steps of 6 chars, starting at 6
        # 9 iterations = nr of possible observation types in 1 line
        for idx in range(6,60,6):  
            obs_types.append(line[idx:idx+6].strip())
    elif 'TIME OF FIRST OBS' in line:
        header.append(line)
    elif 'END OF HEADER' in line:
        break   

obs_nr = int(header[2][:6].strip()) #number of different observation types
header.insert(2,obs_types)

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
    dic['year'] = int(epoch[0])
    dic['month'] = int(epoch[1])
    dic['day'] = int(epoch[2])
    dic['hour'] = int(epoch[3])
    dic['min'] = int(epoch[4])
    dic['sec'] = float(epoch[5])
    dic['flag'] = int(flag)
    dic['satellite_info'] = {}
    for sat_name in sat_names:
        dic['satellite_info'][sat_name] = {}
   
    cur_line = obs_file[i]
    
    ### REACH OBSERVATIONS

    for obs_idx,j in enumerate(range(i,i+nr_sats)):
        cur_line = obs_file[j][:-1]
        for obstype_idx, obs_val in enumerate(range(0,len(cur_line[:-1]),16)):
            dic['satellite_info'][sat_names[obs_idx]][obs_types[obstype_idx]] = cur_line[obs_val:obs_val+14].strip()     

    data['observations'].append(dic) #append the observation data to the dict
    epoch_nr += 1
    i += nr_sats


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

i +=1 # got to next line

for epoch in data['observations']:
    obs_time = GDtoJD(epoch['day'], epoch['month'], epoch['year'], epoch['hour'], epoch['min'], epoch['sec'])

    # minimum time interval in Julian Date
    hour = 1/24.0
    minute = 1/(24*3600.0)
    lst = []
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

        #make first selection
        if delta_time < hour + minute:
            item = delta_time, i
            lst.append(item)

    # remove the duplicate satellites with a smaller delta_time
    for delta_time1, i1 in lst:
        line1 = nav_file[i1]
        nav_sat1 = int(line1[:2])
        for delta_time2, i2 in lst:
            line2 = nav_file[i2]
            nav_sat2 = int(line2[:2])
            if nav_sat1 == nav_sat2 and delta_time2 < delta_time1:
                item = delta_time1,i1      
                lst.remove(item)
            
    for delta_time,i in lst:
        line = nav_file[i]
        nav_sat = int(line[:2])
        print 'NAV SAT', nav_sat
        #iterate over the satellite names of the obervation epoch
        for k in epoch['satellite_info']:
            # it must be a GPS satellite (G)
            # the number of the GPS satellite and navigation satellite must be equal
            if k[0:1] == 'G' and int(k[1:]) == nav_sat:
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

                #populating the dict with the orbital parameters
                params = epoch['satellite_info'][k]['params']
        
                params['SV Clock Bias'] = float(params_lst[0].strip())
                params['SV Clock Drift'] = float(params_lst[1].strip())
                params['SV Clock Drift Rate'] = float(params_lst[2].strip())
                params['IODE'] = float(params_lst[3].strip())
                params['Crs'] = float(params_lst[4].strip())
                params['Delta n'] = float(params_lst[5].strip())
                params['Mo'] = float(params_lst[6].strip())
                params['Cuc'] = float(params_lst[7].strip())
                params['Eccentricity'] = float(params_lst[8].strip())
                params['Cus'] = float(params_lst[9].strip())
                params['Sqrt(a)'] = float(params_lst[10].strip())
                params['TOE'] = float(params_lst[11].strip())
                params['Cic'] = float(params_lst[12].strip())
                params['OMEGA'] = float(params_lst[13].strip())
                params['CIS'] = float(params_lst[14].strip())
                params['Io'] = float(params_lst[15].strip())
                params['Crc'] = float(params_lst[16].strip())
                params['Omega'] = float(params_lst[17].strip())
                params['OMEGA DOT'] = float(params_lst[18].strip())
                params['IDOT'] = float(params_lst[19].strip())
                params['L2 Codes Channel'] = float(params_lst[20].strip())
                params['GPS Week'] = float(params_lst[21].strip())
                params['L2 P Data Flag'] = float(params_lst[22].strip())
                params['SV Accuracy'] = float(params_lst[23].strip())
                params['SV Health'] = float(params_lst[24].strip())
                params['TGD'] = float(params_lst[25].strip())
                params['IODC'] = float(params_lst[26].strip())
                params['Transmission Time'] = float(params_lst[27].strip())
                params['Fit Interval'] = float(params_lst[28].strip())
            else:
                continue
    break #this break is here for testing purposes. Only the first epoch will be used.


for sat_name in data['observations'][0]['satellite_info']:
    if sat_name[0] == 'G':
        print sat_name
        y = data['observations'][0]['year'] +2000
        m = data['observations'][0]['month']
        d = data['observations'][0]['day']
        hh = data['observations'][0]['hour']
        mm = data['observations'][0]['min']
        ss = data['observations'][0]['sec']
        print hh, mm, ss
        JD = GDtoJD(d,m,y,hh,mm,ss)
        print JD
        Trec = JDtoGPS(JD)
        Trec = 135110.000007
        #sec_of_week = 135110.000007
        print 'Trec: ', Trec
        c = 299792458
        Toe = data['observations'][0]['satellite_info'][sat_name]['params']['TOE']
        print 'TOE: ', Toe
        p2 = float(data['observations'][0]['satellite_info'][sat_name]['P2'])
        print 'P2: ', p2
        #TEMIS should be 135109.927325
        print p2/c
        Temis1 = (Trec - (p2/c))-Toe
        print 'TEMIS BEFORE: ', Temis1
        if Temis1 > 302400:
            Temis1 = Temis1 - 604800
        else:
            Temis1 = Temis1 + 604800
        print 'Temis1', Temis1
        a0 = data['observations'][0]['satellite_info'][sat_name]['params']['SV Clock Bias']
        a1 = data['observations'][0]['satellite_info'][sat_name]['params']['SV Clock Drift']
        a2 = data['observations'][0]['satellite_info'][sat_name]['params']['SV Clock Drift Rate']
        Tcorr1 = a0 + a1*Temis1 + a2*(Temis1**2)
        print Tcorr1

        Temis2 = (Trec - (p2/c) - Tcorr1) - Toe
        print Temis2
        if Temis2 > 302400:
            Temis2 = Temis2 - 604800
        else:
            Temis2 = Temis2 + 604800
        print Temis2

        Tcorr2 = a0 + a1*Temis2 + a2*(Temis2**2)

        Temis = (Trec - (p2/c) - Tcorr2) - Toe
        if Temis > 302400:
            Temis = Temis - 604800
        else:
            Temis = Temis + 604800
        print 'Final emission time', Temis
        a = data['observations'][0]['satellite_info'][sat_name]['params']['Sqrt(a)']
        a = a**2

        # compute mean notion
        delta_n = data['observations'][0]['satellite_info'][sat_name]['params']['Delta n']
        n = math.sqrt(3.986005e+14/a**2) + delta_n

        # compute mean anomaly
        Mo = data['observations'][0]['satellite_info'][sat_name]['params']['Mo']
        M = Mo + n * Temis

        # eccentric anomaly
        e = data['observations'][0]['satellite_info'][sat_name]['params']['Eccentricity']
        # first iteration
        E = M + e * math.sin(M)
        for i in range(10):
            E = M + e * math.sin(E)
            print E

        #calculate true anomaly

        v = math.atan((math.sqrt(1-e**2)*math.sin(E))/(math.cos(E)-e))
    
        #argument of latitude
        omega = data['observations'][0]['satellite_info'][sat_name]['params']['Omega']
        arglat = v + omega

        #orbital correction terms
        Cus = data['observations'][0]['satellite_info'][sat_name]['params']['Cus']
        Cuc = data['observations'][0]['satellite_info'][sat_name]['params']['Cuc']
        delta_u = Cus * sin(2*arglat) + Cuc * cos(2*arglat)
        Crs = data['observations'][0]['satellite_info'][sat_name]['params']['Crs']
        Crc = data['observations'][0]['satellite_info'][sat_name]['params']['Crc']
        delta_r = Crs * sin(2*arglat) + Crc * cos(2*arglat)
        Cis = data['observations'][0]['satellite_info'][sat_name]['params']['Cis']
        Cic = data['observations'][0]['satellite_info'][sat_name]['params']['Cic']
        delta_i = Cis * sin(2*arglat) + Cic * cos(2*arglat)

        #argument of latitude, radius and inclination:
        arglat = arglat + delta_u
        r = a*(1-e*cos(E)) + delta_r
        Io = data['observations'][0]['satellite_info'][sat_name]['params']['Io']
        IDOT = data['observations'][0]['satellite_info'][sat_name]['params']['IDOT']
        i = Io + delta_i + IDOT * Temis

        #position of orbital plane
        xop = r * cos(arglat)
        yop = r * sin(arglat)

        #correct longitude of ascending node
        rotation = 7.2921151467e-5
        OMEGA = data['observations'][0]['satellite_info'][sat_name]['params']['OMEGA']
        OMEGA_DOT = data['observations'][0]['satellite_info'][sat_name]['params']['OMEGA DOT']
        lon_asc = OMEGA + (OMEGA_DOT - rotation) * Temis - OMEGA_DOT * Toe

        # final satellite coordinates
        

        break
    else:
        continue

#print data['observations'][0]['satellite_info']['G02']

#for k in data['observations'][0]['satellite_info']['G02']:
    #print k
