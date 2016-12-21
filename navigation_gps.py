#-------------------------------------------#
#   POST-PROCESS GPS NAVIGATION SOFTWARE    #
#   BY XANNDER DEN DUIJN AND TONI PRATS     #
#-------------------------------------------#


import math
import numpy as np
from conversions import *
import csv
import pygmaps
import webbrowser
import os
import matplotlib.pyplot as plt
import datetime

####################################
### READING THE OBSERVATION FILE ###
####################################

obs_input = open('89090590.11o','r')
obs_file = obs_input.readlines()
obs_input.close()

##########################
### READING THE HEADER ###
##########################

obs_header = []
obs_types_temp = [] #will be a list of lists
i = 0 # i = line number
for i,line in enumerate(obs_file):
    if 'APPROX POSITION XYZ' in line:
        aprox_xyz = []
        for idx in range(0,42,14):
            coord = float(line[idx:idx+14].strip())
            aprox_xyz.append(coord)
        obs_header.append(aprox_xyz)
    elif 'ANTENNA: DELTA H/E/N' in line:
        delta_hen = []
        for idx in range(0,42,14):
            hen_value = float(line[idx:idx+14].strip())
            delta_hen.append(hen_value)
        obs_header.append(delta_hen)   
    elif '# / TYPES OF OBSERV' in line: #there could be more lines
        # iterate over line with steps of 6 chars
        # 9 iterations = nr of possible observation types in 1 line           
        for idx in range(0,60,6):
            obs_type = line[idx:idx+6].strip()
            obs_types_temp.append(line[idx:idx+6].strip())     
    elif 'TIME OF FIRST OBS' in line:
        first_obs = []
        for idx in range(0,30,6):
            first_obs.append(line[idx:idx+6].strip())
        first_obs.append(line[30:43].strip())
        first_obs.append(line[43:51].strip())
        obs_header.append(first_obs)
    elif 'END OF HEADER' in line:
        break   

#update list of observations from header
obs_types = [] #final list for observation 
obs_nr = obs_types_temp[0] #the first item is the number of satellites
obs_types_temp.pop(0) #remove nr of satellites from obs_types
for obs_type in obs_types_temp:
    if obs_type: #if the observation is not an empty string..
        obs_types.append(obs_type)
obs_header.insert(2,obs_types) #put the list of observations in the correct spot
        
################################
### READING THE OBSERVATIONS ###
################################
    
i += 1 # go to next line   
data = {} #creating the main datastructure
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

    #string with satellite identification
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

    dic = {} #for every epoch
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
        dic['satellite_info'][sat_name] = {} #for every satellite a dictionary is created
   
    cur_line = obs_file[i]
    
    ### REACH OBSERVATIONS

    for obs_idx,j in enumerate(range(i,i+nr_sats)): #internal loop iterating over lines with observation values
        cur_line = obs_file[j][:-1]
        #map satellite, observation types and values --> put in dictionary
        for obstype_idx, obs_val in enumerate(range(0,len(cur_line[:-1]),16)):
            dic['satellite_info'][sat_names[obs_idx]][obs_types[obstype_idx]] = cur_line[obs_val:obs_val+14].strip()     

    data['observations'].append(dic) #append the observation data to the main dict
    epoch_nr += 1
    i += nr_sats #go to next line with epoch information


####################################
### READING THE NAVIGATION FILE ####
####################################

nav_input = open('brdc0590.11n','r')
nav_file = nav_input.readlines()
nav_input.close()

##########################
### READING THE HEADER ###
##########################

nav_header = []
i = 0 # i can be seen as the line number
for i,line in enumerate(nav_file):
    if 'ION ALPHA' in line:
        ion_alpha = []
        for idx in range(2,50,12):
            value = line[idx:idx+12].strip()
            if 'D' in value:
                value = value.replace('D','e')
            value = float(value)
            ion_alpha.append(value)
        nav_header.append(ion_alpha)
    elif 'ION BETA' in line:
        ion_beta = []
        for idx in range(2,50,12):
            value = line[idx:idx+12].strip()
            if 'D' in value:
                value = value.replace('D','e')
            value = float(value)
            ion_beta.append(value)
        nav_header.append(ion_beta)
    elif 'LEAP SECONDS' in line:
        leap_sec = int(line[0:50].strip())
        nav_header.append(leap_sec)
    elif 'END OF HEADER' in line:
        break

i +=1 # got to next line

start_i = i

########################
### READING THE BODY ###
########################

for obs_idx, epoch in enumerate(data['observations']): #epoch is a dict
    print obs_idx   
    
    obs_time = GDtoJD(epoch['day'], epoch['month'], epoch['year'], epoch['hour'], epoch['min'], epoch['sec'])

    # minimum time interval in Julian Date
    hour = 1/24.0
    minute = 1/(24*3600.0)
    
    lst = []
    # read navigation epoch with orbital parameters
    # get navigation satellite and epoch data
    # jumping over 8 lines
    for i in range(start_i, len(nav_file),8):
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

        #first selection. Add items to list in case the delta_time is less than the minimum time interval
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

    #iterate over the list with the correct linenumber      
    for delta_time,i in lst:
        line = nav_file[i] #go to corresponding line
        nav_sat = int(line[:2])
        #print 'NAV SAT', nav_sat
        
        #iterate over the satellite names of the obervation epoch to enrich the data structure
        for sat_name in epoch['satellite_info']:
            #print sat_name
            # it must be a GPS satellite (G)
            if sat_name[0:1] == 'G':
                # the number of the GPS satellite and navigation satellite must be equal
                if int(sat_name[1:]) == nav_sat:
                    epoch['satellite_info'][sat_name]['params'] = {} #create a dictionary for all the parameters
                    params_lst = []
                    for j in range(i,i+8): #internal iteration over lines of specific epoch
                        line = nav_file[j]
                        for idx in range(3,len(line)-1,19):
                            if 'D' in line[idx:idx+19]:
                                params_lst.append(line[idx:idx+19].replace('D','e'))
                            else:
                                params_lst.append(line[idx:idx+19])
                    # remove first item in list. We only want the orbital parameters
                    params_lst.pop(0)

                    #populating the dict with the orbital parameters
                    params = epoch['satellite_info'][sat_name]['params']
            
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
                    #this is not the satellite we are looking for at the moment/or has not been observed
                    continue
            else:
                # It is not a GPS satellite
                continue

###################################################
### COMPUTE SATELLITE COORDINATE FOR EACH EPOCH ###
###################################################

for obs_idx, epoch in enumerate(data['observations']): #epoch is a dict
    #print epoch
    print obs_idx
    
    for sat_name in epoch['satellite_info']: #iterate over satellites=keys
         
        #print
        #print 'SATELLITE NAME', sat_name
        #print
    
        #it has to be a GPS satellite
        if sat_name[0] == 'G':
            try: #it could be that there is no p2 value
                y = data['observations'][obs_idx]['year'] +2000
                m = data['observations'][obs_idx]['month']
                d = data['observations'][obs_idx]['day']
                hh = data['observations'][obs_idx]['hour']
                mm = data['observations'][obs_idx]['min']
                ss = data['observations'][obs_idx]['sec']

                # Transform Gregorian date to Julian date
                JD = GDtoJD(d,m,y,hh,mm,ss)

                # Get the number of seconds in the GPS week
                Trec = JDtosecofweek(JD)
                
                #print 'Trec: ', Trec
                c = 299792458 # speed of light in m/s
                Toe = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['TOE']
                #print 'TOE: ', Toe

                #there are some cases in which there is no p2 value
                p2 = float(data['observations'][obs_idx]['satellite_info'][sat_name]['P2'])

                # Calculate emission epoch
                Temis1 = (Trec - (p2/c))-Toe
                #print 'TEMIS BEFORE: ', Temis1
                if Temis1 > 302400:
                    Temis1 = Temis1 - 604800
                elif Temis1 < -302400:
                    Temis1 = Temis1 + 604800
                #print 'TEMIS1 AFTER', Temis1

                # Compute clock correction
                a0 = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['SV Clock Bias']
                a1 = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['SV Clock Drift']
                a2 = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['SV Clock Drift Rate']
                Tcorr1 = a0 + a1*Temis1 + a2*(Temis1**2)
                #print Tcorr1

                # Calculate emission epoch again, including clock corrections
                Temis2 = (Trec - (p2/c) - Tcorr1) - Toe
                #print Temis2
                if Temis2 > 302400:
                    Temis2 = Temis2 - 604800
                elif Temis2 < -302400:
                    Temis2 = Temis2 + 604800
                #print Temis2

                # Second calculation of the clock correction
                Tcorr2 = a0 + a1*Temis2 + a2*(Temis2**2)
                #print 'TCORR2', Tcorr2

                # Final emission time
                Temis = (Trec - (p2/c) - Tcorr2) - Toe
                if Temis > 302400:
                    Temis = Temis - 604800
                elif Temis < -302400:
                    Temis = Temis + 604800
                #print 'Final emission time', Temis

                ### START COMPUTATION OF SATELLITE COORDINATES ###

                # Semi-Major axis of the satellite orbit
                a = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Sqrt(a)']
                a = a**2

                # Compute mean notion
                delta_n = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Delta n']
                n = math.sqrt(3.986005e+14/a**3) + delta_n
                #print 'MEAN MOTION: ', n

                # Compute mean anomaly
                M0 = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Mo']
                M = M0 + n * Temis
                #print 'MEAN ANOMALY: ', M

                # Eccentric anomaly
                e = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Eccentricity']
                
                # first iteration
                E = M + e * math.sin(M)
                # iteraration
                for i in range(10):
                    E = M + e * math.sin(E)
                #print 'ECCENTRIC ANOMALY', E

                # Calculate true anomaly with arctan2 (to avoid negative angels)
                v = np.arctan2((math.sqrt(1-e**2)*math.sin(E)),(math.cos(E)-e))
                #print 'TRUE ANOMALY', v
            
                # Argument of latitude
                omega = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Omega']
                phi = v + omega
                #print 'ARGUMENT OF LATITUDE', phi

                # Orbital correction terms
                Cus = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Cus']
                Cuc = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Cuc']
                delta_u = Cus * math.sin(2*phi) + Cuc * math.cos(2*phi)
                Crs = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Crs']
                Crc = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Crc']
                delta_r = Crs * math.sin(2*phi) + Crc * math.cos(2*phi)
                Cis = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['CIS']
                Cic = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Cic']
                delta_i = Cis * math.sin(2*phi) + Cic * math.cos(2*phi)
                #print delta_u, delta_r, delta_i

                # Corrected argument of latitude, radius and inclination:
                phi = phi + delta_u
                #print 'PHI', phi
                r = a*(1-e*math.cos(E)) + delta_r
                #print 'r', r
                Io = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['Io']
                IDOT = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['IDOT']
                i = Io + delta_i + IDOT * Temis
                #print 'i', i

                # Position of orbital plane
                xop = r * math.cos(phi)
                #print xop
                yop = r * math.sin(phi)
                #print yop

                #correct longitude of ascending node
                rotation = 7.2921151467e-5
                OMEGA = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['OMEGA']
                OMEGA_DOT = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['OMEGA DOT']
                lon_asc = OMEGA + (OMEGA_DOT - rotation) * Temis - rotation * Toe
                #print lon_asc

                # Final satellite coordinates (ECEF) !!!at emission time!!!
                sat_x = xop * math.cos(lon_asc)-yop*math.cos(i)*math.sin(lon_asc)
                sat_y = xop * math.sin(lon_asc)+yop*math.cos(i)*math.cos(lon_asc)
                sat_z = yop*math.sin(i)
                #print sat_x
                #print sat_y
                #print sat_z

                #add coordinates to main dictionary
                data['observations'][obs_idx]['satellite_info'][sat_name]['X'] = float(sat_x)
                data['observations'][obs_idx]['satellite_info'][sat_name]['Y'] = float(sat_y)
                data['observations'][obs_idx]['satellite_info'][sat_name]['Z'] = float(sat_z)
                
                ### Final corrections

                # Relavistic effect due to the orbit accentricity
                trel = -2 * (math.sqrt(3.986005e+14*a)/c**2) * e * math.sin(E)
                #print 'TREL', trel
                
                # Total group delay(TGD)/satellite instrumental bias
                TGD_L1 = data['observations'][obs_idx]['satellite_info'][sat_name]['params']['TGD']
                TGD_L2 = 1.65*TGD_L1
                Tcorr = Tcorr2 + trel - TGD_L2
                #print 'Tcorr', Tcorr

                #add correction and sec_of_week to main dictionary
                data['observations'][obs_idx]['satellite_info'][sat_name]['Tcorr'] = Tcorr
                data['observations'][obs_idx]['satellite_info'][sat_name]['sec_of_week'] = Trec
            except:
                #there is no p2 value, thus it is not possible to calculate sat coordinates
                data['observations'][obs_idx]['satellite_info'][sat_name]['X'] = 'None'
                data['observations'][obs_idx]['satellite_info'][sat_name]['Y'] = 'None'
                data['observations'][obs_idx]['satellite_info'][sat_name]['Z'] = 'None'
        else:
            #it is not a GPS Satellite
            data['observations'][obs_idx]['satellite_info'][sat_name]['X'] = 'None'
            data['observations'][obs_idx]['satellite_info'][sat_name]['Y'] = 'None'
            data['observations'][obs_idx]['satellite_info'][sat_name]['Z'] = 'None'
        

##################################################
### START CALCULATION OF THE NAVIGATION POINTS ###
##################################################

#create csv output file
output = open('results.csv', 'wb')
writer = csv.writer(output, delimiter=',')
fieldnames = ['Hour', 'Min', 'Sec', 'X', 'Y', 'Z', 'Xerror', 'Yerror', 'Zerror', 'Lat', 'Lon', 'h']
writer.writerow(fieldnames)

### REFERENCE STATION INFORMATION ###

#ECEF APPROXIMATE XYZ COORDINATES OF REFERENCE STATION
approx_x = obs_header[0][0] #x
approx_y = obs_header[0][1] #y
approx_z = obs_header[0][2] #z

# Convert ECEF (XYZ) of reference station to lat, lon, h / ellipsoidal coordinates / geodetic coordinates
approx_lat, approx_lon, approx_h = xyztolatlonh(approx_x,approx_y, approx_z)
#print approx_lat, approx_lon, approx_h

#convert lon, lat coordinates of reference station to radians
approx_lat_rad = approx_lat*math.pi/180.0
approx_lon_rad = approx_lon*math.pi/180.0

solutions = [] #final results
xerrors = []
yerrors = []
zerrors = []

for obs_idx, epoch in enumerate(data['observations']):
    print
    print obs_idx
    print

    if data['observations'][obs_idx]['flag'] == 0:
        print 'EPOCH NUMBER: ', obs_idx
         
        dTrec = 0
        
        ### ITERATE 6 TIMES TO SOLVE NAVIGATION POINT COORDINATES ###
        for itr in range(6):
            #print 
            #print 'ITERATION NUMBER: ', itr
            #print
            
            #initialize the matrices to calc receiver position
            A_lst = []
            R_lst = []
            W_lst = []
 
            #ITERATE OVER EACH SATELLITE
            for sat in data['observations'][obs_idx]['satellite_info']:
                if sat[0] == 'G':
                    print 'SATELLITE: ', sat
                    
                    #print obs_idx, sat
                    sat_x = data['observations'][obs_idx]['satellite_info'][sat]['X']
                    #print sat_x
                    sat_y = data['observations'][obs_idx]['satellite_info'][sat]['Y']
                    #print sat_y
                    sat_z = data['observations'][obs_idx]['satellite_info'][sat]['Z']
                    #print sat_z

                    # the xyz coordinates of the satellite should be valid
                    if sat_x != 'None' and sat_y != 'None' and sat_z != 'None':

                        ############################################
                        ### CALCULATION OF SATELLITE COORDINATES ###
                        ############################################
 
                        # calculate the distance between the approx station point and the sat (using ECEF coords)
                        d = math.sqrt((sat_x-approx_x)**2+(sat_y-approx_y)**2+(sat_z-approx_z)**2)
                        
                        traveltime = d/c
                        #print traveltime
                        wt = traveltime * rotation

                        ### from xyz of sat at emission time to xyz at reception time
                        
                        R3 = np.array([[math.cos(traveltime*rotation),math.sin(traveltime*rotation),0],
                                      [-math.sin(traveltime*rotation),math.cos(traveltime*rotation),0],
                                      [0,0,1]])
                        sat_x_rot = np.dot(R3,np.array([[sat_x],[sat_y],[sat_z]]))[0]
                        
                        #print float(sat_x_rot)
                        sat_y_rot = np.dot(R3,np.array([[sat_x],[sat_y],[sat_z]]))[1]
                        #print float(sat_y_rot)
                        sat_z_rot = np.dot(R3,np.array([[sat_x],[sat_y],[sat_z]]))[2]
                        #print float(sat_z_rot)

                        # delta distances to satellite
                        dx = [float(sat_x_rot - approx_x)] #x difference with station as array
                        #print dx
                        dy = [float(sat_y_rot - approx_y)] #y difference with station as array
                        #print dy
                        dz = [float(sat_z_rot - approx_z)] #z difference with station as array
                        #print dz
                        dX = np.array([dx,dy,dz])

                        # Convert ECEF (XYZ) of the sats to lat, lon, h
                        lat, lon, h = xyztolatlonh(sat_x_rot,sat_y_rot,sat_z_rot)

                        # Convert lon, lat satellite coordinates to radians
                        lat_rad, lon_rad = lat * math.pi/180.0, lon * math.pi/180.0

                        #matrix/formula to tranform coordinate systems
                        R31 = np.array([[-math.sin(approx_lat_rad)*math.cos(approx_lon_rad), -sin(approx_lat_rad)*sin(approx_lon_rad),cos(approx_lat_rad)],
                                       [-math.sin(approx_lon_rad),math.cos(approx_lon_rad),0],
                                       [math.cos(approx_lat_rad)*math.cos(approx_lon_rad), math.sin(approx_lon_rad)*math.cos(approx_lat_rad),math.sin(approx_lat_rad)]])

                        NEh = np.dot(R31,dX) 
                        NEh = [NEh[0][0],NEh[1][0],NEh[2][0]]
                        #print NEh
                        El_rad, Az_rad = Elev_Az(NEh) #in radians
                        #print El_rad, Az_rad
                        El_degr = El_rad * (180.0/math.pi) #convert to degrees
                        Az_degr = Az_rad * (180.0/math.pi) #convert to degrees
                        #print El_degr, Az_degr
                        
                        ### IN CASE THE SATELLITE ANGLE IS MORE THAN 10 DEGREES...
                        if El_degr > 10:

                            ### TROPOSPHERIC DELAY ###
                            
                            z_degr = 90 - El_degr # degrees
                            z_rad = z_degr * math.pi/180

                            #h_ECEF = float(sat_z_rot) #elipsoidal height
                            h_ECEF = approx_h
                            
                            p = 1013.25*(1-0.000065*h_ECEF)**(5.225)
                            #print p
                            T = 291.15-0.0065*h_ECEF
                            #print T
                            H = 50*math.exp(-0.0006396*h_ECEF)
                            #print H
                            
                            e = (H*0.01)*math.exp(-37.2465+0.213166*T-(0.000256908*T**2))
                            #print e
                            
                            #find smallest diff with height
                            #find corresponding B(mb) parameter

                            B_dic = {0.0:1.156,
                                     0.5:1.079,
                                     1.0:1.006,
                                     1.5:0.938,
                                     2.0:0.874,
                                     2.5:0.813,
                                     3.0:0.757,
                                     4.0:0.654,
                                     5.0:0.563}
                            B = B_dic[min(B_dic, key=lambda x:abs(x-approx_h/1000))]
                            #print B
                            
                            d_tropo=(0.002277/math.cos(z_rad))*(p+((1255/T)+0.05)*e-B*(math.tan(z_rad))**2)
                            #print d_tropo

                            ### IONOSPHERIC DELAY ###

                            #1 calcualte earth centered angle
                            eca = 0.0137/((El_rad/math.pi)+0.11) - 0.022 #semicircles
                            #print eca
                        
                            #2 compute the latitude of the ionospheric Pierce Point (IPP)
                            # lat corresponds with the users position
                            # A is azimuth of the satellite in radians
                            
                            ipp_lat = approx_lat_rad/math.pi + eca*math.cos(Az_rad)
                            if ipp_lat > 0.416:
                                ipp_lat = 0.416
                            elif ipp_lat < -0.416:
                                ipp_lat = -0.416
                            #print ipp_lat
                        
                            #3 compute the longitude of ionospheric Pierce Point (IPP)
                            #lon corresponds with the users position
                            ipp_lon = approx_lon_rad/math.pi + eca*math.sin(Az_rad)/math.cos(ipp_lat*math.pi)
                            #print ipp_lon

                            #4 find the geomagnetic latitude of the IPP
                            geom_lat = ipp_lat + 0.064*math.cos((ipp_lon-1.617)*math.pi)
                            #print geom_lat

                            #5 find the local time at the IPP
                            ipp_time = 43200 * ipp_lon + Trec
                            if ipp_time >= 86400:
                                ipp_time -= 86400
                            elif ipp_time < 0:
                                ipp_time += 86400
                            #print ipp_time

                            #6 compute the amplitude of the ionospheric delay in seconds
                            Ai = []
                            for i in range(4):
                                Ai.append(nav_header[0][i]*geom_lat**i)
                            
                            Ai = sum(Ai)
                            if Ai < 0:
                                Ai = 0
                            #print Ai

                            #7 compute the period of the ionospheric delay in seconds
                            Pi = []
                            for i in range(4):
                                Pi.append(nav_header[1][i]*geom_lat**i)
                            
                            Pi = sum(Pi)
                            #print Pi
                            if Pi < 72000:
                                Ai = 72000
                            #print Pi

                            #8 compute the phase of the ionospheric delay in radians!!!
                            Xi = 2*math.pi*(ipp_time-50400)/Pi
                            #print Xi

                            #9 compute the slant factor (elevation E in semicircles)
                            F = 1.0 + 16.0*(0.53-El_rad/math.pi)**3
                            #print F

                            #10 compute the ionospheric delay time (in seconds) for L1

                            if Xi < 1.57:
                                IL1 = (5e-9+Ai*(1 - Xi**2/2 + Xi**4/24))*F
                            else:
                                IL1 = 5e-9 *F

                            #transform to meters
                            IL1 = IL1 * c
                            #print IL1

                            #11 compute the ionospheric time delay for L2 (in meters)
                            IL2 = 1.65 * IL1
                            #print IL2

                            #actual distance
                            d_ion = IL2

                            ### SYSTEM EQUATION ###
                            dist = math.sqrt((sat_x_rot-approx_x)**2+(sat_y_rot-approx_y)**2+(sat_z_rot-approx_z)**2)

                            
                            A = [float(-(sat_x_rot-approx_x)/dist),
                                 float(-(sat_y_rot-approx_y)/dist),
                                 float(-(sat_z_rot-approx_z)/dist),
                                 1]

                            Tcorr = data['observations'][obs_idx]['satellite_info'][sat]['Tcorr']
                            p2 = float(data['observations'][obs_idx]['satellite_info'][sat]['P2'])

                            R = [float(p2 - dist - dTrec + c*Tcorr - d_tropo - d_ion)]
     
                            W = [float(math.sin(El_rad)**2/(1.5*0.3)**2)]

                            A_lst.append(A)
                            R_lst.append(R)
                            W_lst.append(W)
                            
                        else:
                            print 'elevation is less than 10 degrees'
                            #elevation is less than 10 degrees
                            continue
                    else:
                        print 'satellite coordinate unknown'
                else:
                    print 'this is not a gps satellite'
                    continue
            print 'ITERATED OVER ALL SUITABLE SATELLITES'

            # calculate receivers position when there are at keast 4 proper sats on the sky
            if len(A_lst)>=4: #we use A_lst here, but we could use R_lst and W_lst as well
                ### LEAST SQUARED SOLUTIONS ###
                #print 'yes we can calculate the receivers position'

                #print A_lst
                #print len(A_lst)
                A_matrix = np.matrix(A_lst) #size: nr of sats * 4
                #print A_matrix
                
                R_matrix = np.matrix(R_lst) #size: nr of sats * 1
                #print R_matrix

                W_matrix = np.matrix(W_lst) #size: nr of sats * 1
                W_matrix = np.diagflat(W_matrix) #size: nr of sats * nr of sats (squared matrix)
                #print W_matrix

                k1 = np.dot(A_matrix.transpose(),W_matrix) #(4*nr of sats) * (nr of sats*nr of sats) --> 4 * nr of sats
                #print k1  

                s = np.linalg.inv(np.dot(k1,A_matrix)) #size: (4*nr of sats) * (nr of sats*4) --> 4*4
                #print s
                
                k2 = np.dot(k1,R_matrix) #size: (4*nr of sats) * (nr of sats*1) --> 4*1
                #print k2
                
                xfin = np.dot(s,k2) #size: (4*4) * (4*1) --> 4*1
                #print xfin
                
                dX_coord = xfin[0]
                #print dX_coord
                dY_coord = xfin[1]
                #print dY_coord
                dZ_coord = xfin[2]
                #print dZ_coord
                dT = xfin[3]
                #print dT

                #residuals
                res = np.dot(A_matrix,xfin)-R_matrix #size: (nr of sats*4) * (4*1) --> nr of sats * 1
                #print res

                approx_x = float(approx_x + dX_coord)
                #print 'x coordinate: ', approx_x
                approx_y = float(approx_y + dY_coord)
                #print 'y coordinate: ', approx_y
                approx_z = float(approx_z + dZ_coord)
                #print 'z coordinate: ', approx_z
                dTrec = dT
                #print dTrec
                   
            else:
                print 'There are not enough suitable satellites to calc the receivers position'            
          
        print 'Dame los coordenadas'

        print approx_x, approx_y, approx_z
        lat,lon,h = xyztolatlonh(approx_x,approx_y, approx_z)
        print lat, lon, h
        
        #compute errors
        (m,n)=A_matrix.shape #size: nr of sats * 4
        desv1=np.dot(res.transpose(),W_matrix) #size: (1*nr of sats) * (nr of sats*nr of sats) --> 1*nr of sats
        desv2=np.dot(desv1,res) #size: (1*nr of sats) * (nr of sats*1) --> 1*1
        des=desv2/(m-n) #1/(nr of sats-4)

        desx=math.sqrt(des*s[0,0])
        xerrors.append(desx)
        desy=math.sqrt(des*s[1,1])
        yerrors.append(desy)
        desz=math.sqrt(des*s[2,2])
        zerrors.append(desz)
        
        print desx,desy,desz

        y = data['observations'][obs_idx]['year']
        m = data['observations'][obs_idx]['month']
        d = data['observations'][obs_idx]['day']
        hh = data['observations'][obs_idx]['hour']
        mm = data['observations'][obs_idx]['min']
        ss = data['observations'][obs_idx]['sec']

        solutions.append([y,m,d,hh,mm,ss,approx_x,approx_y,approx_z,desx,desy,desz,lat,lon,h])
        
        #break
    else:
        print 'flag is not 0, this this observation/epoch is not valid'
        #flag is not 0
        continue
    #break


################################################################
### write solutions to file, plot and map navigations points ###
################################################################

mean_xerror = np.mean(xerrors)
mean_yerror = np.mean(yerrors)
mean_zerror = np.mean(zerrors)

mymap = pygmaps.maps(lat, lon, 13)

final_x = []
final_y = []
final_lat = []
final_lon = []
final_h = []
final_epoch = []

final_hfrac = []
final_time = []

for idx, solution in enumerate(solutions):
    #check if error is less than 1.5 * mean error
    year = solution[0]
    month = solution[1]
    day = solution[2]
    hh = solution[3]
    mm = solution[4]
    ss = solution[5]
    sec,dec = divmod(ss,1)
    sec = int(sec)
    dec = int(dec*1e+6)
    
    x = solution[6]
    y = solution[7]
    z = solution[8]
    xerror = solution[9]
    yerror = solution[10]
    zerror = solution[11]
    lat = solution[12]
    lon = solution[13]
    h = solution[14]
    #only use the ones with a error less than 1.5*mean_error
    if xerror < (1.5*mean_xerror) and yerror < (1.5*mean_yerror) and zerror < (1.5*mean_zerror):
        #lists for xy plot
        final_x.append(x)
        final_y.append(y)
        #lists for latlon plot
        final_lat.append(lat)
        final_lon.append(lon)
        
        #lists for height plot
        #this works
        final_h.append(h)
        final_epoch.append(idx)
        #this does not work
        final_time.append(datetime.datetime(year,month,day,hh,mm,sec,dec)) #somehow this does not work...

        #write to csv file
        writer.writerow([hh,mm,ss,x,y,z,xerror,yerror,zerror,lat,lon,h])

        #map circles on google map
        mymap.addradpoint(lat, lon, 3, "#FF0000")

output.close()

### plots

#height plot

plt.plot(final_epoch, final_h) 
plt.xlabel('epoch')
plt.ylabel('height(m)')
plt.title('height')
plt.grid(True)
plt.savefig("height.png")
#plt.show()

###plot
##plt.plot(final_time, final_h)
###plt.gcf().autofmt_xdate()
##plt.xlabel('time')
##plt.ylabel('h')
##plt.title('height')
##plt.grid(True)
##plt.savefig("height.png")
###plt.show()

#xy plot
plt.plot(final_x, final_y)
range_x = max(final_x)-min(final_x)
range_y = max(final_y)-min(final_y)
axes = plt.gca()
axes.set_xlim([min(final_x)-0.1*range_x,max(final_x)+0.1*range_x])
axes.set_ylim([min(final_y)-0.1*range_y,max(final_y)+0.1*range_y])
#plt.gca().invert_yaxis()
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title('x y')
plt.grid(True)
plt.savefig("xy.png")
#plt.show()

#latlon plot
plt.plot(final_lon, final_lat) 
axes = plt.gca()
range_x = max(final_lon)-min(final_lon)
range_y = max(final_lat)-min(final_lat)
axes.set_xlim([min(final_lon)-0.1*range_x,max(final_lon)+0.1*range_x])
axes.set_ylim([min(final_lat)-0.1*range_y,max(final_lat)+0.1*range_y])
plt.xlabel('lon')
plt.ylabel('lat')
plt.title('lat lon')
plt.grid(True)
plt.savefig("latlon.png")
#plt.plot()


### create map

fname = 'solutions'
#create map file
mymap.draw('./'+fname+'.html')

abspath = os.path.abspath(fname+'.html')
    
# open map in browser
webbrowser.open('file:///'+abspath,2,True)




