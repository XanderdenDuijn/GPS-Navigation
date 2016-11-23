from pyproj import Proj, transform
import math
import numpy as np
from conversions import *
import utm

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
        temp_ion_alpha = []
        temp_ion_alpha.append(line[2:14])
        temp_ion_alpha.append(line[14:26])
        temp_ion_alpha.append(line[26:38])
        temp_ion_alpha.append(line[38:50])

        ion_alpha = []
        for value in temp_ion_alpha:
            if 'D' in value:
                value = value.replace('D','e')
                print value
            value = float(value)
            ion_alpha.append(value)
        nav_header.append(ion_alpha)
    elif 'ION BETA' in line:
        temp_ion_beta = []
        temp_ion_beta.append(line[2:14])
        temp_ion_beta.append(line[14:26])
        temp_ion_beta.append(line[26:38])
        temp_ion_beta.append(line[38:50])

        ion_beta = []
        for value in temp_ion_beta:
            if 'D' in value:
                value = value.replace('D','e')
                print value
            value = float(value)
            ion_beta.append(value)       
        nav_header.append(ion_beta)
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
        #print hh, mm, ss
        JD = GDtoJD(d,m,y,hh,mm,ss)
        print JD
        Trec = JDtosecofweek(JD)
        dic['satellite_info'][sat_name]['sec_of_week'] = Trec
        print 'Trec: ', Trec
        c = 299792458
        Toe = data['observations'][0]['satellite_info'][sat_name]['params']['TOE']
        print 'TOE: ', Toe
        p2 = float(data['observations'][0]['satellite_info'][sat_name]['P2'])
        print 'P2: ', p2
        #TEMIS should be 135109.927325
        #print p2/c
        Temis1 = (Trec - (p2/c))-Toe
        print 'TEMIS BEFORE: ', Temis1
        if Temis1 > 302400:
            Temis1 = Temis1 - 604800
        elif Temis1 < -302400:
            Temis1 = Temis1 + 604800
        print 'TEMIS1 AFTER', Temis1
        
        a0 = data['observations'][0]['satellite_info'][sat_name]['params']['SV Clock Bias']
        a1 = data['observations'][0]['satellite_info'][sat_name]['params']['SV Clock Drift']
        a2 = data['observations'][0]['satellite_info'][sat_name]['params']['SV Clock Drift Rate']
        Tcorr1 = a0 + a1*Temis1 + a2*(Temis1**2)
        print Tcorr1
        

        Temis2 = (Trec - (p2/c) - Tcorr1) - Toe
        print Temis2
        if Temis2 > 302400:
            Temis2 = Temis2 - 604800
        elif Temis2 < -302400:
            Temis2 = Temis2 + 604800
        print Temis2
        
        Tcorr2 = a0 + a1*Temis2 + a2*(Temis2**2)
        print 'TCORR2', Tcorr2

        Temis = (Trec - (p2/c) - Tcorr2) - Toe
        if Temis > 302400:
            Temis = Temis - 604800
        elif Temis < -302400:
            Temis = Temis + 604800
        print 'Final emission time', Temis
        a = data['observations'][0]['satellite_info'][sat_name]['params']['Sqrt(a)']
        a = a**2

        # compute mean notion
        delta_n = data['observations'][0]['satellite_info'][sat_name]['params']['Delta n']
        n = math.sqrt(3.986005e+14/a**3) + delta_n
        print 'MEAN MOTION: ', n

        # compute mean anomaly
        Mo = data['observations'][0]['satellite_info'][sat_name]['params']['Mo']
        M = Mo + n * Temis
        print 'MEAN ANOMALY: ', M

        # eccentric anomaly
        e = data['observations'][0]['satellite_info'][sat_name]['params']['Eccentricity']
        # first iteration
        E = M + e * math.sin(M)
        print E
        print
        for i in range(10):
            E = M + e * math.sin(E)
            print E
        #calculate true anomaly

        v = math.atan((math.sqrt(1-e**2)*math.sin(E))/(math.cos(E)-e))
        print v
    
        #argument of latitude
        omega = data['observations'][0]['satellite_info'][sat_name]['params']['Omega']
        arglat = v + omega
        print arglat

        #orbital correction terms
        Cus = data['observations'][0]['satellite_info'][sat_name]['params']['Cus']
        Cuc = data['observations'][0]['satellite_info'][sat_name]['params']['Cuc']
        delta_u = Cus * math.sin(2*arglat) + Cuc * math.cos(2*arglat)
        Crs = data['observations'][0]['satellite_info'][sat_name]['params']['Crs']
        Crc = data['observations'][0]['satellite_info'][sat_name]['params']['Crc']
        delta_r = Crs * math.sin(2*arglat) + Crc * math.cos(2*arglat)
        Cis = data['observations'][0]['satellite_info'][sat_name]['params']['CIS']
        Cic = data['observations'][0]['satellite_info'][sat_name]['params']['Cic']
        delta_i = Cis * math.sin(2*arglat) + Cic * math.cos(2*arglat)
        print delta_u, delta_r, delta_i

        #argument of latitude, radius and inclination:
        arglat = arglat + delta_u
        print arglat
        r = a*(1-e*math.cos(E)) + delta_r
        print r
        Io = data['observations'][0]['satellite_info'][sat_name]['params']['Io']
        IDOT = data['observations'][0]['satellite_info'][sat_name]['params']['IDOT']
        i = Io + delta_i + IDOT * Temis
        print i

        #position of orbital plane
        xop = r * math.cos(arglat)
        print xop
        yop = r * math.sin(arglat)
        print yop

        #correct longitude of ascending node
        rotation = 7.2921151467e-5
        OMEGA = data['observations'][0]['satellite_info'][sat_name]['params']['OMEGA']
        OMEGA_DOT = data['observations'][0]['satellite_info'][sat_name]['params']['OMEGA DOT']
        lon_asc = OMEGA + (OMEGA_DOT - rotation) * Temis - rotation * Toe
        print lon_asc

        # final satellite coordinates
        sat_x = xop * math.cos(lon_asc)-yop*math.cos(i)*math.sin(lon_asc)
        sat_y = xop * math.sin(lon_asc)+yop*math.cos(i)*math.cos(lon_asc)
        sat_z = yop*math.sin(i)
        print sat_x
        print sat_y
        print sat_z

        dic['satellite_info'][sat_name]['X'] = float(sat_x)
        dic['satellite_info'][sat_name]['X'] = float(sat_y)
        dic['satellite_info'][sat_name]['X'] = float(sat_z)
        
        #final corrections
        ### SENOR EL DOCUMENTO NO ES CORRECTO ###
        trel = -2 * (math.sqrt(3.986005e+14*a)/c**2) * e * math.sin(E)
        print 'TREL', trel
        # total group delay/satellite instrumental bias
        TGD_L1 = data['observations'][0]['satellite_info'][sat_name]['params']['TGD']
        TGD_L2 = 1.65*TGD_L1
        print TGD_L2
        Tcorr = Tcorr2 + trel - TGD_L2
        print 'Tcorr', Tcorr

        if data['observations'][0]['flag'] == 0:
            
            traveltime = p2/c
            
            R3 = np.array([[math.cos(traveltime*rotation),math.sin(traveltime*rotation),0],
                          [-math.sin(traveltime*rotation),math.cos(traveltime*rotation),0],
                          [0,0,1]])
            sat_x_rot = np.dot(R3,np.array([[sat_x],[sat_y],[sat_z]]))[0]
            print float(sat_x_rot)
            sat_y_rot = np.dot(R3,np.array([[sat_x],[sat_y],[sat_z]]))[1]
            print float(sat_y_rot)
            sat_z_rot = np.dot(R3,np.array([[sat_x],[sat_y],[sat_z]]))[2]
            print float(sat_z_rot)

            ecef = Proj(proj='geocent', ellps='WGS84', datum='WGS84')
            lla = Proj(proj='latlong', ellps='WGS84', datum='WGS84')
            sat_lon, sat_lat, sat_alt = transform(ecef,lla,sat_x_rot,sat_y_rot, sat_z_rot)
            print sat_lon, sat_lat, sat_alt

            Nutm, Eutm, zone_nr, zone_letter,  = utm.from_latlon(54.9878505, 52.22551822)
            print Nutm, Eutm

            ### WE ALSO NEED THE ELEVATION AND AZIMUTH ###
            
            if elevation > 10:
                z = 90 - elevation
                h_ECEF = sat_z_rot #elipsoidal height
                p = 1013.25*(1-0.000065*h_ECEF)**5.225
                T = 291.15-0.0065*h_ECEF
                H = 50*np.exp(-0.0006396*h_ECEF)
                e = (H*0.01)*np.exp(-37.2465+0.213166*T*(0.000256908*T**2))
                
                #find smallest diff with height
                #find corresponding B(mb)
                
                B_dic = {0.0:1.156,
                         0.5:1.079,
                         1.0:1.006,
                         1.5:0.938,
                         2.0:0.874,
                         2.5:0.813,
                         3.0:0.757,
                         4.0:0.654,
                         5.0:0.563}
                B = B_dict[min(B_dic, key=lambda x:abs(x-2.2))]
                    
                dtropo = (0.002277/math.cos(z))*(p+(1255/T+0.05)*e-B*(math.tan(z))**2 )

                #ionospheric delay with klobuchar model

                #1 calcualte earth centered angle
                eca = 0.0137/E+0.11 - 0.022 #semicircles

                #2 compute the latitude of the ionospheric Pierce Point (IPP)
                # u_lat has to be in degrees!!! and corresponds with the users position
                # A is azimuth of the satellite in radians
                ipp_lat = u_lat/180 + eca*math.cos(A)
                if ipp_lat > 0.416:
                    ipp_lat = 0.416
                elif ipp_lat < -0.416:
                    ipp_lat = -0.416

                #3 compute the longitude of ionospheric Pierce Point (IPP)
                # u_lat has to be in degrees!!! and corresponds with the users position
                ipp_lon = u_lon + eca*math.sin(A)/math.cos(ipp_lat)

                #4 find the geomagnetic latitude of the IPP
                geom_lat = ipp_lat + 0.064*math.cos((ipp_lon-1.617)*math.pi)

                #5 find the local time at the IPP
                ipp_time = 43200 * ipp_lon + sec_of_week
                if ipp_time >= 86400:
                    ipp_time -= 86400
                elif ipp_time < 0:
                    ipp_time += 86400

                #6 compute the amplitude of the ionospheric delay in seconds
                Ai = []
                for i in range(4):
                    Ai.append(nav_header[0][i]*geom_lat**i)
                
                Ai = sum(Ai)
                if Ai < 0:
                    Ai = 0

                #7 compute the period of the ionospheric delay in seconds
                Pi = []
                for i in range(4):
                    Pi.append(nav_header[0][i]*geom_lat**i)
                
                Pi = sum(Pi)
                if Pi < 72000:
                    Ai = 72000                    
                #8 compute the phase of the ionospheric delay in radians!!!
                Xi = 2*math.pi*(ipp_time-50400)/Pi

                #9 compute the slant factor (elevation E in semicircles)
                F = 1.0 + 16.0*(0.53-E)**3

                #10 compute the ionospheric delay time
                if Xi < 1.57:
                    IL1 = 5e-9+Ai*(1-(Xi**2/2)+(Xi**4/24))*F
                else:
                    IL1 = 5e-9 *F

                #11 compute the ionospheric time delay for L2 in seconds
                IL2 = 1.65 * IL1

                #12 transform to meters
                IL2 = IL2 * c

                

                

            break
        else:
            #it not a GPS Satellite
            #add 'None' items to dictionary?
            continue
            
    else:
        continue

#print data['observations'][0]['satellite_info']['G02']

#for k in data['observations'][0]['satellite_info']['G02']:
    #print k
