from math import *
import numpy as np

def GDtoJD(d,m,y,hh,mm,ss):
    #representing time as a real number
    hour = hh+mm/60.0+ss/3600.0

    if m <= 2:
        y = y-1
        m = m+12
    JD = int(365.25*y)+int(30.6001*(m+1))+d+hour/24.0+1720981.5
    return JD

def JDtosecofweek(JD):
    week = int((JD-2444244.5)/7)
    day_week = int(JD+0.5+1)%7
    decimal = JD - int(JD)
    hour = (JD + 0.5 - int(JD + 0.5))*24
    #seconds + seconds in days
    secs = ((hour/24)%1+day_week)*86400 
    return secs

def xyztolatlonh(x,y,z):
    a=6378137.0
    b=6356752.3142
    f=1/298.257223563
    e=sqrt((a**2-b**2)/(a**2))
    lon=atan(y/x)
    p=sqrt(x**2+y**2)
    
    lat=atan((z/p)*(1/(1-e**2)))

    #iteration
    for i in range(10):
        n=(a**2)/sqrt((a**2)*(cos(lat)**2)+(b**2)*(sin(lat)**2))
        h=(p/(cos(lat)))-n
        lat=atan((z/p)*(1/(1-e**2*(n/(n+h)))))

    lat = lat * (180/pi)
    lon = lon * (180/pi)
    return lat,lon,h

def Elev_Az(NEh):
    # HD : Distancia horizontal
    # NEh : Coordenadas cartesianas (N, E, h)
    
    # Horizontal distance
    HD = (NEh[0]**2 + NEh[1]**2)**0.5


    if abs(HD/NEh[2]) <= 0.0001:
        El = pi/2
    else:
        El = np.arctan2(NEh[2], HD)

    Az = np.arctan2(NEh[1], NEh[0])
    if Az<0:
        Az += 2*pi

    return El, Az

