from math import *

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
###from xyz lat lon to enu
##def xyzlatlon(x,y,z,lat,lon):
##    e = -math.sin(lon)*x+math.cos(lon)*y
##    n = -math.sin(lat)*math.cos(lon)*x-math.sin(lat)*math.sin(lon)*y+math.cos(lat)*z
##    u = math.cos(lat)*math.cos(lon)*x+math.cos(lat)*math,sin(lon)*y+mat.sin(lat)*z
##    return e,n,u
##
##def enulatlon(e,n,u,lat,lon):
##    x=-math.sin(lon)*e-math.sin(lat)*math.cos(lon)*n+math.cos(lat)*math.cos(lon)*u;
##    y=math.cos(lon)*e-math.sin(lat)*math.sin(lon)*n+math.cos(lat)*math.sin(lon)*u;
##    z=0+math.cos(lat)*n+math.sin(lat)*u
##    return x,y,z
