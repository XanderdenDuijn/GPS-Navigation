import math

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

#from xyz lat lon to enu
def xyzlatlon(x,y,z,lat,lon):
    e = -math.sin(lon)*x+math.cos(lon)*y
    n = -math.sin(lat)*math.cos(lon)*x-math.sin(lat)*math.sin(lon)*y+math.cos(lat)*z
    u = math.cos(lat)*math.cos(lon)*x+math.cos(lat)*math,sin(lon)*y+mat.sin(lat)*z
    return e,n,u

def enulatlon(e,n,u,lat,lon):
    x=-math.sin(lon)*e-math.sin(lat)*math.cos(lon)*n+math.cos(lat)*math.cos(lon)*u;
    y=math.cos(lon)*e-math.sin(lat)*math.sin(lon)*n+math.cos(lat)*math.sin(lon)*u;
    z=0+math.cos(lat)*n+math.sin(lat)*u
    return x,y,z
