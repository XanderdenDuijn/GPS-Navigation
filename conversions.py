import math

def GDtoJD(d,m,y,hh,mm,ss):
    #representing time as a real number
    hour = hh+mm/60.0+ss/3600.0

    if m <= 2:
        y = y-1
        m = m+12
    JD = int(365.25*y)+int(30.6001*(m+1))+d+hour/24.0+1720981.5
    return JD

def JDtoGPS(JD):
    week = int((JD-2444244.5)/7)
    day_week = int(JD+0.5+1)%7
    decimal = JD - int(JD)
    hour = (JD + 0.5 - int(JD + 0.5))*24
    #seconds + seconds in days
    secs = ((hour/24)%1+day_week)*86400 
    return secs
