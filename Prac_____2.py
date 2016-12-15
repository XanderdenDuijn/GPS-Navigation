#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import itertools
import re
import json
import numpy as np
import pygmaps
import webbrowser
import utm

def calcB(h):
    km=[0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0]
    Bmb=[1.156,1.079,1.006,0.938,0.874,0.813,0.757,0.654,0.563]

    h=h/1000
    if h<=0:
        B=Bmb[0]
        return B
    else:
     for i in range(0,8):
        if km[i]<h and h<km[i+1]:
            dfm=Bmb[i+1]-Bmb[i]
            dfa=km[i+1]-km[i]
            B=Bmb[i]+h*dfm/dfa
            return B
def geocent2geod(X,Y,Z):
    a=6378137.0
    b=6356752.3142
    f=1/298.257223563
    e=math.sqrt((a**2-b**2)/a**2)
    lon=math.atan(Y/X)
    p=math.sqrt(X**2+Y**2)
    lat=math.atan((Z/p)*(1/(1-e**2)))

    for i in range(0,10):
        N=a**2/(math.sqrt(a**2*math.cos(lat)**2+b**2*math.sin(lat)**2))
        h=p/(math.cos(lat))-N
        lat=math.atan((Z/p)*(1/(1-e**2*(N/(N+h)))))

    return lat,lon,h

def gre2JD(year,month,day,hh,mm,ss):

    hour=hh+mm/60.0+ss/3600.0
    if month <=2:
        y=year-1.0
        m=month+12.0
    if month > 2:
        y=year
        m=month
    JD=int(365.25*y)+int(30.6001*(m+1))+day+hour/24.0+1720981.5
    return JD



def JD2gre(JD):

    b=int(JD+0.5)+1537.0
    c=int((b-122.1)/365.25)
    d=int(365.25*c)
    e=int((b-d)/30.6001)
    hour= (JD+0.5-int(JD+0.5))*24
    day=b-d-int(30.6001*e)
    month= e-1-12*int(e/14.0)
    year=c-4715.0-int((7+month)/10.0)
    return hour, day, month, year

def modulo(x,y):
    resto=x%y
    return resto
def JD2GPS(JD):

    week=int((JD-2444244.5)/7.0)

    day_week=modulo(int(JD+0.5+1),7)
    hour, day, month, year=JD2gre(JD)
    secs= float(((hour/24.0)+day_week)*86400.0)
    return secs



coc=open('brdc0590.11n','r')
coche=open('brdc0590.11n','r').read()
cou=0
observaciones=json.loads(open('test.json').read())

mapa=pygmaps.maps(39.46, -0.37, 12)

count=0
#CONTAR EL NUMERO TOTAL DE LINEAS DEL ARCHIVO
with coc as y:
    for line in y:
        cou=cou+1
print cou
with open('brdc0590.11n','r') as f:
    for line in f:
        count=count+1
        if 'ION ALPHA' in line:
            alpha= line[3:50].replace('D','e').split()
        if 'ION BETA' in line:
            beta= line[3:50].replace('D','e').split()
        if 'LEAP SECONDS' in line:
            leap_seconds= line[3:6].replace(' ','')
        if 'END OF HEADER' in line:
            break

fk=open('brdc0590.11n')
bucle=fk.readlines()



#Valores universales

c=299792458
mu=3.986005e14
gammae=7.2921151467e-5
ww=7.2921151467e-5

Aprox_x_o=float(observaciones['Cabecera']['Aporx_position'][0])
Aprox_y_o=float(observaciones['Cabecera']['Aporx_position'][1])
Aprox_z_o=float(observaciones['Cabecera']['Aporx_position'][2])
P2_pos=observaciones['Cabecera']['P2']

archivo_final = open("Archivo_Final.txt","a")
archivo_final.write('Hora;Min;Seg;X;Y;Z;Num_Sat;eX;eY;eZ;Lat;Lon;h'+'\n')

dant=0

cnt=0
for key in observaciones['Observables']:
    for iii in key:
     matriz_ayuda_sat=[]
       # Epoca_obs=(key[i]['Ano']+key[i]['Mes']+key[i]['Dia']+key[i]['Hora']+key[i]['Minuto']+key[i]['Segundo'][:-8]).replace(' ','')
     if key[iii]['OK_Flag']=='0':
        for kk in key[iii]['Datos']:
           kaala=99999999
           if 'G' in kk:
          #   Ep_Ob_Sat=(kk[1:3]+Epoca_obs)
          #   print Ep_Ob_Sat
             licz=count
           #  print key[i]['Hora']+key[i]['Minuto']+key[i]['Segundo']
             while licz<cou:
                for line2 in bucle[licz:licz+1]:
                    hmso=int(key[iii]['Hora'])+int(key[iii]['Minuto'])/60.0
                    hmsn=int(line2[12:14])+int(line2[15:17])/60.0
                    if int(kk[1:3])==int(line2[0:2]) and int(key[iii]['Ano'])==int(line2[3:5]) and int(key[iii]['Mes'])==int(line2[6:8]) and int(key[iii]['Dia'])==int(line2[9:12]) and abs(hmso-hmsn)<1.0 :
                       # print line2[0:20].replace(' ','')
                 #OBSERVACIONES

                        C1=22527689.086
                        C2='nan'
                        L1=-578290.975
                        L2=-410051.963
                        P2= float(key[iii]['Datos'][kk][P2_pos])
                        if key[iii]['Datos'][kk][P2_pos]<>'nan':
                          #  print P2
                            Ano_o=float(key[iii]['Ano'])
                            Mes_o=float(key[iii]['Mes'])
                            Dia_o=float(key[iii]['Dia'])
                            Hora_o=float(key[iii]['Hora'])
                            Minuto_o=float(key[iii]['Minuto'])
                            Segundo_o=float(key[iii]['Segundo'])
    ##                        print P2
    ##                        print Ano_o
    ##                        print Mes_o
    ##                        print Dia_o
    ##                        print Hora_o
    ##                        print Minuto_o
    ##                        print Segundo_o
    ##                        print Aprox_x_o
    ##                        print Aprox_y_o
    ##                        print Aprox_z_o

                            licz0=0
                            for line22 in bucle[licz:licz+7]:
                                if licz0==0:
                                    Ano_n=int(line22[3:5])
                                    Mes_n=int(line22[6:8])
                                    Dia_n=int(line22[9:12])
                                    Hora_n=int(line22[12:14])
                                    Minuto_n=int(line22[15:17])
                                    Segundo_n=line22[19:22]
                                    a0=float(line22[22:41].replace('D','e'))
                                    a1=float(line22[41:60].replace('D','e'))
                                    a2=float(line22[60:79].replace('D','e'))
                                    #print Ano_n,Mes_n,Dia_n,Hora_n,Minuto_n,Segundo_n,a0,a1,a2
                                if licz0==1:
                                    crs=float(line22[22:41].replace('D','e'))
                                    An=float(line22[41:60].replace('D','e'))
                                    M0=float(line22[60:79].replace('D','e'))
                                   # print crs, An,M0
                                if licz0==2:
                                    cuc=float(line22[3:22].replace('D','e'))
                                    e=float(line22[22:41].replace('D','e'))
                                    cus=float(line22[41:60].replace('D','e'))
                                    raizSemMay=float(line22[60:79].replace('D','e'))
                                   # print cuc,e,cus,raizSemMay
                                if licz0==3:
                                    TOE=float(line22[3:22].replace('D','e'))
                                    cic=float(line22[22:41].replace('D','e'))
                                    gamma0=float(line22[41:60].replace('D','e'))
                                    cis=float(line22[60:79].replace('D','e'))
                                   # print TOE,cic,gamma0,cis
                                if licz0==4:
                                    i0=float(line22[3:22].replace('D','e'))
                                    crc=float(line22[22:41].replace('D','e'))
                                    w=float(line22[41:60].replace('D','e'))
                                    gammap=float(line22[60:79].replace('D','e'))
                                  #  print i0,crc,w,gammap
                                if licz0==5:
                                    ip=float(line22[3:22].replace('D','e'))
                                   # print ip
                                if licz0==6:
                                    TGDL1=float(line22[41:60].replace('D','e'))
                                  #  print TGDL1
                                licz0=licz0+1


                            hora_ref= float(Hora_n)-(float(Hora_o)+Minuto_o/60.0+Segundo_o/3600.0)
                            jde= gre2JD(2000+Ano_o,Mes_o,Dia_o,Hora_o,Minuto_o,Segundo_o)
                            sec_of_week=float(JD2GPS(jde))
                            Temis1=(sec_of_week-P2/c)
                            dt_eph0=Temis1-TOE
                            Tcorr1=a0+a1*(dt_eph0)+a2*(dt_eph0**2)
                            dt_eph1=((Temis1-Tcorr1)-TOE)
                            Temis2=(sec_of_week-P2/c-Tcorr1)-TOE
                            Tcorr2=a0+a1*Temis2+a2*(Temis2**2)
                            dt_eph2=((Temis1-Tcorr2)-TOE)
                            a=raizSemMay**2
                            n=math.sqrt(mu/a**3)+An
                            M=M0+n*dt_eph2
                            jj=0
                            E=M
                            while jj <8:
                                E=M+e*math.sin(E)
                                jj=jj+1
                            v=np.arctan2((math.sqrt(1-e**2)*math.sin(E)),(math.cos(E)-e))
                            phi=v+w
                            du=cus*math.sin(2*phi)+cuc*math.cos(2*phi)
                            dr=crs*math.sin(2*phi)+crc*math.cos(2*phi)
                            di=cis*math.sin(2*phi)+cic*math.cos(2*phi)
                            u=phi+du
                            r=a*(1-e*math.cos(E))+dr
                            ieph=i0+di+ip*dt_eph2
                            xop=r*math.cos(u)
                            yop=r*math.sin(u)
                            gamma=gamma0+(gammap-gammae)*dt_eph2-gammae*TOE
                            mk=[]
                            Xsat=xop*math.cos(gamma)-yop*math.cos(ieph)*math.sin(gamma)
                            Ysat=xop*math.sin(gamma)+yop*math.cos(ieph)*math.cos(gamma)
                            Zsat=yop*math.sin(ieph)

                            trel=-2*(math.sqrt(mu*a)/c**2)*e*math.sin(E)
                            TGDL2=TGDL1*1.65
                            Tcorr=Tcorr2+trel-TGDL2

                            if kk<>kaala:
                               kaala=kk
                               mk=[Xsat, Ysat, Zsat,Tcorr,sec_of_week,kk,P2]
                               matriz_ayuda_sat.append(mk)

                         #   print  kk
    ##                        print TGDL1,TGDL2
    ##                        print hora_ref
    ##                        print 'jde','{0:.15f}'.format(jde)
    ##                        print 'sec_of_week',sec_of_week
    ##                        print 'Temis1',Temis1
    ##                        print 'TOE',TOE
    ##                        print 'dt_eph0',dt_eph0
    ##                        print 'Tcorr1',Tcorr1
    ##                        print 'dt_eph1',dt_eph1
    ##                        print 'Tcorr2',Tcorr2
    ##                        print 'dt_eph2', dt_eph2
    ##                        print 'mean motion', '{0:.15f}'.format(n)
    ##                        print 'mean anomaly', M
    ##                        print 'ano_ecc', E
    ##                        print 'true_ano', v
    ##                        print 'phi', phi
    ##                        print 'du', du
    ##                        print 'dr', dr
    ##                        print 'di', di
    ##                        print 'u', u
    ##                        print 'r', r
    ##                        print 'i', ieph
    ##                        print 'xop', xop
    ##                        print 'yop', yop
    ##                        print 'Omega', gamma
    ##                        print 'Xsat', Xsat
    ##                        print 'Ysat', Ysat
    ##                        print 'Zsat', Zsat
    ##                        print 'trel', trel
    ##                        print 'Tcorr', Tcorr

                licz=licz+8
  #  print matriz_ayuda_sat
  #  print iii
    for iki in range(0,6):
        A_ayuda=[]
        A=[]
        P=[]
        K=[]
        PP=[]
        licz_sat=0
        for kil in matriz_ayuda_sat:
            #print kil
            Xsat=kil[0]
            Ysat=kil[1]
            Zsat=kil[2]
            Tcorr=kil[3]
            sec_of_week=kil[4]
            n_sat=kil[5]
            P2=kil[6]
            travel_time=(math.sqrt((Xsat-Aprox_x_o)**2+(Ysat-Aprox_y_o)**2+(Zsat-Aprox_z_o)**2))/c
            omegatau=ww*travel_time
            R3=[[math.cos(omegatau),math.sin(omegatau),0],[-math.sin(omegatau),math.cos(omegatau),0],[0,0,1]]
            calcsat=[[Xsat],[Ysat],[Zsat]]
            calcsat2=np.dot(R3,calcsat)
            Xsat_rot=calcsat2[0]
            Ysat_rot=calcsat2[1]
            Zsat_rot=calcsat2[2]
            lat_sat,lon_sat,h_sat= geocent2geod(Xsat_rot,Ysat_rot,Zsat_rot)
            lat_pap,lon_pap,h_pap= geocent2geod(Aprox_x_o,Aprox_y_o,Aprox_z_o)
            A_ENU=[[-math.sin(lat_pap)*math.cos(lon_pap), -math.sin(lat_pap)*math.sin(lon_pap),math.cos(lat_pap)],[-math.sin(lon_pap),math.cos(lon_pap),0],[math.cos(lat_pap)*math.cos(lon_pap),math.cos(lat_pap)*math.sin(lon_pap),math.sin(lat_pap)]]
          #  print A_ENU
            Corrd_dif=[[float(Xsat_rot-Aprox_x_o)],[float(Ysat_rot-Aprox_y_o)],[float(Zsat_rot-Aprox_z_o)]]
           # print 'Coord dif', Corrd_dif

            ENU=np.dot(A_ENU,Corrd_dif)
          #  print 'Uppin', ENU

            dis_hor=np.sqrt(ENU[0]*ENU[0]+ENU[1]*ENU[1])
      #      print dis_hor
            if dis_hor<=0.000001:
                eleva_sat=float(math.pi/2)
            else:
                eleva_sat=np.arctan2(ENU[2],dis_hor)
            azi_sat=np.arctan2(ENU[1],ENU[0])
      #      print eleva_sat*180/math.pi
            if azi_sat<0:
                azi_sat=azi_sat+2*math.pi
            if eleva_sat*180/math.pi>10:
                Pres=1013.25*(1-0.000065*h_pap)**(5.225)
                Temp=291.15-0.0065*h_pap
                H=50*np.exp(-0.0006396*h_pap)
                e_wp=(H*0.01)*np.exp(-37.2465+0.213166*Temp-(0.000256908*Temp*Temp))
                z=math.pi/2-eleva_sat
             #   print z,Pres,Temp,e_wp,calcB(h_pap)
                d_tropo=(0.002277/math.cos(z))*(Pres+(1255/Temp+0.05)*e_wp-calcB(h_pap)*math.tan(z)**2)
                kak=(Pres+((1255/Temp)+0.05)*e_wp-calcB(h_pap)*math.tan(z)**2)
                ang_tierra=float(0.0137/((eleva_sat*180/math.pi)/180+0.11)-0.022)
                lat_iono=(lat_pap*180/math.pi)/180+ang_tierra*math.cos(azi_sat)
                if lat_iono>0.416:
                    lat_iono=0.416
                if lat_iono<-0.416:
                    lat_iono=-0.416
                lon_iono=(lon_pap*180/math.pi)/180+(ang_tierra*math.sin(azi_sat))/math.cos(lat_iono*math.pi)
                lat_geom=lat_iono+0.064*math.cos((lon_iono-1.617)*math.pi)
                t_local=43200*lon_iono+sec_of_week
                if t_local<0:
                    t_local=t_local+86400
                if t_local>=86400:
                    t_local=t_local-86400
                A_iono=(float(alpha[0]))+(float(alpha[1])*lat_geom)+(float(alpha[2])*lat_geom**2)+(float(alpha[3])*lat_geom**3)
                P_ret=(float(beta[0]))+(float(beta[1])*lat_geom)+(float(beta[2])*lat_geom**2)+(float(beta[3])*lat_geom**3)
                if P_ret<72000:
                    P_ret=72000
                fase_ret_iono=(2*math.pi*(t_local-50400))/P_ret
                #INCLINACION EN SEMICIRCULOS
                Fac_incl=1.0+16.0*(0.53-eleva_sat/math.pi)**3
                if fase_ret_iono<=1.57:
                    Ret_iono=(5e-9+A_iono*(1-(fase_ret_iono**2/2)+(fase_ret_iono**4/24)))*Fac_incl
                if fase_ret_iono>1.57:
                    Ret_iono=5e-9*Fac_incl
                dlon1=Ret_iono*c
                dlon2=dlon1*1.65
                dist=math.sqrt((Xsat_rot-Aprox_x_o)**2+(Ysat_rot-Aprox_y_o)**2+(Zsat_rot-Aprox_z_o)**2)
               # print P2
                K_ayuda=[float(P2-dist-dant+c*Tcorr-d_tropo-dlon2)]
              #  print dlon2
                K.append(K_ayuda)
    ##            dant=d_tropo
                P_ayuda=[float(math.sin(eleva_sat)**2/(1.5*(0.3))**2)]
                PP.append(P_ayuda)
                A_ayuda=[float(-(Xsat_rot-Aprox_x_o)/dist),float(-(Ysat_rot-Aprox_y_o)/dist),float(-(Zsat_rot-Aprox_z_o)/dist),1]
                A.append(A_ayuda)
                licz_sat=licz_sat+1
       #     print 'SAT Y PUNTO', Zsat_rot,Zsat
      #      print eleva_sat
##            K_ayuda=
           # print kil
            #print dis_hor
##            if iki==5 and n_sat=='G02':
##                print 'omegatau', omegatau
##                print 'Xsat_rot', Xsat_rot
##                print 'Ysat_rot', Ysat_rot
##                print 'Zsat_rot', Zsat_rot
##                print 'Lat sat', lat_sat
##                print 'Lon sat', lon_sat
##                print 'h sat', h_sat
##                print 'Lat p', lat_pap
##                print 'Lon p', lon_pap
##                print 'h p', h_pap
##                print 'Elevacion Satelite', eleva_sat*180/math.pi
##                print 'Azimut Satelite', azi_sat*180/math.pi
##                print 'Presion', Pres
##                print 'Temperatura', Temp
##                print 'Humedad', H
##                print 'Humedad relativa', e_wp
##                print 'Correcion Troposfercia',d_tropo
##                print 'Psi', ang_tierra
##                print 'Lat_i', lat_iono
##                print 'Lon_i', lon_iono
##                print 'Lat_geom', lat_geom
##                print 'Tiempo local', t_local
##                print 'A_iono', A_iono
##                print 'P_ret', P_ret
##                print 'Retraso Ionosferico', fase_ret_iono
##                print 'Factor Inclinacion', Fac_incl
##                print 'Delay Longitud L1',dlon1
##                print 'Delay Longitud L2',dlon2

##            print 'omegatau', omegatau
##            print 'Xsat_rot', Xsat_rot
##            print 'Ysat_rot', Ysat_rot
##            print 'Zsat_rot', Zsat_rot
##            print 'Lat sat', lat_sat
##            print 'Lon sat', lon_sat
##            print 'h sat', h_sat
##            print 'Lat p', lat_pap
##            print 'Lon p', lon_pap
##            print 'h p', h_pap
##            print 'Elevacion Satelite', eleva_sat*180/math.pi
##            print 'Azimut Satelite', azi_sat*180/math.pi
##            print 'Presion', Pres
##            print 'Temperatura', Temp
##            print 'Humedad', H
##            print 'Humedad relativa', e_wp
##            print 'Correcion Troposfercia',d_tropo
##            print 'Psi', ang_tierra
##            print 'Lat_i', lat_iono
##            print 'Lon_i', lon_iono
##            print 'Lat_geom', lat_geom
##            print 'Tiempo local', t_local
##            print 'A_iono', A_iono
##            print 'P_ret', P_ret
##            print 'Retraso Ionosferico', fase_ret_iono
##            print 'Factor Inclinacion', Fac_incl
##            print 'Delay Longitud L1',dlon1
##            print 'Delay Longitud L2',dlon2
        #print kil
        A=np.matrix(A)
    #    print 'A', A
##        P=np.diagflat(PP)
##        print P
        K=np.matrix(K)
      #  print 'K', K
        P=np.matrix(PP)
       # print P
        P=np.diagflat(P)
 #       print 'P', P
        k=np.dot(A.transpose(),P)

        s=np.linalg.inv(np.dot(k,A))
        k2=np.dot(A.transpose(),P)
        k3=np.dot(k2,K)
        xfin=np.dot(s,k3)
        dX=xfin[0]
        dY=xfin[1]
        dZ=xfin[2]
        dT=xfin[3]

        res=np.dot(A,xfin)-K
##
        Aprox_x_o=Aprox_x_o+dX
        Aprox_y_o=Aprox_y_o+dY
        Aprox_z_o=Aprox_z_o+dZ
        dant=dT



##        print dX
##        print dY
##        print dZ
##        print dT
##        ##
##        print Aprox_x_o
##        print Aprox_y_o
##        print Aprox_z_o

        iki=iki+1
##    print res
    (m,n)=A.shape
##    print m,n
    desv1=np.dot(res.transpose(),P)
    desv2=np.dot(desv1,res)
    des=desv2/(m-n)
##    print des
##    print s
    desx=math.sqrt(des*s[0,0])
    desy=math.sqrt(des*s[1,1])
    desz=math.sqrt(des*s[2,2])
    print iii
    print desx,desy,desz

    (lat1,lon1,h)=geocent2geod(Aprox_x_o,Aprox_y_o,Aprox_z_o)
##    print lat1*180/math.pi,lon1*180/math.pi,h
##    print utm.from_latlon(lat1*180/math.pi,lon1*180/math.pi)
    mapa.addradpoint(lat1*180/math.pi, lon1*180/math.pi,5,'#808000')
##    print licz_sat
    akasa=float(iii[0:2]),float(iii[2:4]),float(iii[4:6]),float(Aprox_x_o[0,0]),float(Aprox_y_o[0,0]),float(Aprox_z_o[0,0]),licz_sat,desx,desy,desz,lat1*180/math.pi,lon1*180/math.pi,h
##    print aka
    archivo_final.write(str(akasa).replace('(','').replace(')','').replace(',',';').replace(' ','')+'\n')

    cnt=cnt+1
##    if cnt==1:
##        break


##
archivo_final.close()
mapa.draw('mapa2.html')


f=open('mapa2.html','r')

webbrowser.open('file://C:\Users\Konrad\Desktop\Programacion\Angel\mapa2.html', 2, True)