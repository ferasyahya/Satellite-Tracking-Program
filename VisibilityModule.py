import Fileio as stk
import dateAndTimeCalculations as dt
import datetime as datetime
from OMPython import ModelicaSystem
import numpy as np
import math
import LinkCalculation as link

def AOSLOS(AOSLOSfile, startTime, stopTime, TimeStep, F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath,Latm):

    #Extract the year, month, date, hour, minute and second from the start time
    (YR, MO, D, HR, MIN, SEC) = dt.formatDate(startTime)

    #Compute the interval between the two dates to use for the OM simulation
    interval = dt.TimeBetweenTwoDates(startTime, stopTime)

    #Open the TLE file
    f = open('gps-ops.txt','r')

    #Determine the number of lines in the file to use for the while loop for reading and storing all satellites
    num_lines = sum(1 for line in open('gps-ops.txt'))

    #Create the neccessary lists required for storing the number, name, AOS and LOS of each satellite
    name_list = []
    satNum_list = []
    AOS_list = []
    LOS_list = []
    dBm_list = []
    
    #Reading OM files from python
    mod=ModelicaSystem("SatellitePackage.mo","SatellitePackage.test","Modelica.SIunits.Conversions")

    #Tracking hour for GMST
    hr0 = (HR + MIN/60 + SEC/3600)
    
    
    #time since J2000 in days
    MICS = 0
    timeSinceJ2000 = dt.TimeSinceJ2000(YR,MO,D,HR,MIN,SEC,MICS)
    
    #Julian date of current time in days
    JulDay = float(dt.gregToJulian(YR,MO,D,HR,MIN,SEC));            

    for i in range ((int(num_lines/3))):

        #Reading Each Satellite Lines/Storing Number
        line0 = f.readline()
        line1 = f.readline()
        line2 = f.readline()
        sat = stk.ReadNoradTLE(line0, line1, line2)
        name_list.append(sat.name)
        satNum_list.append(i+1)

        #Print out the name of each satellite
        print(name_list[i])
        
        #extracting orbital elements from TLE
        RefEpoch = sat.refepoch
        eccn = float("0." + sat.eccn)
        incl = float(sat.incl)
        raan = float(sat.raan)
        argper = float(sat.argper)
        meanan = float(sat.meanan)
        meanmo = float(sat.meanmo)
        ndot = float(sat.ndot)
        ndot6 = 0
        bstar = 0
        orbitnum = float(sat.orbitnum)
        
        #time since epoch in seconds
        epsec = dt.timeSinceEpoch(RefEpoch,YR,MO,D,HR,MIN,SEC,MICS);    
        
        #Set the parameters in OM
        mod.setParameters(**{"GPS_Test.ecc":eccn,"GPS_Test.M0":meanan,"GPS_Test.N0":meanmo,
                             "GPS_Test.Ndot2":ndot,"GPS_Test.Nddot6":ndot6,"GPS_Test.tstart":epsec,
                             "ARO.elevation":2604.2,"ARO.longitude":281.927,"ARO.latitude":45.9555,
                             "JulianDate":timeSinceJ2000, "hr0":hr0, "raan":raan, "inc":incl, "argper":argper})

        #Obtain parameters from OM
        print(mod.getParameters())
        
        #Simulate
        mod.setSimulationOptions(startTime=0, stopTime=interval, stepSize=TimeStep)
        mod.simulate()
        
        #Get Elevation for satellite index
        (Elevation,time) = mod.getSolutions("Elevation", "time")
        
        #time = time+epsec

        (Px,Py,Pz,Vx,Vy,Vz) = mod.getSolutions("PositionTOPO.x","PositionTOPO.y",
                                               "PositionTOPO.z","VelocityTOPO.x",
                                               "VelocityTOPO.y","VelocityTOPO.z")

        Position = [Px, Py, Pz]
        Velocity = [Vx,Vy,Vz]

        Elevation = Elevation*(180/(math.pi))

        #print(Elevation[0]);

        index = 0;
        AOS = 0;
        LOS = 0;
        
        #Acquire the first acquisition, else 0[no view]
        while index < (len(time)):
            if 9 > Elevation[index] or 89 < Elevation[index]:
                AOS = 0; #no acquisition
                index += 1
            elif 9 < Elevation[index] < 89:
                AOS = 1; #first acquisition
                break;

        if AOS == 1:
            daysSinceTracking = time[index]/86400
            currTimeInJulians = daysSinceTracking + JulDay
            timeFormatted = dt.ep2dat(currTimeInJulians)
            #print(timeFormatted)
            #AOS_list.append(time[index]);
            AOS_list.append(timeFormatted);

            """
            #Doppler Calculations
            dopplerTop = Px[index]*Vx[index] +Py[index]*Vy[index] + Pz[index]*Vz[index]
            dopplerBottom = np.sqrt(Px[index]*Vx[index] +Py[index]*Vy[index] + Pz[index]*Vz[index])
            doppler = dopplerTop/dopplerBottom 
            c = 300000  
            fr = -1227000 
            fD = (doppler/c)*fr
            
            dBm_list.append(fD)
            """
            rangeTOPO = (Px[index]**2 + Py[index]**2 + Pz[index]**2)**(1/2)
            power = link.LinkDesign(F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath,Latm, rangeTOPO)
            dBm_list.append("{0:.2f}".format(power))
        else:
            AOS_list.append('none');
            dBm_list.append('none');
        
        
        #Acquire Loss of Sight, time
        if index == len(time):
            LOS = 0; # no acquision, no loss
        else:
            while index < len(time):
                if 9 < Elevation[index] < 89:
                    LOS = 1; # no loss
                    index = index + 1;
                elif 9 > Elevation[index] or 89 < Elevation[index]:
                    LOS = 2; #first loss
                    break;
                
        if LOS == 0:
            LOS_list.append('none');
        elif LOS == 1:
            LOS_list.append('none');
        else:
            daysSinceTrackingLOS = time[index]/86400
            currTimeInJuliansLOS = daysSinceTrackingLOS + JulDay
            timeFormattedLOS = dt.ep2dat(currTimeInJuliansLOS)
            #LOS_list.append(time[index]);
            LOS_list.append(timeFormattedLOS);
            
        #print(name_list[i])
        print("AOS: " + AOS_list[i]) #acqusition time = Seconds after start time, 'none' = no acquisition
        print("LOS: " + LOS_list[i]) #loss time = Seconds after start time, 'none' = no loss

        #break;

    #NEED to calculate DBM.(???)
    #Not sure about time format, left at seconds after start time
    #Need to extract PRN number and Satellite Number from Each name List

    d = open(AOSLOSfile + '.txt','w')
    d.write('Sat No.\tName\tAOS\t\t\tLOS\t\t\tMin. Expected Level(dBm)\n');
    for j in range (0,len(name_list)):
        if AOS_list[j] == 'none':
            d.write(str(satNum_list[j]) + '\t' + name_list[j] + '\t' + AOS_list[j] + '\t\t\t' +  LOS_list[j] + '\t\t\t'+ str(dBm_list[j])+'\n')
        elif LOS_list[j] == 'none':
            d.write(str(satNum_list[j]) + '\t' + name_list[j] + '\t' + AOS_list[j] + '\t' +  LOS_list[j] + '\t\t\t'+ str(dBm_list[j])+'\n')
        else:
            d.write(str(satNum_list[j]) + '\t' + name_list[j] + '\t' + AOS_list[j] + '\t' +  LOS_list[j] + '\t'+ str(dBm_list[j])+'\n')
    d.close()

    return


