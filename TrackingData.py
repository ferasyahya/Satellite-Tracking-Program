import Fileio as stk
import dateAndTimeCalculations as dt
import datetime as datetime
from OMPython import ModelicaSystem
import LinkCalculation as LC
import math

def ComputePositionAndVelocity(SatName, Px, Py, Pz, Vx, Vy, Vz, startTime, stopTime, TimeStep,F_cnt, AE, D, B, RG, RNT, Pt, Gt,Ltranspath, Latm):

    #Extract the year, month, date, hour, minute and second from the start time
    (YR, MO, D, HR, MIN, SEC) = dt.formatDate(startTime)

    #Compute the interval between the two dates to use for the OM simulation
    interval = dt.TimeBetweenTwoDates(startTime, stopTime)
    f = open('gps-ops.txt', 'r')
    num_lines = sum(1 for line in open('gps-ops.txt'))

    #To track a certain satellite, choose the required TLE file, then, input the satellite name in the while loop below.
    #After choosing the satellite to track, choose the tracking time and then set the tracking duration

    line0=''
    while SatName not in line0:

        line0=f.readline()
        line1=f.readline()
        line2=f.readline()

    f.close()

    sat = stk.ReadNoradTLE(line0,line1, line2)

    RefEpoch = sat.refepoch;
    MICS = 0
    epsec = dt.timeSinceEpoch(RefEpoch,YR,MO,D,HR,MIN,SEC,MICS);    #time since epoch in seconds
    timeSinceJ2000 = dt.TimeSinceJ2000(YR,MO,D,HR,MIN,SEC,MICS);    #time since J2000 in days
    JulDay = float(dt.gregToJulian(YR,MO,D,HR,MIN,SEC));            #Julian date of current time in days

    #Tracking hour for GMST
    hr0 = (HR + MIN/60 + SEC/3600);
    #extracting orbital elements from TLE
    eccn = float("0." + sat.eccn);
    incl = float(sat.incl);
    raan = float(sat.raan);
    argper = float(sat.argper);
    meanan = float(sat.meanan);
    meanmo = float(sat.meanmo);
    ndot = float(sat.ndot);
    ndot6 = 0;
    bstar = 0;
    orbitnum = float(sat.orbitnum);
    #Reading OM files from python
    mod=ModelicaSystem("SatellitePackage.mo","SatellitePackage.test","Modelica.SIunits.Conversions")

    #Obtain parameters from OM
    print(mod.getParameters())

    #Set the parameters in OM
    mod.setParameters(**{"GPS_Test.ecc":eccn,"GPS_Test.M0":meanan,"GPS_Test.N0":meanmo,
                         "GPS_Test.Ndot2":ndot,"GPS_Test.Nddot6":ndot6,"GPS_Test.tstart":epsec,
                         "ARO.elevation":2604.2,"ARO.longitude":281.927,"ARO.latitude":45.9555,
                         "JulianDate":timeSinceJ2000, "hr0":hr0, "raan":raan, "inc":incl, "argper":argper})

    #Print new parameters
    print(mod.getParameters())

    #Simulate
    mod.setSimulationOptions(startTime=0, stopTime=interval, stepSize=TimeStep)
    mod.simulate()

    #Used for the STK out file
    (time,posx,posy,posz,velx,vely,velz) = mod.getSolutions("time",Px,Py,Pz,Vx,Vy,Vz)

    #Used for making the tracking file
    (posx1,posy1,posz1,velx1,vely1,velz1) = mod.getSolutions('PositionTOPO.x','PositionTOPO.y','PositionTOPO.z','VelocityTOPO.x','VelocityTOPO.y','VelocityTOPO.z')

    (Azimuth,Elevation) = mod.getSolutions("Azimuth","Elevation")
    
    Azimuth = Azimuth*(180/math.pi)
    Elevation = Elevation*(180/math.pi)

    #Used for ephemeris file
    time = time+epsec
    Position = [posx, posy, posz]
    Velocity = [velx,vely,velz]
    
    doy = dt.doy(YR,MO,D)
    simTime = datetime.datetime.utcnow()
    simTime = simTime.replace(year = YR, month=MO, day=D,hour=HR,minute=MIN,second=SEC,microsecond=0)

    f = open('Tracking.txt', 'w')
    f.write("---------------------------------------------------------------------------------------------\n")
    f.write("# UTC\t\t\tAz\tEl\tAz-vel   El-vel   Range   Range Rate\tDoppler\tLevel\n")
    f.write("# UTC\t\t\tDeg\tdeg\tdeg/sec  deg/sec   km\t   km/sec\tkHz\tdBm\n")
    f.write("---------------------------------------------------------------------------------------------\n")

    for i in range (0,len(Azimuth)):
        YR = simTime.year
        HR = simTime.hour
        MIN = simTime.minute
        SEC = simTime.second
        D = simTime.day
        MO = simTime.month
        angle = (Azimuth[i])
        angleEl = (Elevation[i])
        Range = (posx1[i]**2 + posy1[i]**2 + posz1[i]**2)**(1/2)
        RangeRate = (velx1[i]**2 + vely1[i]**2 + velz1[i]**2)**(1/2)
        minPower = LC.LinkDesign(F_cnt, AE, D, B, RG, RNT, Pt, Gt,Ltranspath, Latm, Range) - 30     #dBm
        doppler = (RangeRate/299792.458)*(1575.42*1000)
        
        f.write(("{0:4d}.{1:03d}-{2:02d}:{3:02d}:{4:02d}\t{5:.2f}\t{6:.2f}\t{7:1.1f}\t {8:1.1f}\t{9:.2f}   {10:.2f}\t\t{11:.2f}\t{12:.2f}\n").format(YR,doy, HR, MIN, SEC, angle, angleEl,0,0, Range,RangeRate,doppler,minPower))
        #f.write(("{0:4d}.{1:03d}-{2:02d}:{3:02d}:{4:02d}\t{5:.2f}  {6:.2f}  {7:.2f}  {8:.2f}  {9:.2f}  {10:.2f}\n").format(YR,doy, HR, MIN, SEC,posx[i],posy[i],posz[i],velx[i],vely[i],velz[i]))
        simTime = simTime + datetime.timedelta(seconds=TimeStep)
        doy = dt.doy(YR,MO,D)
    
    f.close()

    return RefEpoch, time, Position, Velocity, Azimuth, Elevation










    
