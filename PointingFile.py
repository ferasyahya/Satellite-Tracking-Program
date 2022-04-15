import Fileio as stk
import dateAndTimeCalculations as dt
import datetime as datetime
from OMPython import ModelicaSystem
import math


def PointingFile(PointingFileName, SatToTrack, startTime, stopTime, TimeStep, Azimuth, Elevation):

    #Extract the year, month, date, hour, minute and second from the start time
    (YR, MO, D, HR, MIN, SEC) = dt.formatDate(startTime)
    """
    f = open('gps-ops.txt', 'r')
    num_lines = sum(1 for line in open('gps-ops.txt'))

    #To track a certain satellite, choose the required TLE file, then, input the satellite name in the while loop below.
    #After choosing the satellite to track, choose the tracking time and then set the tracking duration

    #Extract the year, month, date, hour, minute and second from the start time
    (YR, MO, D, HR, MIN, SEC) = dt.formatDate(startTime)
    MICS = 0;
    #Compute the interval between the two dates to use for the OM simulation
    interval = dt.TimeBetweenTwoDates(startTime, stopTime)

    line0=''
    while SatToTrack not in line0:

        line0=f.readline()
        line1=f.readline()
        line2=f.readline()

    f.close()

    sat = stk.ReadNoradTLE(line0,line1, line2)

    RefEpoch = sat.refepoch;
   

    epsec = dt.timeSinceEpoch(RefEpoch,YR,MO,D,HR,MIN,SEC,MICS);    #time since epoch in seconds
    timeSinceJ2000 = dt.TimeSinceJ2000(YR,MO,D,HR,MIN,SEC,MICS);    #time since J2000 in days
    JulDay = float(dt.gregToJulian(YR,MO,D,HR,MIN,SEC));            #Julian date of current time in days

    #print(timeSinceJ2000)

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

    (Azimuth,Elevation) = mod.getSolutions("Azimuth", "Elevation")


    Azimuth = Azimuth*(180/math.pi)
    Elevation = Elevation*(180/math.pi)
    """
    #writing pointing file

    doy = dt.doy(YR,MO,D)

    simTime = datetime.datetime.utcnow()
    simTime = simTime.replace(year = YR, month=MO, day=D,hour=HR,minute=MIN,second=SEC,microsecond=0)

    f = open(PointingFileName + '.txt', 'w')
    f.write("# UTC Date/Time   Azimuth and AZ_Velocity   Elevation and EL_Velocity\n")
    for i in range (0,len(Azimuth)):
        YR = simTime.year
        HR = simTime.hour
        MIN = simTime.minute
        SEC = simTime.second
        D = simTime.day
        MO = simTime.month
        angle = int(Azimuth[i])
        left = (Azimuth[i]-angle)*60
        minute = int(left)
        seconds = (left - minute)*60

        angleEl = int(Elevation[i])
        leftE = (Elevation[i]-angleEl)*60
        minuteE = int(leftE)
        secondsE = (leftE - minuteE)*60
        
        f.write(("{0:4d}.{1:3d}.{2:02d}:{3:02d}:{4:02d} {5:03d} {6:02d} {7:04.1f} {8:1.1f} {9:03d} {10:02d} {11:04.1f} {12:1.1f}\n").format(YR,doy, HR, MIN, SEC, angle, minute, seconds, 0, angleEl, minuteE, secondsE, 0))
        simTime = simTime + datetime.timedelta(seconds=TimeStep)
        doy = dt.doy(YR,MO,D)

    f.close()
        
