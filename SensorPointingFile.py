import dateAndTimeCalculations as dt
import PointingFile as point
#import CalculatePositionAndVelocity as posvel
import TrackingData as posvel
import Fileio as io
import numpy as np

def SensorPointingFile(FileName, time, Azimuth, Elevation):
    

    f = open(FileName + '.txt', 'w')
    f.write("stk.v.4.1.1\n")
    f.write("Begin   Attitude\n")
    f.write("NumberofAttitudePoints  {0:4d}\n".format(np.size(time)))
    f.write("Squence\t\t323\n")
    f.write("AttitudeTimeAzElAngles\n")

    for i in range(np.size(time)):
        f.write("{0:.2f}\t{1:.2f}  {2:.1f}\n".format(time[i],Azimuth[i],Elevation[i]))
                
    f.write("End   Attitude")         
    f.close()

    
