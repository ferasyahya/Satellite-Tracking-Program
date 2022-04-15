import dateAndTimeCalculations as dt
import VisibilityModule as stk
import PointingFile as point
import TrackingData as posvel
import Fileio as io
import SensorPointingFile as spf
import os

os.chdir("D:\My Files\School Work\Winter 2018\ENG 4350 WINTER\P5")


io.banner()

print(io.ReadStationFile())

startTime = input("\nInput the start time of tracking (Format: YYYY-MM-DD HR:MM:SS)\n")
startTime = str(startTime)
startTime = dt.validate(startTime)

stopTime = input("Input the stop time of tracking (Format: YYYY-MM-DD HR:MM:SS)\n")
stopTime = str(stopTime)
stopTime = dt.validateStopTime(stopTime,startTime)

timeStep = -1
while timeStep < 0:
    timeStep = input("Input time step:")
    timeStep = float(timeStep)

print("\nEnter the following RF characteristics of the GPS signal\n")
F_cnt = input("Enter the Frequency band center (MHz): ")
F_cnt = float(F_cnt)
AE = input("\nEnter the Antenna efficiency: ")
AE = float(AE)
Diameter = input("\nEnter the Antenna diameter (m): ")
Diameter = float(Diameter)
Bandwidth = input("\nEnter the bandwidth (MHz): ")
Bandwidth = float(Bandwidth)
RCV = input("\nEnter the receiving antenna gain RCV (dB): ")
RCV = float(RCV)
RNT = input("\nEnter the receiving antenna noise temperature RNT (deg K): ")
RNT = float(RNT)
Pt = input("\nEnter the power transmitted (W): ")
Pt = float(Pt)
Gt = input("\nEnter the gain of the transmitting antenna (dBi): ")
Gt = float(Gt)
Ltranspath = input("\nEnter the RF Losses in trasmitter path (dB): ")
Ltranspath = float(Ltranspath)
Latm = input("\nEnter any Atmospheric and polarization losses (dB): ")
Latm = float(Latm)

#power = link.LinkDesign(1575.42, 0.5, 46, 2, 56, 200, 25,5, 1.25, 0.5, rangeTOPO)

#generate AOSLOS file
stk.AOSLOS('AOSLOS', startTime,stopTime, timeStep, F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath,Latm)

#Prompt user to enter the satellite to track
satToTrack = input("Choose your satellite to Track from the list above(Format: PRN XX): ")

satToTrack = str(satToTrack)

ephemCoord = input("Choose your coordinate system for ephemeris file\n1) Perifocal\n2) ECI\n3) ECF\n4)Topocentric\n ")
invalidInput = True
#Calculate the position and velocity to be used for ephem file

while invalidInput != False:
    if ephemCoord == 'Perifocal':
        (RefEpoch, time, Position, Velocity,Azimuth, Elevation) = posvel.ComputePositionAndVelocity(satToTrack, 'GPS_Test.p_sat_pf.x','GPS_Test.p_sat_pf.y','GPS_Test.p_sat_pf.z','GPS_Test.v_sat_pf.x','GPS_Test.v_sat_pf.y','GPS_Test.v_sat_pf.z',startTime,stopTime,timeStep,F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath, Latm)
        #generate ephem file
        io.STKout('ephemPERInew.e', 'empty', dt.FormatEpoch(RefEpoch),time,"Custom Perifocal_2 CentralBody/Earth",Position, Velocity)
        invalidInput = False
    elif ephemCoord == 'ECI':
        (RefEpoch, time, Position, Velocity,Azimuth, Elevation) = posvel.ComputePositionAndVelocity(satToTrack, 'Position.x','Position.y','Position.z','Velocity.x','Velocity.y','Velocity.z',startTime,stopTime,timeStep,F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath, Latm)
        io.STKout('ephemECInew.e', 'empty', dt.FormatEpoch(RefEpoch),time,"J2000",Position, Velocity)
        invalidInput = False
    elif ephemCoord == 'ECF':
        (RefEpoch, time, Position, Velocity,Azimuth, Elevation) = posvel.ComputePositionAndVelocity(satToTrack, 'PositionECF.x','PositionECF.y','PositionECF.z','VelocityECF.x','VelocityECF.y','VelocityECF.z',startTime,stopTime,timeStep,F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath, Latm)
        #generate ephem file
        io.STKout('ephemECFnew.e', 'empty', dt.FormatEpoch(RefEpoch),time,"Fixed",Position, Velocity)
        invalidInput = False
    elif ephemCoord == 'Topocentric':
        (RefEpoch, time, Position, Velocity,Azimuth, Elevation) = posvel.ComputePositionAndVelocity(satToTrack, 'PositionTOPO.x','PositionTOPO.y','PositionTOPO.z','VelocityTOPO.x','VelocityTOPO.y','VelocityTOPO.z',startTime,stopTime,timeStep,F_cnt, AE, Diameter, Bandwidth, RCV, RNT, Pt, Gt,Ltranspath, Latm)
        #generate ephem file
        io.STKout('ephemTOPOnew.e', 'empty', dt.FormatEpoch(RefEpoch),time,"Custom TOPOCENTRIC_1 Facility/Algonquin",Position, Velocity)
        invalidInput = False
    else:
        ephemCoord = input("Invalid coordinate system\n1) Perifocal\n2) ECI\n3) ECF\n4)Topocentric\n ")
        invalidInput = True


#generate pointing file
point.PointingFile('point', satToTrack, startTime,stopTime,timeStep,Azimuth, Elevation)

spf.SensorPointingFile('SensorPointing',time,Azimuth, Elevation)

