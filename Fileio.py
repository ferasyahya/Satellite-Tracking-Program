#this module defines all functions needed to display outputs and read in inputs
import winsound
from collections import namedtuple
from OMPython import ModelicaSystem
import numpy as np
import os

os.chdir("D:\My Files\School Work\Winter 2018\ENG 4350 WINTER\P5")

#banner with group info
def banner():
    print ("Group 2\n Rajika Pati Arambage\n Feras Yahya\n Mark Lopez\n March 6 2018\n Welcome!")
    if (open("STNFIL.txt", "r")):
        print("Here")


#error message with beep
def errmsg(ERRMSG):
    print("%s" % (ERRMSG))
    winsound.Beep(2000,700)  #frequency of 2000Hz and duration of 700ms
    
#Reading station file
def ReadStationFile():
    Station = ""
    while Station == "":
        STNFILE = input("Input file name to open: ")
        try:
            with open(STNFILE, "r") as file:
                data = file.readlines()
            file.close()

            #array of required elements for structure    
            elements = ["name","stnlat","stnlong","stnalt","utc_offset","az_el_nlim", "az_el_lim","st_az_speed_max","st_el_speed_max"]
            #array of required elements for angles structure
            az = ["az", "elmin", "elmax"]

            Elements_tuple = namedtuple ('station_info', elements)
            Angles_tuple = namedtuple ('angle_limits', az)

            #a_list = []

            #for i in range (6, 6)
            lims = data[6].split(",")
            lims = [float(j) for j in lims]
            limits = (Angles_tuple(lims[0], lims[1],lims[2]))
            #creating final structure (tuple) called station 
            Station = Elements_tuple(str(data[0]),float(data[1]),float(data[2]),float(data[3]) \
                                        ,float(data[4]), float(data[5]), limits,float(data[7]),float(data[8]))


        except:
            print(f"File does not exist: {STNFILE}")
            Station =""

    return Station

#Reading NORAD TLE file
def ReadNoradTLE (line0, line1, line2):
    
    #array of required elements for structure
    tle_elements = ["name","refepoch","incl","raan","eccn","argper", \
                 "meanan","meanmo","ndot", "nddot6", "bstar", "orbitnum"]

    TLE_tuple = namedtuple('TLE_info', tle_elements)

    Satellite = TLE_tuple(str(line0[13:19]), line1[18:32], line2[8:16], line2[17:25], line2[26:33], line2[34:42], \
                          line2[43:51], line2[52:63], line1[33:43], line1[44:52], line1[53:61], line2[63:68])
    
    #str(line0[0:70])
    return Satellite



#TESTING ReadStationFile
#print(ReadStationFile('STNFIL.txt'))


"""
#TESTING ReadNoradTLE
f = open('gps-ops.txt', 'r')
num_lines = sum(1 for line in open('gps-ops.txt'))
#print(num_lines)

satellite_list = [];

for i in range ((int(num_lines/3))):
    line0 = f.readline();
    line1 = f.readline();
    line2 = f.readline(); 
    Satellite = ReadNoradTLE(line0, line1, line2)
    satellite_list.append(Satellite)

#print(satellite_list[29])
"""

"""
The function generates an Ephermis file with a format similar to that generated using STK

Arguments:
1. outfile(string) File writing to (includes file ty: name.txt)
2. EphemFile (string)
3. StartString. Format: 24 Apr 2005 00:00:00.00 (string) Tracking time
4. time    (double array) Time results generated from OM
5. Coord   (string) Coordinate in which calculation are done for
6. position  (3D double array) holds x,y,z position values generated from OM
7. velocity  (3D double array) holds x,y,z velocity values generated from OM

Returns:
- Generates an Ephemeris File with the specified name (outfile)
    
Identification: 
    author: Feras Yahya
"""
def STKout(outfile,EphemFile,StartString,time,
Coord,position,velocity):

   
    #path = "C:\\Users\\Owner\\Desktop\\" + outfile 
    fid = open(outfile, 'w')
    fid.write('stk.v.4.3\n\n')
    fid.write('BEGIN Ephemeris\n\n')
   
    
    fid.write('NumberOfEphemerisPoints {0:4d}\n\n'.format(np.size(time)))
    fid.write('ScenarioEpoch \t\t {}\n'.format(StartString))
    fid.write('InterpolationMethod \t Lagrange\n')
    fid.write('InterpolationOrder \t 7\n')
    fid.write('CentralBody \t\t Earth\n')
    fid.write('CoordinateSystem \t {}\n\n'.format(Coord))

    fid.write('EphemerisTimePosVel\n\n')

    PositioninTopox = position[0]*1000
    PositioninTopoy = position[1]*1000
    PositioninTopoz = position[2]*1000

    VelocityinTopox = velocity[0]*1000
    VelocityinTopoy = velocity[1]*1000
    VelocityinTopoz = velocity[2]*1000
  
    for i in range(np.size(time)): 
        fid.write('{0:-16.14E} {1:-16.14E} {2:-16.14E} {3:-16.14E} {4:-16.14E} {5:-16.14E} {6:-16.14E}\n'.format(time[i],PositioninTopox[i],PositioninTopoy[i],PositioninTopoz[i],VelocityinTopox[i],VelocityinTopoy[i],VelocityinTopoz[i]))

    fid.write('\n\nEND Ephemeris')

    fid.close()

