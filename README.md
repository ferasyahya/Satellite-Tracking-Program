# Satellite-Tracking-Program
The Satellite Tracking Program was created throughout the year for a course called Space Hardware (ENG 4350)
The purpose of this project is to read a file containing Two Line Elements of different GPS satellites. The program then analyzes the TLE for each satellite in the file and outputs human-readable information about the satellite such as its orbital elements. Once the file is read and analyzed, the program prompts the user to choose a satellite to track. Note that the ability to track a satellite depends on the location of the radio observatory the user is currently located at. This infomation is also supplied to the program in a text file. For the purpose of this project, the radio observatory used was Algoquin Radio Observatory. Following choosing the satellite, the program outputs a series of timings which indicate when the radio observatory will have a clear line of site to the satellite within the specified time range (The time range is a variable that is perviously chosen by the user).

dateAndTimeCalculations.py: A file that performs neccessary time calculations,
Fileio.py: The file that reads and analyzes the gps-ops.txt file and breaks down the TLE of each satellite into a readable form,
gps-ops.txt: The file that contains the TLE for several GPS satellites,
LinkCalculations.py: A file that performs the neccessary calculations for the link budget,
P5.pdf: Project's document and description,
SatellitePackage.mo: OpenModelica files that are responsible for performing the main calculations of converting between the various coordinate systems,
SensorPointingFile.py: Produces a file that contains information on the azimuith and elevation angles required to point towards the chosen satellite,
STNFIL.txt: File containing the required information of the current radio observatory,
TrackingData.py: Utilizes the openModelica files to aid in converting between the coordinate systems,
UserInputParser.py: Main file to run the program,
VisiblityModule.py: Generates the AOSLOS file that contains the information of AOS and LOS of each of the satellites in the gps-ops.txt file.
