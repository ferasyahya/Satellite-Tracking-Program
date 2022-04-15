"""
This module library contains a collection of conversion methods used to change the time format and its units.

Created on March 9, 2018
authors: Feras Yahya, Rajika Pati Arambage, Mark Lopez
        
"""
from datetime import datetime
from math import floor
import Fileio as io

"""
The function returns the fraction of the given value. 

Arguments:
1. value - Decimal number
Returns:
- The decimal portion of the value
    
Identification: 
    author: Feras Yahya
"""
def frac(value):
    return value-int(value)

"""
The function determines whether a year is a leap year or not.

Arguments:
1. YR - Integer value of the year
Returns:
- True if YR is leap, false otherwise
    
Identification: 
    author: Feras Yahya
"""
def isLeapYear(YR):
    isLeap = False
    if (YR % 4) == 0 and (YR % 100) !=0:
        isLeap = True
    elif (YR % 400) == 0:
        isLeap= True
    else:
        isLeap= False
    return isLeap

"""
The function adjusts the number of days in a month for a leap year and for a non-leap year.

Arguments:
1. YR - Integer value of the year
Returns:
- An array of values where N-1 correspondes to the month and the value at that index correspondes
number of days (N is the index number).
    
Identification: 
    author: Feras Yahya
"""
def daysInMonth(YR):
    months = [31,28,31,30,31,30,31,31,30,31,30,31]
    mList = list(range(1,13)) #List of months 
    #Checking for a leap year
    if isLeapYear(YR) == True:
        months[1] = 29
    return months

"""
The function determines the number of days for a given date.

Arguments:
1. YR - Year (int) 
2. MO - Month (int)
3. D - Day (int)
Returns:
- Days of the year
    
Identification: 
    author: Feras Yahya
"""
def doy(YR,MO,D):
    if type(YR) == type(''):
        YR = int(YR);MO = int(MO);D = int(D)
    day =0
    #Adding up the number of days 
    months= daysInMonth(YR)
    for mm in range(len(months)):
        if mm == MO-1:
            break
        day = day + months[mm]
    day = day+D
    
    return day

"""
The function calculates the fraction of day at the specified input time.

Arguments:
1. HR - Hour (int)
2. MI - Minutes (int)
3. SEC - Seconds (float)
Returns:
- Hours (float)
    
Identification: 
    author: Feras Yahya
"""
def frcofd(HR,MI,SEC):
    hours = HR+(MI/60)+(SEC/3600)
    return hours

"""
The function converts the epoch date, given in Julians, to standard Gregorian date.

Arguments:
1. JulianDate (int)
Returns:
- Gregorian Date (String)
    
Identification: 
    author: Feras Yahya
"""
def ep2dat(JulainDate):
    I = floor(JulainDate+0.5)
    Fr = abs( I - ( JulainDate + 0.5) )

    if I >= 2299160: 
         A = floor( ( I- 1867216.25 ) / 36524.25 )
         a4 = floor( A / 4 )
         B = I + 1 + A - a4
    else:
         B = I
     

    C = B + 1524
    D = floor( ( C - 122.1 ) / 365.25 )
    E = floor( 365.25 * D )
    G = floor( ( C - E ) / 30.6001 )
    day = floor( C - E + Fr - floor( 30.6001 * G ) )

    if G <= 13.5: 
        month = G - 1
    else:
        month = G - 13
    

    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    

    hour = floor( Fr * 24 )
    minu = floor( abs( hour -( Fr * 24 ) ) * 60 )
    minufrac = ( abs( hour - ( Fr * 24 ) ) * 60 )
    sec = ( abs( minu - minufrac ) * 60)
    AA = ( JulainDate + 1.5 ) / 7
    nd = floor( (abs( floor(AA) - AA ) ) * 7 )
    str1 = "%d-%02d-%02d %02d:%02d:%02d" % (year,month,day,hour,minu,sec)
    return str1

"""
The function returns the current day and date at the time the function is called.


Returns:
- UTC as a string
    
Identification: 
    author: Feras Yahya
""" 
def curday():    
    return str(datetime.utcnow())

"""
The function returns the time (in seconds) since epoch given a specific tracking time

Arguments:
1. Epoch(string)
2. YR   (int) Year of tracking time
3. MO   (int) Month of tracking time
4. D    (int) Day of tracking time
5. HR   (int) Hour of tracking time
6. MIN  (int) Minute of tracking time
7. SEC  (int) Seconds of tracking time
8. MICS (int) Microseconds of tracking time

Returns:
- String consisting of epoch date, tracking time date and time since epoch in seconds
    
Identification: 
    author: Feras Yahya
"""    
def timeSinceEpoch(JulianDate,YR,MO,D,HR,MIN,SEC,MICS):
    year = '2'+'0'+JulianDate[0]+JulianDate[1]
    year = int(year)
    month = float(JulianDate[0:13])
    
    
    date = datetime.strptime(JulianDate[:5],'%y%j')
    dfrac = month - int(month)
    hr = int(dfrac*24)
    rem = dfrac*24 - hr
    mins = int(60*rem)
    rem2 = rem*60 - mins
    secs = int(60*rem2)
    mics = int((rem2*60-secs)*10**6)
    date = date.replace(hour=hr, minute=mins, 
                                second=secs, microsecond=mics)

    trackingtime = datetime.utcnow()
    trackingtime = trackingtime.replace(year=YR,month=MO,day=D,hour=HR, minute=MIN, 
                                second=SEC, microsecond=MICS)
    deltat = trackingtime - date
    epsec = deltat.total_seconds()
    #str1 = "Epoch date:\t" + str(date) + "\n Tracking time:\t" + str(trackingtime) + "\n Time since Epoch (s):\t%f" % (epsec)
    return epsec

"""
The function returns the Julian date given the gregorian date

Arguments:
1. YR   (int) Year of gregorian time
2. MO   (int) Month of gregorian time
3. D    (int) Day of gregorian time
4. HR   (int) Hour of gregorian time
5. MIN  (int) Minute of gregorian time
6. SEC  (int) Seconds of gregorian time

Returns:
- Julain date
    
Identification: 
    author: Feras Yahya
"""  
def gregToJulian(YR,MO,D,HR,MIN,SEC):
    timeUT = HR + (MIN/60) + (SEC/3600)
    JD = ( 367 * YR ) - floor ( 7 * ( YR + floor( ( MO + 9 ) / 12 ) ) / 4 ) - floor( 3 * ( floor( ( YR + ( MO - 9 ) / 7 ) / 100 ) + 1 ) / 4 ) + floor( ( 275 * MO ) / 9 ) + D + 1721028.5 + ( timeUT / 24 );

    return JD

"""
The function returns the time (in seconds) since J2000 given a specific tracking time

Arguments:
1. YR   (int) Year of tracking time
2. MO   (int) Month of tracking time
3. D    (int) Day of tracking time
4. HR   (int) Hour of tracking time
5. MIN  (int) Minute of tracking time
6. SEC  (int) Seconds of tracking time
7. MICS (int) Microseconds of tracking time

Returns:
- Time since J2000 in days
    
Identification: 
    author: Feras Yahya
"""  
def TimeSinceJ2000(YR,MO,D,HR,MIN,SEC,MICS):
    #date = datetime.utcnow()
    #date = date.replace(year=YR,month=MO,day=D,hour=HR, minute=MIN, 
    #                           second=SEC, microsecond=MICS)
    
    #J2000 = datetime.utcnow()
    #J2000 = J2000.replace(year=2000,month=1,day=1,hour=12, minute=0, 
    #                            second=0, microsecond=0)
    #deltat = date - J2000
    #timeSinceJ2000 = deltat.total_seconds()/3600

    J2000 = 2451545.0
    Now = gregToJulian(YR,MO,D,HR,MIN,SEC)

    timeSinceJ2000 = Now-J2000
    
    
    return timeSinceJ2000


"""
The function returns the year, month, date, hours, minutes and seconds of a date of the following format
YYYY-MM-DD HH:MM:SS

Arguments:
1. Date

Returns:
- YR
- MO
- D
- HR
- MIN
- SEC
    
Identification: 
    author: Feras Yahya
"""  
def formatDate(date):
    YR = int(date[0:4])
    MO = int(date[5:7])
    D  = int(date[8:10])
    HR = int(date[11:13])
    MIN= int(date[14:16])
    SEC= int(date[17:19])


    return YR, MO, D, HR, MIN, SEC

"""
The function returns the interval in seconds between two dates of the following format
YYYY-MM-DD HH:MM:SS

Arguments:
1. Date1
2. Date2

Returns:
- Interval (in seconds)

Identification: 
    author: Feras Yahya
""" 
def TimeBetweenTwoDates(date1, date2):

    (YR1, MO1, D1, HR1, MIN1, SEC1) = formatDate(date1)
    (YR2, MO2, D2, HR2, MIN2, SEC2) = formatDate(date2)

    #Jul1 = gregToJulian(YR1,MO1,D1,HR1,MIN1,SEC1);
    #Jul2 = gregToJulian(YR2,MO2,D2,HR2,MIN2,SEC2);

    #interval = Jul2 - Jul1      #time difference in days
    #interval = interval*86400   #time difference in seconds

    
    Dt1 = datetime.utcnow()
    Dt1 = Dt1.replace(year=YR1,month=MO1,day=D1,hour=HR1, minute=MIN1, second=SEC1, microsecond = 0)
    
    Dt2 = datetime.utcnow()
    Dt2 = Dt2.replace(year=YR2,month=MO2,day=D2,hour=HR2, minute=MIN2, second=SEC2,microsecond = 0)
    deltat = Dt2 - Dt1
    interval = int(deltat.total_seconds())

    return interval

"""
The function formats the epoch in TLE as DD MON YYYY HR:MM:SS.MIC for ephemeris file

Arguments:
1. TLE Epoch

Returns:
- String date 

Identification: 
    author: Feras Yahya
""" 
def FormatEpoch(epochTLE):
    year = '2'+'0'+epochTLE[0]+epochTLE[1];
    year = int(year);
    month = float(epochTLE[0:13]);
    
    
    date = datetime.strptime(epochTLE[:5],'%y%j');
    dfrac = month - int(month);
    hr = int(dfrac*24)
    rem = dfrac*24 - hr
    mins = int(60*rem)
    rem2 = rem*60 - mins
    secs = int(60*rem2)
    mics = int((rem2*60-secs)*10**6)
    date = date.replace(hour=hr, minute=mins, 
                                second=secs, microsecond=mics)

    if int(date.month)==1:
        monthString = 'Jan'
    elif int(date.month)==2:
        monthString = 'Feb'
    elif int(date.month)==3:
        monthString = 'Mar'
    elif int(date.month)==4:
        monthString = 'Apr'
    elif int(date.month)==5:
        monthString = 'May'
    elif int(date.month)==6:
        monthString = 'Jun'
    elif int(date.month)==7:
        monthString = 'Jul'
    elif int(date.month)==8:
        monthString = 'Aug'
    elif int(date.month)==9:
        monthString = 'Sep'
    elif int(date.month)==10:
        monthString = 'Oct'
    elif int(date.month)==11:
        monthString = 'Nov'
    elif int(date.month)==12:
        monthString = 'Dec'

    stringEpoch = str(date.day) + ' ' + monthString + ' ' + str(year) + ' {0:02d}'.format(date.hour) + ':{0:02d}'.format(date.minute) + ':{0:02d}'.format(date.second) + '.' + str(date.microsecond)
    return stringEpoch

"""
The function determines whether the date entered by the user is in a valid format or not. The function
also determines whether the year, month, date, hour, minute or seconds entered are valid. The requirement
for each argument is as follows:
Year between 2010 and 2020
Month between 1 and 12
Day between 1 and 31
Hour between 0 and 23
Minute between 0 and 59
Second between 0 and 59

Arguments:
1. Date

Returns:
- Validated date
    
Identification: 
    author: Feras Yahya
"""  

def validate(date):

    buffer = False
    tries = 0
    notValidated = True
    while buffer != True:
        if tries > 0:
            date = input("TRY AGAIN (Format: YYYY-MM-DD HR:MM:SS)\n")
            notValidated = True
            
        if notValidated != False:
            try:
                YR = int(date[0:4])
                MO = int(date[5:7])
                D  = int(date[8:10])
                HR = int(date[11:13])
                MIN= int(date[14:16])
                SEC= int(date[17:19])
                notValidated = False
            except:
                io.errmsg("Invalid format!")
                tries = tries+1
                buffer = False
                notValidated = True
                
        if notValidated != True:
            if YR > 2020 or YR < 2010:
                io.errmsg("Enter a valid year please! Year must be between 2010 and 2020")
                tries = tries+1
                buffer = False
            elif MO > 12 or MO < 1:
                io.errmsg("Enter a valid month please! Month must be between 1 and 12")
                tries = tries+1
                buffer = False
            elif D > 31 or D < 1:
                io.errmsg("Enter a valid day please! Day must be between 1 and 31")
                tries = tries+1
                buffer = False
            elif HR > 23 or HR < 0:
                io.errmsg("Enter a valid hour please! Hour must be between 0 and 23")
                tries = tries+1
                buffer = False
            elif MIN > 59 or MIN < 0:
                io.errmsg("Enter a valid minute please! Minute must be between 0 and 59")
                tries = tries+1
                buffer = False
            elif SEC > 59 or SEC < 0:
                io.errmsg("Enter a valid second please! Second must be between 0 and 59")
                tries = tries+1
                buffer = False
            else:
                buffer = True
                return date

"""
The function determines whether the stopping date entered by the user is in a valid format or not as well as whether its
before or after the starting date. The function also determines whether the year, month, date, hour, minute or seconds entered are valid.
The requirement for each argument is as follows:
Year between 2010 and 2020 and after the year of the starting date
Month between 1 and 12 and after the month of the starting date
Day between 1 and 31 and after the day of the starting date
Hour between 0 and 23 and after the hour of the starting date
Minute between 0 and 59 and after the minute of the starting date
Second between 0 and 59 and after the second of the starting date

Arguments:
1. Date  (Stopping date)
2. Date2 (Starting date)

Returns:
- Validated stopping date
    
Identification: 
    author: Feras Yahya
"""  
def validateStopTime(date,date2):
    buffer = False
    tries = 0
    notValidated = True
    YR2 = int(date2[0:4])
    MO2 = int(date2[5:7])
    D2  = int(date2[8:10])
    HR2 = int(date2[11:13])
    MIN2= int(date2[14:16])
    SEC2= int(date2[17:19])
    
    while buffer != True:
        if tries > 0:
            date = input("TRY AGAIN (Format: YYYY-MM-DD HR:MM:SS)\n")
            notValidated = True
            
        if notValidated != False:
            try:
                YR = int(date[0:4])
                MO = int(date[5:7])
                D  = int(date[8:10])
                HR = int(date[11:13])
                MIN= int(date[14:16])
                SEC= int(date[17:19])
                notValidated = False
            except:
                io.errmsg("Invalid format!")
                tries = tries+1
                buffer = False
                notValidated = True
                
        if notValidated != True:
            if YR > 2020 or YR < 2010 or YR < YR2:
                io.errmsg("Enter a valid year please! Year must be between 2010 and 2020 and must be after starting date")
                tries = tries+1
                buffer = False
            elif MO > 12 or MO < 1 or MO < MO2:
                io.errmsg("Enter a valid month please! Month must be between 1 and 12 and must be after starting date")
                tries = tries+1
                buffer = False
            elif D > 31 or D < 1 or D < D2:
                io.errmsg("Enter a valid day please! Day must be between 1 and 31 and must be after starting date")
                tries = tries+1
                buffer = False
            elif HR > 23 or HR < 0 or HR < HR2:
                if HR > 23 or HR < 0:
                    io.errmsg("Enter a valid hour please! Hour must be between 0 and 23")
                    tries = tries+1
                    buffer = False
                elif YR == YR2 and MO == MO2 and D == D2:
                    io.errmsg("Enter a valid hour please! Hour must be after the hour of the starting date")
                    tries = tries+1
                    buffer = False
                else:
                    buffer = True
                    return date 
            elif MIN > 59 or MIN < 0 or MIN < MIN2:
                if MIN > 59 or MIN < 0:
                    io.errmsg("Enter a valid minute please! Minute must be between 0 and 59 and must be after starting date")
                    tries = tries+1
                    buffer = False
                elif YR == YR2 and MO == MO2 and D == D2 and HR == HR2:
                    io.errmsg("Enter a valid minute please! Minutes must be after the number of minutes of the starting date")
                    tries = tries+1
                    buffer = False
                else:
                    buffer = True
                    return date
            elif SEC > 59 or SEC < 0 or SEC < SEC2:
                if SEC > 59 or SEC < 0:
                    io.errmsg("Enter a valid second please! Second must be between 0 and 59 and must be after starting date")
                    tries = tries+1
                    buffer = False
                elif YR == YR2 and MO == MO2 and D == D2 and HR == HR2 and MIN == MIN2:
                    io.errmsg("Enter a valid number of seconds please! Seconds must be after the number of seconds of the starting date")
                    tries = tries+1
                    buffer = False
                else:
                    buffer = True
                    return date
            else:
                buffer = True
                return date          
                






