import numpy as np
import math

"""
The function calculates the link budget of the satellite

Arguments:
1. F_cnt - Frequency band center (MHz)
2. AE    - Antenna Efficieny
3. D     - Diameter (m)
4. B     - Bandwidth (MHz)
5. RG    - RCV gain (dB)
6. RNT   - Noise temperature (K)
7. Pt    - Power Transmitted
8. Gt    - Gain of transmitting antenna (Onboard satellite)
9. Ltranspath - RF Losses in trasmitter path
10.Latm - Atmospheric and polarization losses
11.rangeTopo - Topocentric range

Returns:
- Power received by satellite
    
Identification: 
    author: Feras Yahya
"""
def LinkDesign(F_cnt, AE, D, B, RG, RNT, Pt, Gt,Ltranspath, Latm, rangeTopo):
    #Calculate power in decibels
    Ptdb = 10*np.log10(Pt)

    #Calculate antenna gain
    FGHz = F_cnt/1000
    Gr = 10*np.log10(110*AE*(FGHz)**2*(D)**2)
    
    #EIRP
    EIRP = Ptdb + Gt - Ltranspath
   
    #Find wavelength
    c = 299792458
    lam = c/(F_cnt*1000000)     #In meters
    

    #Find free space loss
    Lf = ((4*math.pi*(rangeTopo*1000))/lam)**2
    LfdBw = 10*np.log10(Lf)
    
    #Power received
    Pr = EIRP - LfdBw + Gr - Latm   #In dBW
    PrdBm = Pr + 30
   
    #Noise power density
    N0 = -228.6 + 10*np.log10(RNT)

    #Carrier to Noise Density ratio
    C_N0 = Pr - N0

    #Bandwidth
    B = 10*np.log10(B*1000000)

    #Carrier to noise temperature
    C_T = C_N0 - 228.6

    #Carrier to Noise ratio
    C_N = C_N0 - B
 

    return PrdBm
