package SatellitePackage//228.868

  model Satellite
    //import Modelica.SIunits.Conversions.*;
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;
    import Modelica.Constants.*;
              /*
              parameter Real ecc = 0.0034018"Eccentricity";
              parameter Real M0 = 267.4987"Mean anomaly at Epoch (deg)";
              parameter Real N0 = 2.005652661"Mean motion at Epoch (rev/d)";
              parameter Real Ndot2 = 0.00000029"TLE drag parameter rev/d^2";
              parameter Real Nddot6 = 0"TLE parameter rev/d^3";
              parameter Real tstart = 937743.43392"Simulation start time, seconds since Epoch (s)";
              */
    parameter Real ecc "Eccentricity";
    parameter Real M0 "Mean anomaly at Epoch (deg)";
    parameter Real N0 "Mean motion at Epoch (rev/d)";
    parameter Real Ndot2 "TLE drag parameter rev/d^2";
    parameter Real Nddot6 "TLE parameter rev/d^3";
    parameter Real tstart "Simulation start time, seconds since Epoch (s)";
    Real M "Mean Anomaly (deg)";
    Real E "Eccentric Anomaly (deg)";
    Real theta "True Anomaly (deg)";
    Real r "Range (km)";
    Real a "Semi-major axis (km)";
   
    constant Real mu = 398600.4418 "Gravitational Parameter for Earth";
    Vector p_sat_pf "Position, Perifocal coords";
    Vector v_sat_pf "Velocity Peerifocal coords";
  equation
  
//Compute M
    M = mod(M0 + N0 * (360 * (time + tstart) / 86400) + 360 * Ndot2 * ((time + tstart) / 86400) ^ 2 + 360 * Nddot6 * ((time + tstart) / 86400) ^ 3, 360);
  
  //Compute E
    E = from_deg(M) + ecc * sin(E);
//compute theta
    tan(theta / 2) = ((1 + ecc) / (1 - ecc)) ^ (1 / 2) * tan(E / 2);
//find a
    a = (mu / (N0 * (360 / 86400) * pi / 180) ^ 2) ^ (1 / 3);
//find range
    r = a * (1 - ecc ^ 2) / (1 + ecc * cos(theta));

//find position x ,y ,z
    p_sat_pf.x = r * cos(theta);
    p_sat_pf.y = r * sin(theta);
    p_sat_pf.z = 0;
//find velocity xdot, ydot, zdot
    v_sat_pf.x = der(p_sat_pf.x);
    v_sat_pf.y = der(p_sat_pf.y);
    v_sat_pf.z = 0;
  end Satellite;








  model GndStn
    import Modelica.Constants.*;
          /*
          parameter Real longitude = 281.927 "Station longitude (deg)";
          parameter Real latitude = 45.9555 "Station latitude (deg)";
          parameter Real elevation = 2604.2 "Station elevation (m)";
          */
    parameter Real longitude "Station longitude (deg)";
    parameter Real latitude "Station latitude (deg)";
    parameter Real elevation "Station elevation (m)";
    Real N, e2;
    Real a = 6378137 "Ellipsoid EquatorialRadius (m) WGS84";
    Real f = 1 / 298.257223563 " Ellipsoid Flattening";
    output Vector p_stn_ecf "Station ECF coordinates";
  equation
    e2 = 2 * f - f ^ 2 "ellipsoidal eccentricity ^2";
    N = a / (1 - e2 * sin(latitude * (pi / 180)) ^ 2) ^ (1 / 2) "Ellisoidal radius of curvature of meridian";
    p_stn_ecf.x = (N + elevation) * cos(latitude * (pi / 180)) * cos(longitude * (pi / 180)) / 1000;
    p_stn_ecf.y = (N + elevation) * cos(latitude * (pi / 180)) * sin(longitude * (pi / 180)) / 1000;
    p_stn_ecf.z = ((1 - e2) * N + elevation) * sin(latitude * (pi / 180)) / 1000;
  end GndStn;


  function range_ECF2topo
    import Modelica.Constants.*;
    input Vector p_sat_ecf "Satellite position in ECF coords (km)";
    input Vector v_sat_ecf "Satellite velocity in ECF coords (km/s)";
    input Vector p_stn_ecf "Station position in ECF coords (km)";
    input Real longitude "Station longitude (degE)";
    input Real latitude "Station latitude (deg)";
    output Vector p_sat_topo "Satellite position in topo coords (km)";
    output Vector v_sat_topo "Satellite velocity in topo coords (km/s)";
  protected
    /*Real RM[3,3] = [-sin(longitude*pi/180), cos(longitude*pi/180), 0; -cos(longitude*pi/180)*sin(latitude*pi/180), -sin(longitude*pi/180)*sin(latitude*pi/180), cos(latitude*pi/180); cos(longitude*pi/180)*cos(latitude*pi/180), sin(longitude*pi/180)*cos(latitude*pi/180), sin(latitude*pi/180)]; */
    Real RM[3, 3] = [-sin(longitude * pi / 180), cos(longitude * pi / 180), 0; -cos(longitude * pi / 180) * sin(latitude * pi / 180), -sin(longitude * pi / 180) * sin(latitude * pi / 180), cos(latitude * pi / 180); cos(longitude * pi / 180) * cos(latitude * pi / 180), sin(longitude * pi / 180) * cos(latitude * pi / 180), sin(latitude * pi / 180)];
    Real xecf = p_sat_ecf.x;
    Real yecf = p_sat_ecf.y;
    Real zecf = p_sat_ecf.z;
    Real xdpecf = v_sat_ecf.x;
    Real ydpecf = v_sat_ecf.y;
    Real zdpecf = v_sat_ecf.z;
    Real xstn = p_stn_ecf.x;
    Real ystn = p_stn_ecf.y;
    Real zstn = p_stn_ecf.z;
  algorithm
    p_sat_topo.x := (RM[1, 1] * (xecf-xstn) + RM[1, 2] * (yecf-ystn) + RM[1, 3] * (zecf-zstn));
    p_sat_topo.y := (RM[2, 1] * (xecf-xstn) + RM[2, 2] * (yecf-ystn) + RM[2, 3] * (zecf-zstn));
    p_sat_topo.z := (RM[3, 1] * (xecf-xstn) + RM[3, 2] * (yecf-ystn) + RM[3, 3] * (zecf-zstn));
    
    v_sat_topo.x := (RM[1, 1] * xdpecf + RM[1, 2] * ydpecf + RM[1, 3] * zdpecf);
    v_sat_topo.y := (RM[2, 1] * xdpecf + RM[2, 2] * ydpecf + RM[2, 3] * zdpecf);
    v_sat_topo.z := (RM[3, 1] * xdpecf + RM[3, 2] * ydpecf + RM[3, 3] * zdpecf);
  end range_ECF2topo;






function range_topo2look_angles
  import Modelica.Math.Vectors.*;
  import Modelica.Constants.*;
  input Real az_vel_lim "Azimuth velocity limit (deg/s)";
  input Real el_vel_lim "Elevation velocity limit (deg/s)";
  input Vector p_sat_topo "Satellite position in topo coords (km)";
  input Vector v_sat_topo "Satellite velocity in topo coords (km/s)";
  output Real azimuth "Azimuth angle (deg)";
  output Real elevation "Elevation angle (deg)";
protected
  Real Rxy[1, 3] = [p_sat_topo.x, p_sat_topo.y, 0];
  Real Vxy[1, 3] = [v_sat_topo.x, v_sat_topo.y, 0];
  Real norm_Rxy = (Rxy[1, 1] ^ 2 + Rxy[1, 2] ^ 2) ^ 0.5;
  Real norm_psatopo = (p_sat_topo.x ^ 2 + p_sat_topo.y ^ 2 + p_sat_topo.z ^ 2) ^ 0.5;
  Real cross_VR[1, 3] = [Vxy[1, 2] * Rxy[1, 3] - Vxy[1, 3] * Rxy[1, 2], Vxy[1, 3] * Rxy[1, 1] - Vxy[1, 1] * Rxy[1, 3], Vxy[1, 1] * Rxy[1, 2] - Vxy[1, 2] * Rxy[1, 1]];
  Real dAz[1,3] = 1 / norm_Rxy ^ 2 * cross_VR;
  Real dotProduct = Rxy[1, 1] * Vxy[1, 1] + Rxy[1, 2] * Vxy[1, 2];
  Real dEl = 1 / norm_psatopo ^ 2 * (norm_Rxy * v_sat_topo.z - p_sat_topo.x / norm_Rxy * dotProduct);
algorithm
  if dAz[1,1] > az_vel_lim or dEl > el_vel_lim then
    azimuth := 0;
    elevation := 0;
  else
    azimuth := atan2(p_sat_topo.x, p_sat_topo.y);
    if azimuth <0 then
      azimuth := (2*pi)+azimuth;
    end if;
    elevation := atan(p_sat_topo.z / (p_sat_topo.x ^ 2 + p_sat_topo.y ^ 2) ^ 0.5);
  end if;
end range_topo2look_angles;


























  
  


  function sat_ECF
    import Modelica.Constants.*;
    input Vector p_sat_eci "Satellite position, ECI coords (km)";
    input Vector v_sat_eci "Satellite velocity, ECI coords (km/s)";
    /*
          input Real xpf;
          input Real ypf;
          input Real zpf;
            
          input Real xdpf;
          input Real ydpf;
          input Real zdpf;
          */
    input Real theta_t "GMST (deg)";
    output Vector p_sat_ecf "Satellite position, ECF coords (km)";
    output Vector v_sat_ecf "Satellite velocity, ECF coords (km/s)";
  protected
    /*
          Real side_time = 23*3600 + 56*60 + 4.091; // sidereal time
          Real sidereal_rot_rate = 2*pi/side_time; //sidereal rotation rate
          */
    //Real omega_e = 7.29211585275553 * 10 ^ (-5);
    Real omega_e = (360/86164.091) *(pi/180);
    Real rotT[3, 3] = [cos(theta_t * pi / 180), sin(theta_t * pi / 180), 0; -sin(theta_t * pi / 180), cos(theta_t * pi / 180), 0; 0, 0, 1];
    Real rot2[3, 3] = [-sin(theta_t * pi / 180), cos(theta_t * pi / 180), 0; -cos(theta_t * pi / 180), -sin(theta_t * pi / 180), 0; 0, 0, 0];
    //Real v[3,1];
    Real xpf = p_sat_eci.x;
    Real ypf = p_sat_eci.y;
    Real zpf = p_sat_eci.z;
    Real xdpf = v_sat_eci.x;
    Real ydpf = v_sat_eci.y;
    Real zdpf = v_sat_eci.z;
  algorithm
// v := rotT* [xdpf;ydpf;zdpf] - (-sidereal_rot_rate)*rot2* [xpf;ypf;zpf];
    p_sat_ecf.x := rotT[1, 1] * xpf + rotT[1, 2] * ypf + rotT[1, 3] * zpf;
    p_sat_ecf.y := rotT[2, 1] * xpf + rotT[2, 2] * ypf + rotT[2, 3] * zpf;
    p_sat_ecf.z := rotT[3, 1] * xpf + rotT[3, 2] * ypf + rotT[3, 3] * zpf;
  /*
  v_sat_ecf.x := (rotT[1,1]*xdpf + rotT[1,2]*ydpf + rotT[1,3]*zdpf) - (rot2[1,1]*xpf*(-sidereal_rot_rate) + rot2[1,2]*ypf*(-sidereal_rot_rate) + rot2[1,3]*zpf*(-sidereal_rot_rate));
  v_sat_ecf.y := (rotT[2,1]*xdpf + rotT[2,2]*ydpf + rotT[2,3]*zdpf) - (rot2[2,1]*xpf*(-sidereal_rot_rate) + rot2[2,2]*ypf*(-sidereal_rot_rate) + rot2[2,3]*zpf*(-sidereal_rot_rate));
  v_sat_ecf.z := (rotT[3,1]*xdpf + rotT[3,2]*ydpf + rotT[3,3]*zdpf) - (rot2[3,1]*xpf*(-sidereal_rot_rate) + rot2[3,2]*ypf*(-sidereal_rot_rate) + rot2[3,3]*zpf*(-sidereal_rot_rate));
  */
    v_sat_ecf.x := rotT[1, 1] * xdpf + rotT[1, 2] * ydpf + rotT[1, 3] * zdpf - (rot2[1, 1] * xpf * (-omega_e) + rot2[1, 2] * ypf * (-omega_e) + rot2[1, 3] * zpf * (-omega_e));
    v_sat_ecf.y := rotT[2, 1] * xdpf + rotT[2, 2] * ydpf + rotT[2, 3] * zdpf - (rot2[2, 1] * xpf * (-omega_e) + rot2[2, 2] * ypf * (-omega_e) + rot2[2, 3] * zpf * (-omega_e));
    v_sat_ecf.z := rotT[3, 1] * xdpf + rotT[3, 2] * ydpf + rotT[3, 3] * zdpf - (rot2[3, 1] * xpf * (-omega_e) + rot2[3, 2] * ypf * (-omega_e) + rot2[3, 3] * zpf * (-omega_e));
/*
 v_sat_ecf.x := v[1,1];
  v_sat_ecf.y := v[2,1];
  v_sat_ecf.z := v[3,1];
  */
  end sat_ECF;



  function sat_ECI
    import Modelica.Constants.*;
    input Vector p_sat_pf "Satellite position, PF coords (km)";
    input Vector v_sat_pf "Satellite velocity, PF coords (km/s)";
    /*
          input  Real xpf;
          input Real ypf;
          input Real zpf;
          
          input Real xdpf;
          input Real ydpf;
          input Real zdpf;
          */
    input Real ecc "Eccentricity";
    input Real raan "Right Ascension of ascending node (deg)";
    input Real inc "Inclination angle (deg)";
    input Real argper "Argument of Perigee (deg)";
    input Real N "Mean Motion (rev/d)";
    output Vector p_sat_eci "Satellite position ECI coords (km)";
    output Vector v_sat_eci "Satellite velocity ECI coords (km/s)";
  protected
    Real rotw[3, 3] = [cos(-argper * pi / 180), sin(-argper * pi / 180), 0; -sin(-argper * pi / 180), cos(-argper * pi / 180), 0; 0, 0, 1];
    Real roti[3, 3] = [1, 0, 0; 0, cos(-inc * pi / 180), sin(-inc * pi / 180); 0, -sin(-inc * pi / 180), cos(-inc * pi / 180)];
    Real rotr[3, 3] = [cos(-raan * pi / 180), sin(-raan * pi / 180), 0; -sin(-raan * pi / 180), cos(-raan * pi / 180), 0; 0, 0, 1];
    Real rotT_p[3, 3] = rotr * roti * rotw;
    //this is the transformation matrix for ECI to perifocal.
    //To get transformation matrix for perifocal to ECI, we take the transpose
    Real rotT[3, 3] = rotT_p;
    Real xpf = p_sat_pf.x;
    Real ypf = p_sat_pf.y;
    Real zpf = p_sat_pf.z;
    Real xdpf = v_sat_pf.x;
    Real ydpf = v_sat_pf.y;
    Real zdpf = v_sat_pf.z;
  algorithm
    p_sat_eci.x := rotT[1, 1] * xpf + rotT[1, 2] * ypf + rotT[1, 3] * zpf;
    p_sat_eci.y := rotT[2, 1] * xpf + rotT[2, 2] * ypf + rotT[2, 3] * zpf;
    p_sat_eci.z := rotT[3, 1] * xpf + rotT[3, 2] * ypf + rotT[3, 3] * zpf;
    v_sat_eci.x := rotT[1, 1] * xdpf + rotT[1, 2] * ydpf + rotT[1, 3] * zdpf;
    v_sat_eci.y := rotT[2, 1] * xdpf + rotT[2, 2] * ydpf + rotT[2, 3] * zdpf;
    v_sat_eci.z := rotT[3, 1] * xdpf + rotT[3, 2] * ydpf + rotT[3, 3] * zdpf;
  end sat_ECI;

  /*function theta_t
        import Modelica.Constants.*;
        import Modelica.SIunits.Conversions.*;
        input Real t "Time since J2000 (d)";
        output Real theta_t "GMST angle (deg)";
        
        
        
      protected
        Real t_mid = mod(t,0.5);
        Real Du = t - t_mid;
        
        Real Tu = Du/36525;
        
        Real GMST = mod((24110.5484 + (8640184*Tu) + (0.093104*Tu^2) - (6.2*10^(-6)*Tu^3)), 86400);
        
      
         
        
        Real theta_mid = 360*(GMST/86400);
        
        Real r = 1.002737909350795 + (5.9006*10^(-11)*Tu) - 5.9*10^(-15)*Tu^2;
        
      algorithm
        theta_t := mod(280.46061837 + (360.98564736629*t)+(0.0003875*Tu^2)-(2.6*10^(-8)*Tu^3), 360); 
        
        
        
      end theta_t;*/
  /*
      function theta_t
        input Real din "JD of the time of interst (not time difference since J2000)";
        //input Real t "time of interest";
        input Real tm "time of mindnight of day of interest";
        output Real theta_t "GMST angle (deg)";
      protected
        Real d = din - 2451545.0;
       // Real d = 1879.386107;
        Real T = d/36525;
        Real theta_mid;
        Real r;
      algorithm
        theta_mid := mod(24110.54841 + 8640184.812866*T +0.093104*T^2 - 0.0000062*T^3,86400)*360/86400;
        r:=1.002737909350795 + 5.9006*10^(-11)*T -5.9*10^(-15)*T^2;
        //theta_t := theta_mid + 360*r*(2453424.416366-2453423.500000)/86400;//(time of interest ... in this case start of tracking time - tmidnight of start of tracking period)
        theta_t := mod(theta_mid + 360*r*((din-tm)*86400)/86400,360);//(time of interest ... in this case start of tracking time - tmidnight of start of tracking period)
      end theta_t;
      */

  record Vector
    Real x, y, z;
  end Vector;

  model test
  import Modelica.Constants.*;
    parameter Real JulianDate "time since J2000= 2458241.8445833335 ";
    parameter Real hr0 "Time since midnight in hours";     
    parameter Real raan;
    parameter Real inc;
    parameter Real argper;
    
    Vector Position;
    Vector Velocity;
    Vector PositionECF;
    Vector VelocityECF;
    Vector PositionTOPO;
    Vector VelocityTOPO;
    
    Real Azimuth;
    Real Elevation;
    Real GMST;
  
    Real hr;
    Satellite GPS_Test "(ecc=0.0032278,M0=268.978,N0=2.00565,Ndot2=0.00000002,Nddot6=0,tstart=188680.32288)";
    GndStn ARO "(longitude = 281.927 , latitude = 45.9555 , elevation = 2604.2)";
  equation
    (Position, Velocity) = sat_ECI(GPS_Test.p_sat_pf, GPS_Test.v_sat_pf, GPS_Test.ecc, raan, inc, argper, GPS_Test.N0);
    hr = hr0 + time/3600;
    GMST = theta_t(JulianDate, hr);
  
    (PositionECF, VelocityECF) = sat_ECF(Position, Velocity, GMST);
    (PositionTOPO, VelocityTOPO) = range_ECF2topo(PositionECF, VelocityECF, ARO.p_stn_ecf, ARO.longitude, ARO.latitude);
    (Azimuth,Elevation)=range_topo2look_angles(6,7,PositionTOPO,VelocityTOPO);
    annotation(
      startTime = 0,
      stopTime = 86400,
      stepSize = 500);
  end test;































function theta_t

  //input Real t_start "Time since J2000 (d)";
  //input Real dt "time";
  input Real days "Time since J2000";
  input Real hours "Hours since midnight UTC";
  output Real GMST "GMST angle (deg)";
protected
  //Real t_since_mid, Du, Tu, theta0, theta0_rad, r;
  Real Tu;
  
algorithm
// Time since midnight on the same date (in days):
 // t_since_mid := mod(t_start, 0.5);
// Time (d) passed between J2000 and midnight on the same date (in days);
  //Du := ((dt/60/60/24)+t_start) - t_since_mid;
// Time in Julian centuries:
  Tu := days / 36525;
// GMST at 00h UT on the day of start of observations (sec):
 // theta0 := mod(23991.84441 + 8640185.812866 * Tu + 0.093104 * Tu ^ 2 - 6.2 * 10 ^ (-6) * Tu ^ 3, 86400);
// Convert it to Degrees:
  //theta_mid := 360 * theta0 / 86400;
//Ratio of an interval in Sid seconds to sae interval in UT seconds:
//  r := 1.002737909350795 + 5.9006 * 10 ^ (-11) * Tu - 5.9 * 10 ^ (-15) * Tu ^ 2;
// GMST at current time:
//theta := (180/pi)*mod(theta0_rad + 360*r*(Du),2*pi);
//  theta := mod(theta_mid + 360 *r*((dt+(tstart*24*60*60))-t_since_mid), 360);
GMST := mod(mod(23991.84441 + 8640185.812866 * Tu + 0.093104 * Tu ^ 2 - 6.2 * 10 ^ (-6) * Tu ^ 3, 86400)*(360/86400) + 4.178074163e-3*hours*3600., 360.);

end theta_t;















end SatellitePackage;
