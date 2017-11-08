function []  = site_DATA(TMY_file_name,latitude,altitude,optimization)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global  DATA...
%     INPUT...
%     field_table...
%     dishNreceiver_table...
%     ducts_size_table...
%     annulus_insulation_table....
%     pipes_insulation_table...
%     insulation_table...
%     HEX_table...
%     plant_table... % optimization V
%     field_Xb_table...
%     annulus_Xb_table... 
%     pipes_Xb_table...
%     HTS_Xb_table...
%     constraints_table...
%     model_type...
%     location...
%     plant_type...
%     system_pressure...
%     HTS_type...
%     file_name...
    
dUTC = 0;
longitude = 0;

% TMY_file_name=uigetfile('*.xls','select a TMY excel sheet data to run'); % target file aquazition
[num,txt,raw]=xlsread(TMY_file_name);% TMY matrix rendered from xl

TMY_length = size(num,1);

TMY1=txt(2:TMY_length+1,2);
TMY2=num(1:TMY_length,3);
TMY3=num(1:TMY_length,4);
TMY4=num(1:TMY_length,5);
TMY5=num(1:TMY_length,6);
TMY6=num(1:TMY_length,7);

Time = time(TMY1);
[Azimuth,Elevation] = sun_position(Time,dUTC,latitude,longitude,altitude);
    
for i = 1:TMY_length-2
    dev_time(i,:) = (Time(i,:)-Time(i+1,:));

%     if i>1 && dev_time(i-1)~=dev_time(i) && dev_time(i)~=22
%         msgbox(['TMY time steps is inconsistant - please see line number' num2str(i) ' : ' TMY1(i) ', of TMY (source) file.']);
%         return
%     end
    
end

step_type = dev_time(1,:)~=0;
seconding = [(3600.*24.*365),(3600.*24.*30),(3600.*24),3600,60,1];
seconds = seconding(step_type);
common = abs(mode(dev_time(dev_time~=0)));
step_duration = common.*seconds;
step_duration(size(step_duration,2)<=1) = 3600 % for DNI short run - default is 1 hour

    if optimization == 0

        DATA.time = TMY1;
        DATA.Azimuth = Azimuth;
        DATA.Elevation = Elevation;
        DATA.DNI = TMY2;
        DATA.wind_velocity = TMY3;
        DATA.wind_direction = TMY4;
        DATA.ambient_temp = TMY5;
        DATA.humidity = TMY6;
        DATA.step_duration = step_duration;

    else
        
        %
        histogram_accuracy = 5;%
        %
        Azimuth_21_6 = (Azimuth(Time(:,2)==6 & Time(:,3)==21));
        Elevation_21_6 = (Elevation(Time(:,2)==6 & Time(:,3)==21));
        Azimuth_21_12 = (Azimuth(Time(:,2)==12 & Time(:,3)==21));
        Elevation_21_12 = (Elevation(Time(:,2)==12 & Time(:,3)==21));
        [n,bin] = hist(TMY2(TMY2>350),histogram_accuracy);
        common_DNI = bin'; 
        DNI = common_DNI;
        DNI_commonality = n'; 
        DNI_common = DNI_commonality;       
        
        Elevation1 = Elevation_21_6(Elevation_21_6>0);
        ind = find(Elevation1==max(Elevation1));
        Elevation1 = Elevation1(1:ind);
        Azimuth1 = Azimuth_21_6(Elevation_21_6>0);
        Azimuth1 = Azimuth1(1:ind);
        
        Elevation2 = Elevation_21_12(Elevation_21_12>0);
        ind = find(Elevation2==max(Elevation2));
        Elevation2 = Elevation2(1:ind);
        Azimuth2 = Azimuth_21_12(Elevation_21_12>0);
        Azimuth2 = Azimuth2(1:ind);
        
        AzimuthX = [Azimuth1;Azimuth2];
        ElevationX = [Elevation1;Elevation2];

        for i = 1:length(common_DNI).*length(ElevationX)
            XAzimuth(i) = AzimuthX(ceil(i./length(common_DNI)));
            XElevation(i) = ElevationX(ceil(i./length(common_DNI)));
            if i./length(ElevationX)>1 && i./length(ElevationX)<2
                DNI = [DNI;common_DNI];
                DNI_common = [DNI_common;DNI_commonality];
            end
        end     
        
        DATA.Azimuth = XAzimuth';
        DATA.Elevation = XElevation';
        DATA.DNI = DNI;
        DATA.DNI_commonality = DNI_common;        
        DATA.wind_velocity = mean(TMY3(TMY2>350)).*ones(length(XAzimuth),1);
        DATA.wind_direction = mean(TMY4(TMY2>350)).*ones(length(XAzimuth),1);
        DATA.ambient_temp = mean(TMY5(TMY2>350)).*ones(length(XAzimuth),1);
        DATA.humidity = mean(TMY6(TMY2>350)).*ones(length(XAzimuth),1);
        DATA.step_duration = step_duration;
        
    end

end

function Time = time(time_string)

        tt1=datevec(time_string,'dd/mm/yyyy');
        tt2=datevec(time_string); 
            year = tt1(:,1);
            months = tt1(:,2);
            days = tt1(:,3);
            hours = tt2(:,4);
            minutes = tt2(:,5);
            seconds = tt2(:,6);

Time = [year,months,days,hours,minutes,seconds];

end


%% SUN POSITION:

function [Az,El] = sun_position(Time,dUTC,Lat,Lon,Alt) %UTC

% Programed by Darin C. Koblick 2./17./2009

% External Function Call Sequence:
%[Az El] = SolarAzEl('1991./05./19 13:00:00',50,10,0)

% Function Description:
% SolarAzEl will ingest a Universal Time, and specific site location on earth
% it will then output the solar Azimuth and Elevation angles relative to that
% site.

%Input Description:
% UTC (Coordinated Universal Time YYYY./MM./DD hh:mm:ss)
% Lat (Site Latitude in degrees -90:90 -> S(-) N(+))
% Lon (Site Longitude in degrees -180:180 W(-) E(+))
% Altitude of the site above sea level (km)

%Output Description:
%Az (Azimuth location of the sun in degrees)
%El (Elevation location of the sun in degrees)

%Source References:
%Solar Position obtained from:
%http:././stjarnhimlen.se./comp./tutorial.html#5

% Code Sequence

%compute JD
jd = juliandate(Time,dUTC); %UTC,'yyyy./mm./dd HH:MM:SS'
d = jd-2451543.5;

% Keplerian Elements for the Sun (geocentric)
w = 282.9404+4.70935e-5.*d; %    (longitude of perihelion degrees)
a = 1.000000;%                  (mean distance, a.u.)
e = 0.016709-1.151e-9.*d;%       (eccentricity)
M = mod(356.0470+0.9856002585.*d,360);%   (mean anomaly degrees)
L = w + M;                     %(Sun's mean longitude degrees)
oblecl = 23.4393-3.563e-7.*d;  %(Sun's obliquity of the ecliptic)

%auxiliary angle
E = M+(180./pi).*e.*sin(M.*(pi./180)).*(1+e.*cos(M.*(pi./180)));

%rectangular coordinates in the plane of the ecliptic (x axis toward
%perhilion)
x = cos(E.*(pi./180))-e;
y = sin(E.*(pi./180)).*sqrt(1-e.^2);

%find the distance and true anomaly
r = sqrt(x.^2 + y.^2);
v = atan2(y,x).*(180./pi);

%find the longitude of the sun
lon = v + w;

%compute the ecliptic rectangular coordinates
xeclip = r.*cos(lon.*(pi./180));
yeclip = r.*sin(lon.*(pi./180));
zeclip = 0.0;

%rotate these coordinates to equitorial rectangular coordinates
xequat = xeclip;
yequat = yeclip.*cos(oblecl.*(pi./180))+zeclip.*sin(oblecl.*(pi./180));
zequat = yeclip.*sin(23.4406.*(pi./180))+zeclip.*cos(oblecl.*(pi./180));

%convert equatorial rectangular coordinates to RA and Decl:
r = sqrt(xequat.^2 + yequat.^2 + zequat.^2)-(Alt./149598000); %roll up the altitude correction
RA = atan2(yequat,xequat).*(180./pi);
delta = asin(zequat./r).*(180./pi);

%Following the RA DEC to Az Alt conversion sequence explained here:
%http:././www.stargazing.net./kepler./altaz.html

%Find the J2000 value
J2000 = jd - 2451545.0;
hourvec = Time;%datevec(UTC,'yyyy./mm./dd HH:MM:SS');
UTH = hourvec(:,4) + hourvec(:,5)./60 + hourvec(:,6)./3600;

%Calculate local siderial time
GMST0=mod(L+180,360)./15;
SIDTIME = GMST0 + UTH + Lon./15;

%Replace RA with hour angle HA
HA = (SIDTIME.*15 - RA);

%convert to rectangular coordinate system
x = cos(HA.*(pi./180)).*cos(delta.*(pi./180));
y = sin(HA.*(pi./180)).*cos(delta.*(pi./180));
z = sin(delta.*(pi./180));

%rotate this along an axis going east-west.
xhor = x.*cos((90-Lat).*(pi./180))-z.*sin((90-Lat).*(pi./180));
yhor = y;
zhor = x.*sin((90-Lat).*(pi./180))+z.*cos((90-Lat).*(pi./180));

%Find the h and AZ 
Az = (atan2(yhor,xhor)).*(180./pi) + 180;
El = asin(zhor).*(180./pi);

% R = (6400-6380).*Lat./90+6380;
% Az =  180.*(1 - sind(15.*(UTH-12)).*cosd(delta)./cosd(El));%.*(sind(El));%.*((R+Alt)./(2.*R)));
%
% Az(Az>360) = Az(Az>360)-360;
% Az(Az<0) = 360-Az(Az<0);

% 
% month = Time(:,2);
% day = Time(:,3) + floor( 30.6001.*(month + 1.0));
% hour = Time(:,4) - 12;
% ds = 23.45.*sind(360.*(284+day)./365);
% hs = 15.*hour;
% Az = 180-(180.*(cosd(ds).*sind(hs)./cosd(El)));
%
% if 12-abs(month)<=6
%     Az = abs(180-Az);
% end
% 
% if latitude>=0
%     Az  = 180 - asind(cosd(ds).*sind(hs)./cosd(elevation));
% else
%     Az  = - 180 - asind(cosd(ds).*sind(hs)./cosd(elevation));
% end

%sun = [Az,El];

end
function jd = juliandate(Time,dUTC) %varargin
    % This sub function is provided in case juliandate does not come with your 
    % distribution of Matlab

    %[year,month,day,hour,min,sec] = datevec(datenum(varargin{:}));
    year = Time(:,1);
    month = Time(:,2);
    day = Time(:,3);
    hour = Time(:,4);
    min = Time(:,5);
    sec = Time(:,6);

    
    
    
    for k = length(month):-1:1
        if ( month(k) <= 2 ) % january & february
            year(k)  = year(k) - 1.0;
            month(k) = month(k) + 12.0;
        end
    end

    jd = floor( 365.25.*(year + 4716.0)) + floor( 30.6001.*(month + 1.0)) + 2.0 - ...
        floor(year./100.0) + floor( floor(year./100.0 )./4.0) + day - 1524.5 + ...
        (hour+ dUTC + min./60 + sec./3600)./24;

    end

