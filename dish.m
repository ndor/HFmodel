function DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA)
%DISH4OPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here

% time = DATA.time;
azimuth = DATA.Azimuth;
elevation = DATA.Elevation;
DNI = DATA.DNI;
wind_velocity = DATA.wind_velocity;
wind_direction = DATA.wind_direction;
% ambient_temp = DATA.ambient_temp;
% humidity = DATA.humidity;
% step_duration = DATA.step_duration;

a = 22;%sqrt(2*17.3^2); % single dish facet length [m]
N_dishes = N_dishes_per_cluster.*N_clusters_in_field;

% cleaning parameters:
    time_length = size(DATA.DNI,1);
    min_dirt_eff = 0.96;
    t_to_cleaning = 7*24;
    t_to_primary_cleaning = 4*7*24;
    
Slope_Error_intercept = 0.98.*(-0.0525.*(Slope_Error+0.2).^4 + 1.3333.*(Slope_Error+0.2).^3 - 11.372.*(Slope_Error+0.2).^2 + 27.422.*(Slope_Error+0.2) + 78.08)./100; % only for the range of: 1.5<=Slope_Error<=10

%dish: x-number of dishes in field
COST_steel_construction = @(xi) interp1([0,4,8,100,10000],[0,341000,217800,176000,176000],xi,'pchip');
COST_reciever = @(xi) interp1([0,4,8,100,10000],[0,31300,34750,29000,29000],xi,'pchip');
COST_hydraulics = @(xi) interp1([0,4,8,100,10000],[0,70000,70000,60000,60000],xi,'pchip');
COST_steel_components = @(xi) interp1([0,4,8,100,10000],[0,150000,70000,63000,63000],xi,'pchip');
COST_mirrors = @(xi) interp1([0,4,8,100,10000],[0,100000,100000,72250,72250],xi,'pchip');

    if N_clusters_in_field>1
        Nns = N_dishes_per_cluster;
        New = N_clusters_in_field./2;
    elseif N_clusters_in_field.*N_dishes_per_cluster==1
        Nns = 1;
        New = 1;
    else
        Nns = N_dishes_per_cluster;
        New = 1;
    end                

    El2 = 0; El1 = 0;
    Az2 = 0; Az1 = 0;
    
    for i=1:time_length    
        

        
        if elevation(i)>0

            non_shaded_percent(i) = shading(lns,lew,alpha,Nns,New,field_azimuth_or,field_elevation_or,azimuth(i),elevation(i),a);
            effstruct(i) = structural(azimuth(i),elevation(i),wind_velocity(i),wind_direction(i));
            effdust(i) = dust(i,min_dirt_eff,t_to_cleaning,t_to_primary_cleaning,reflectivity);
            efficiency(i) = 0.96.*non_shaded_percent(i).*effstruct(i).*effdust(i).*Slope_Error_intercept;
            dish_power_pr_sqr_m(i) = efficiency(i).*DNI(i);

        else
            
            non_shaded_percent(i) = 0;
            effstruct(i) = 0;
            effdust(i) = 0;
            efficiency(i) = 0;
            dish_power_pr_sqr_m(i) = 0;
            
        end

                            % Dish hydraulics Module:
                            if i/2 == round(i/2)
                                El2 = elevation(i);
                                Az2 = azimuth(i);
                            else
                                El1 = elevation(i);
                                Az1 = azimuth(i);
                            end

                            if elevation(i)<=0 && wind_velocity(i)<14
                                dEl = 0;
                            else
                                dEl = abs(El2-El1); % elevation movement per hour 
                            end

                            dAz = abs(Az2-Az1); % azimuth movement per hour

                            if abs(abs(Az1)-abs(Az2))>30
                                dAz = 15;
                            end

                            if abs(abs(El1)-abs(El2))>30
                                dAz = 15;
                            end
                            
                            consumption(i) = (1.*dEl./60+0.833.*dAz./60); % to [KW] in [deg/min], times 1 hour is [KWh]
        
    end
    
    
    
    
    
    DISH.Slope_Error_intercept = Slope_Error_intercept; % [%]
    DISH.effdust = effdust'; % [%]
    DISH.non_shaded = non_shaded_percent'; % [%]
    DISH.effstruct = effstruct'; % [%]
    DISH.efficiency = efficiency'; % [%]
    DISH.dish_power_pr_sqr_m = dish_power_pr_sqr_m'; % [W/m2]
% %     [n,bin] = hist(DISH.dish_power_pr_sqr_m(DISH.dish_power_pr_sqr_m>350),30);
% %     DISH.output_statistics = [n;bin]';
    DISH.consumption = consumption'; % [KW] (one dish)
    DISH.COST.steel_construction = COST_steel_construction(N_dishes); % [$/sish]
    DISH.COST.reciever = COST_reciever(N_dishes); % [$/sish]
    DISH.COST.hydraulics = COST_hydraulics(N_dishes); % [$/sish]
    DISH.COST.steel_components = COST_steel_components(N_dishes); % [$/sish]
    DISH.COST.mirrors = COST_mirrors(N_dishes); % [$/sish]
    DISH.COST.total = DISH.COST.steel_construction+DISH.COST.reciever+DISH.COST.hydraulics+DISH.COST.steel_components+DISH.COST.mirrors; % [$/sish]
    
end

%% dish efficiency functions:

function shading_percent = shading(lns,lew,alpha,Nns,New,field_azimuth_or,field_elevation_or,azimuth,elevation,a)
%%  shading calculates the percent of non-shaded dish area in the field (for the given parameters)
%++++++++++++++++    parameters   ++++++++++++++++++++++++++++++++++++++++++++++++++++
%   alpha – field alpha angle (field tilt of parallelogram) [deg]
%   field_azimuth_or - field azimuth orientation (north axis tilt) [Deg]
%   field_elevation_or - field elevation orientation (horizontal axis tilt) [Deg]
%   Nns - number of rows along north south 
%   New - number of rows along west east
%   lns: North to South centers distance [m]
%   lew: East to West dish centers distance [m]
%   azimuth: the sun's Azimuth angle [Deg]
%   elevation: the sun's Elevation angle [Deg]
%   a: quad dish facet length [m]
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if elevation>0

    elevation = elevation + field_elevation_or;
    azimuth = azimuth + field_azimuth_or;

    tg_alpha = tand(alpha);

    S = a.^2;


    B(1) = lew.*0;
    A(1) = lns.*1 - B(1).*tg_alpha;

    B(2) = lew.*0.5;
    A(2) = lns.*1.5 - B(2).*tg_alpha;

    B(3) = lew.*1;
    A(3) = lns.*2 - B(3).*tg_alpha;

    B(4) = lew.*0.5;
    A(4) = lns.*0.5 - B(4).*tg_alpha;

    B(5) = lew.*2;
    A(5) = lns.*1 - B(5).*tg_alpha;

    B(6) = lew.*1.5;
    A(6) = lns.*0.5 - B(6).*tg_alpha;

    B(7) = lew.*1;
    A(7) = lns.*0 - B(7).*tg_alpha;

    B(8) = lew.*1.5;
    A(8) = -lns.*0.5 + B(8).*tg_alpha;

    B(9) = lew.*2;
    A(9) = -lns.*1 + B(9).*tg_alpha;

    B(10) = lew.*0.5;
    A(10) = -lns.*0.5 + B(10).*tg_alpha;

    B(11) = lew.*1;
    A(11) = -lns.*2 + B(11).*tg_alpha;

    B(12) = lew.*0.5;
    A(12) = -lns.*1.5 + B(12).*tg_alpha;

    B(13) = lew.*0;
    A(13) = -lns.*1 + B(13).*tg_alpha;

    Dto = sqrt(A.^2+B.^2);
    Azto = atand(B./A);
    Azto(Azto<0) = 180+Azto(Azto<0);
    Azto(end) = 180;

    D = a./sind(elevation);
    El_shade = (D-Dto).*sind(elevation);
    El_shade(El_shade<0) = 0;

    if azimuth>180; azimuth = azimuth-180; end
    
    Az_shade = a - Dto.*abs(sind(Azto-azimuth));
    
    if abs(Azto-azimuth)>=90; Az_shade(1:end) = 0; end
    
    Az_shade(Az_shade<0) = 0;

            % each shading dish azimuthial / width span on target dish:
        for k = 1:length(Azto)
            
            if  Az_shade(k)>0
                
                if sign(sind(Azto(k)-azimuth))>0 && abs(Azto(k)-azimuth)>=90
                    Az_shade_span(1,k) = a./2;
                    Az_shade_span(2,k) = a./2 - Az_shade(k);
                elseif sign(sind(Azto(k)-azimuth))<0 && abs(Azto(k)-azimuth)>=90
                    Az_shade_span(1,k) = Az_shade(k) - a./2;
                    Az_shade_span(2,k) = -a./2;
                elseif sign(sind(Azto(k)-azimuth))==0 && abs(Azto(k)-azimuth)>=90
                    Az_shade_span(1,k) = a./2;
                    Az_shade_span(2,k) = -a./2;
                else
                    Az_shade_span(1,k) = 0;
                    Az_shade_span(2,k) = 0;
                end
                
            else
                
                Az_shade_span(1,k) = 0;
                Az_shade_span(2,k) = 0;
                
            end
            
        end

        % 1st order shade overlap:
        if Az_shade(1)>0 && Az_shade(2)>0
            p1 = Az_shade_span(1,1);
            p2 = Az_shade_span(2,1);
            s1 = Az_shade_span(1,2);
            s2 = Az_shade_span(2,2);
            if s1<p2; s1 = p2; L = s1-s2; elseif s2<p2; s2 = p1; L = s1-s2; else L = Az_shade(2); end
            Az_shade(2) = L;
        end
        if Az_shade(2)>0 && Az_shade(3)>0
            p1 = Az_shade_span(1,2);
            p2 = Az_shade_span(2,2);
            s1 = Az_shade_span(1,3);
            s2 = Az_shade_span(2,3);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(3); end
            Az_shade(3) = L;
        end
        if Az_shade(4)>0 && Az_shade(3)>0
            p1 = Az_shade_span(1,4);
            p2 = Az_shade_span(2,4);
            s1 = Az_shade_span(1,3);
            s2 = Az_shade_span(2,3);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(3); end
            Az_shade(3) = L;
        end
        if Az_shade(4)>0 && Az_shade(5)>0
            p1 = Az_shade_span(1,4);
            p2 = Az_shade_span(2,4);
            s1 = Az_shade_span(1,5);
            s2 = Az_shade_span(2,5);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(5); end
            Az_shade(5) = L;
        end
        if Az_shade(6)>0 && Az_shade(5)>0
            p1 = Az_shade_span(1,6);
            p2 = Az_shade_span(2,6);
            s1 = Az_shade_span(1,5);
            s2 = Az_shade_span(2,5);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(5); end
            Az_shade(5) = L;
        end
        if Az_shade(7)>0 && Az_shade(6)>0
            p1 = Az_shade_span(1,7);
            p2 = Az_shade_span(2,7);
            s1 = Az_shade_span(1,6);
            s2 = Az_shade_span(2,6);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(6); end
            Az_shade(6) = L;
        end
        if Az_shade(7)>0 && Az_shade(8)>0
            p1 = Az_shade_span(1,7);
            p2 = Az_shade_span(2,7);
            s1 = Az_shade_span(1,8);
            s2 = Az_shade_span(2,8);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(8); end
            Az_shade(8) = L;
        end
        if Az_shade(8)>0 && Az_shade(9)>0
            p1 = Az_shade_span(1,8);
            p2 = Az_shade_span(2,8);
            s1 = Az_shade_span(1,9);
            s2 = Az_shade_span(2,9);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(9); end
            Az_shade(9) = L;
        end
        if Az_shade(10)>0 && Az_shade(9)>0
            p1 = Az_shade_span(1,10);
            p2 = Az_shade_span(2,10);
            s1 = Az_shade_span(1,9);
            s2 = Az_shade_span(2,9);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(9); end
            Az_shade(9) = L;
        end
        if Az_shade(10)>0 && Az_shade(11)>0
            p1 = Az_shade_span(1,10);
            p2 = Az_shade_span(2,10);
            s1 = Az_shade_span(1,11);
            s2 = Az_shade_span(2,11);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(11); end
            Az_shade(11) = L;
        end
        if Az_shade(12)>0 && Az_shade(11)>0
            p1 = Az_shade_span(1,12);
            p2 = Az_shade_span(2,12);
            s1 = Az_shade_span(1,11);
            s2 = Az_shade_span(2,11);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(11); end
            Az_shade(11) = L;
        end
        if Az_shade(13)>0 && Az_shade(12)>0
            p1 = Az_shade_span(1,13);
            p2 = Az_shade_span(2,13);
            s1 = Az_shade_span(1,12);
            s2 = Az_shade_span(2,12);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(12); end
            Az_shade(12) = L;
        end
        
        % 2nd order shade overlap:
        if Az_shade(1)>0 && Az_shade(3)>0
            p1 = Az_shade_span(1,1);
            p2 = Az_shade_span(2,1);
            s1 = Az_shade_span(1,3);
            s2 = Az_shade_span(2,3);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(3); end
            Az_shade(3) = L;
        end
        if Az_shade(4)>0 && Az_shade(2)>0
            p1 = Az_shade_span(1,4);
            p2 = Az_shade_span(2,4);
            s1 = Az_shade_span(1,2);
            s2 = Az_shade_span(2,2);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(2); end
            Az_shade(2) = L;
        end
        if Az_shade(4)>0 && Az_shade(6)>0
            p1 = Az_shade_span(1,4);
            p2 = Az_shade_span(2,4);
            s1 = Az_shade_span(1,6);
            s2 = Az_shade_span(2,6);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(6); end
            Az_shade(6) = L;
        end
        if Az_shade(7)>0 && Az_shade(5)>0
            p1 = Az_shade_span(1,7);
            p2 = Az_shade_span(2,7);
            s1 = Az_shade_span(1,5);
            s2 = Az_shade_span(2,5);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(5); end
            Az_shade(5) = L;
        end
        if Az_shade(7)>0 && Az_shade(9)>0
            p1 = Az_shade_span(1,7);
            p2 = Az_shade_span(2,7);
            s1 = Az_shade_span(1,9);
            s2 = Az_shade_span(2,9);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(9); end
            Az_shade(9) = L;
        end
        if Az_shade(10)>0 && Az_shade(8)>0
            p1 = Az_shade_span(1,10);
            p2 = Az_shade_span(2,10);
            s1 = Az_shade_span(1,8);
            s2 = Az_shade_span(2,8);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(8); end
            Az_shade(8) = L;
        end
        if Az_shade(10)>0 && Az_shade(12)>0
            p1 = Az_shade_span(1,10);
            p2 = Az_shade_span(2,10);
            s1 = Az_shade_span(1,12);
            s2 = Az_shade_span(2,12);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(12); end
            Az_shade(12) = L;
        end
        if Az_shade(13)>0 && Az_shade(11)>0
            p1 = Az_shade_span(1,13);
            p2 = Az_shade_span(2,13);
            s1 = Az_shade_span(1,11);
            s2 = Az_shade_span(2,11);
            if s1>p2; s1 = p2; L = s1-s2; elseif s2>p2; s2 = p1; L = s1-s2; else L = Az_shade(11); end
            Az_shade(11) = L;
        end
        
        Az_shade(Az_shade<0) = 0;
        
        % avarage shading with respect to exterior rows:
          if (azimuth<90 && azimuth>=0) || (azimuth<270 && azimuth>=180)
              A_outer_rim = (Nns-1).*(Az_shade(1).*El_shade(1))+(New-1).*(Az_shade(7).*El_shade(7));
              A_secondary_rim = (Nns-2).*(Az_shade(4).*El_shade(4))+(New-2).*(Az_shade(4).*El_shade(4));
          else
              A_outer_rim = (Nns-1).*(Az_shade(13).*El_shade(13))+(New-1).*(Az_shade(7).*El_shade(7));
              A_secondary_rim = (Nns-2).*(Az_shade(10).*El_shade(10))+(New-2).*(Az_shade(10).*El_shade(10));
          end
        
        S_shade = sum(Az_shade.*El_shade); 
        A_inner_field_shade = S_shade.*((Nns-2).*(New-2));
        A_shade = abs(A_outer_rim+A_secondary_rim+A_inner_field_shade);
        
        shading_percent = 1-A_shade./(S.*Nns.*New);

else
    
     shading_percent =0;
     
end

end

function Feq = structural(azimuth,elevation,wind_velocity,wind_direction)
%%  structural calculates the optical efficiency of HelioFocus's dish by wind and statical deformation
%++++++++++++++++    parameters   ++++++++++++++++++++++++++++++++++++++++++++++++++++
% azimuth - sun's Azimuth angle [Deg]
% elevation - sun's Elevation angle [Deg]
% wind_velocity - from getTMY => TMY(:,3) [m./sec]
% wind_direction - from getTMY => TMY(:,4) [Deg.Azimuth]
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%              the function: 

wel=elevation;
waz=abs(azimuth-wind_direction);

while waz>180
    waz=abs(180-waz);   
end
    
    Az = [0.*ones(1,4),15.*ones(1,4),30.*ones(1,4),45.*ones(1,4),60.*ones(1,4),75.*ones(1,4),90.*ones(1,4),105.*ones(1,4),120.*ones(1,4),135.*ones(1,4),150.*ones(1,4),165.*ones(1,4),180.*ones(1,4)]';
    El = [0,10,45,90]; El = [El,El,El,El,El,El,El,El,El,El,El,El,El]';

    w0=0.01.*[96.99,97.22,97.54,98.05]';

    w7=0.01.*[96.88,96.91,97.01,97.99...
    96.86,96.94,97.05,98.00...
    96.85,96.91,97.09,98.00...
    96.83,96.99,97.16,98.02...
    96.83,97.03,97.24,98.03...
    96.81,97.10,97.38,98.08...
    96.79,97.15,97.49,98.01...
    96.81,97.19,97.57,97.92...
    96.92,97.19,97.57,97.88...
    96.88,97.21,97.57,97.79...
    96.87,97.33,97.59,97.83...
    96.87,97.34,97.53,97.74...
    96.87,97.33,97.59,97.80]';

    w14=0.01.*[96.17,95.47,93.10,95.84...
    96.13,95.52,93.47,95.93...
    96.07,95.34,93.85,96.10...
    96.01,95.73,94.51,96.33...
    96.12,95.78,95.34,96.94...
    96.09,96.21,96.29,97.49...
    95.99,96.53,96.92,97.37...
    96.33,96.76,97.19,97.02...
    96.35,96.76,97.03,96.65...
    96.26,97.03,96.99,96.21...
    96.21,97.20,97.06,95.75...
    96.18,97.18,97.02,95.56...
    96.20,97.14,97.06,95.39]';

    F00 = @(x) interp1([0,10,45,90],w0,x);%0.000000037257496.*x.^3-0.000005128527337.*x.^2+0.000277559523809.*x+0.969900000000001;
    F07 = TriScatteredInterp(El,Az,w7,'natural');
    F14 = TriScatteredInterp(El,Az,w14,'natural');


    if wind_velocity==0    
        Feq = F00(wel);
    end

    if  (wind_velocity>0) && (wind_velocity<=7)
        q = wind_velocity./7;
        Feq = (1-q).*F00(wel) + q.*F07(wel,waz);
    end

    if  (wind_velocity>7) && (wind_velocity<=14)
        q = (wind_velocity-7)./7;
        Feq = (1-q).*F07(wel,waz) + q.*F14(wel,waz);
    end

    if  wind_velocity>14
        q = wind_velocity./14;
        Feq = F14(wel,waz).*(2-q);
    end

    Feq(isnan(Feq)) = 0;
    Feq(Feq<0) = 0;
    
end

function effdust = dust(time_length,min_dirt_eff,t_to_cleaning,t_to_primary_cleaning,reflectivity)
%%  dirt calculates the reflectance of the HelioFocus dish under dirt cover effect
%++++++++++++++++    parameters   ++++++++++++++++++++++++++++++++++++++++++++++++++++
% hours_vec - vector of timing - hourly intervaled - i.e. TMY(:,1) (i)[hour]
% min_dirt_eff - the minimal reflective efficiency - at max dirt
% t_to_cleaning - time periode betwine dish cleaning [hours]
% t_to_primary_cleaning - time periode betwine mechanical dish cleaning [hours]
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%                  the function:
dirt_time=0;dirt_time_prime=0;dirt_effect=1;
%reflectivity = 0.943; % reflectance coeeficcient (reflectancy)
dirt_increment=((1-min_dirt_eff)./2)./t_to_cleaning;
prime_dirt_increment=((1-min_dirt_eff)./2)./t_to_primary_cleaning;

    for i=1:time_length

     if dirt_time>t_to_cleaning
         dirt_time=0;
     else
         dirt_time=dirt_time+1;
     end

     if dirt_time_prime>t_to_primary_cleaning
         dirt_time_prime=0;
     else
         dirt_time_prime=dirt_time_prime+1;
     end

       dirt_effect=1-dirt_increment.*dirt_time-prime_dirt_increment.*dirt_time_prime;


    effdust = dirt_effect.*reflectivity;
    
    end
    
end
