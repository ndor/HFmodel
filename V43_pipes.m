function [] = V43_pipes()
%V21 Summary of this function goes here
%   Detailed explanation goes here
global  DATA...
    INPUT...
    field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    HTS_type...
    file_name...
    stop_flag
    

message = ['Optomizing...'];
refreshdata
drawnow expose update
optimization = 1;

FIELD.OPT_USER = INPUT.OPT_USER;

TMY_file_name = INPUT.TMY_file_name;

% file_name = INPUT.file_name;

% Site Parameters:
% site_name = INPUT.site_name;
% TMY_file = INPUT.TMY_file;
% dUTC = 0; % delta UTC (i.e. UTC1)
% longitude = 105.696716;%34.77499; % (decimal)
latitude = INPUT.latitude;%31.024694; % (decimal)
altitude = INPUT.altitude;%1081.9; % meters over (or under) sea level [m]

% Field Parameters:
field_azimuth_or = INPUT.field_azimuth_or; % rectangular field shape azimuthial angle tilt from north [Deg]
field_elevation_or = INPUT.field_elevation_or; % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
N_dishes_per_cluster = INPUT.N_dishes_per_cluster; % number of dishes at each cluster (minimum 1)
N_clusters_in_field = INPUT.N_clusters_in_field; % number of clusters in the field (minimum 1)
lns = INPUT.lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
lew = INPUT.lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
alpha = INPUT.alpha; % field parralelogram angle (from E-W line, +CW) [deg]
% HEX_distance_from_field_center = INPUT.HEX_distance_from_field_center; % HEX distance from center field - HEX linkage [m]

% Dish Parameters:
Dish_effective_Area = INPUT.Dish_effective_Area; % effective reflective area of each dish [m2]
Slope_Error = INPUT.Slope_Error;%sqrt(2.*2.5.^2); % mirror slope error (production quality) [mili.radians]
reflectivity = INPUT.reflectivity; % mirror max reflectivity [.%]

% Receiver Parameters:
FIELD.RECEIVER.concentration_ratio = Dish_effective_Area/(pi*(0.275.^2));
FIELD.RECEIVER.receiver_peak_efficiency = INPUT.receiver_peak_efficiency;

% HTS Parameters:
        if HTS_type>1
            hot.duct_R = INPUT.HTS_hot.duct_R;
            cold.duct_R = INPUT.HTS_cold.duct_R;
            
            hot.insulation_occupation(1) = INPUT.HTS_hot.insulation_occupation(1);
            hot.insulation_occupation(2) = INPUT.HTS_hot.insulation_occupation(2);
            hot.insulation_occupation(3) = INPUT.HTS_hot.insulation_occupation(3);
%             hLm = INPUT.HTS_hot.insulation_material;
%             hLq = INPUT.HTS_hot.insulation_quality;

            cold.insulation_thickness(1) = INPUT.HTS_cold.insulation_thickness(1);
            cold.insulation_thickness(2) = INPUT.HTS_cold.insulation_thickness(2);
            cold.insulation_thickness(3) = INPUT.HTS_cold.insulation_thickness(3);
%             cLm = INPUT.HTS_cold.insulation_material;
%             cLq = INPUT.HTS_cold.insulation_quality;
        else
            hot.duct_R = INPUT.HTS_hot.duct_R;
            cold.duct_R = INPUT.HTS_cold.duct_R;
            
            hot.insulation_thickness(1) = INPUT.HTS_hot.insulation_thickness(1);
            hot.insulation_thickness(2) = INPUT.HTS_hot.insulation_thickness(2);
            hot.insulation_thickness(3) = INPUT.HTS_hot.insulation_thickness(3);
%             hLm = INPUT.HTS_hot.insulation_material;
%             hLq = INPUT.HTS_hot.insulation_quality;

            cold.insulation_thickness(1) = INPUT.HTS_cold.insulation_thickness(1);
            cold.insulation_thickness(2) = INPUT.HTS_cold.insulation_thickness(2);
            cold.insulation_thickness(3) =INPUT.HTS_cold.insulation_thickness(3);
%             cLm = INPUT.HTS_cold.insulation_material;
%             cLq = INPUT.HTS_cold.insulation_quality;
        end

% HEX Parameters:
FIELD.HEX.HP.approach = INPUT.HP_approach;
FIELD.HEX.HP.pinch = INPUT.HP_pinch;
FIELD.HEX.HP.pressure = INPUT.HP_pressure;
FIELD.HEX.HP.temperature = INPUT.HP_temperature;
FIELD.HEX.MP.approach = INPUT.MP_approach;
FIELD.HEX.MP.pinch = INPUT.MP_pinch;
FIELD.HEX.MP.pressure = INPUT.MP_pressure;
FIELD.HEX.MP.temperature = INPUT.MP_temperature;
FIELD.HEX.LP.approach = INPUT.LP_approach;
FIELD.HEX.LP.pinch = INPUT.LP_pinch;
FIELD.HEX.LP.pressure = INPUT.LP_pressure;
FIELD.HEX.LP.temperature = INPUT.LP_temperature;
FIELD.HEX.WS.pressure = INPUT.WS_pressure;
FIELD.HEX.WS.temperature = INPUT.WS_temperature;
FIELD.HEX.power_block_peak = INPUT.power_block_peak;
FIELD.HEX.design_inlet_temperature = INPUT.design_inlet_temperature;
FIELD.HEX.plant_type = INPUT.plant_type;

%============================%
%       V3 - Model *V3.3
%============================%
% transient_state = 0;
%% INPUT PARAMETERS	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MODEL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

site_DATA(TMY_file_name,latitude,altitude,optimization);
%         DATA.time
%         DATA.Azimuth
%         DATA.Elevation
%         DATA.DNI
%         DATA.wind_velocity
%         DATA.wind_direction
%         DATA.ambient_temp
%         DATA.humidity
%         DATA.step_duration

%% Optimization Border
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
field_azimuth_or_0 = field_azimuth_or;
field_elevation_or_0 = field_elevation_or;
N_dishes_per_cluster_0 = N_dishes_per_cluster; % number of dishes at each cluster (minimum 1)
N_clusters_in_field_0 = N_clusters_in_field; % number of clusters in the field (minimum 1)
lns_0 = lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
lew_0 = lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
alpha_0=alpha; % field parralelogram angle (from E-W line, +CW) [deg]
% HEX_distance_from_field_center_0 = 0; % HEX distance from center field - HEX linkage [m]
            hot_duct_R_0 = INPUT.HTS_hot.duct_R;
            cold_duct_R_0 = INPUT.HTS_cold.duct_R;
            hot_insulation_thickness_0(1) = INPUT.HTS_hot.insulation_thickness(1);
            hot_insulation_thickness_0(2) = INPUT.HTS_hot.insulation_thickness(2);
            hot_insulation_thickness_0(3) = INPUT.HTS_hot.insulation_thickness(3);
            cold_insulation_thickness_0(1) = INPUT.HTS_cold.insulation_thickness(1);
            cold_insulation_thickness_0(2) = INPUT.HTS_cold.insulation_thickness(2);
            cold_insulation_thickness_0(3) = INPUT.HTS_cold.insulation_thickness(3);


X0 = [field_azimuth_or_0;field_elevation_or_0;N_dishes_per_cluster_0;N_clusters_in_field_0;lns_0;lew_0;alpha_0;hot_duct_R_0;cold_duct_R_0;...
    hot_insulation_thickness_0(1);hot_insulation_thickness_0(2);hot_insulation_thickness_0(3);...
    cold_insulation_thickness_0(1);cold_insulation_thickness_0(2);cold_insulation_thickness_0(3)];

%% Optimization Function
OPT = @ (X) field_opt_pipes(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD); %X,DATA,Slope_Error,reflectivity,system_pressure,FIELD
nonlcon = @ (X) nonlinear_constraints(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD);

% A+B<=1
A = []; % [x1,x2,....xn]
b = []; % [max]

Aeq = [];%[0,0,1,0,0,0,0,0,0,0,0,0,0;...
           %0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [];%[;];
ub = [90,10,30,30,150,250,45,0.1,1,1,1,1,1,1,1]';
lb = [-90,-10,1,1,39,sqrt(39.^2-(39./2).^2),0,0.03,0.1,0.03,0,0,0.03,0,0]';

%% Constraints
%     if INPUT.OPT_USER.ublb == 1
        
        field_azimuth_or__min = INPUT.LB.field_azimuth_or;
        field_elevation_or__min = INPUT.LB.field_elevation_or;
        N_dishes_per_cluster__min = INPUT.LB.N_dishes_per_cluster; % number of dishes at each cluster (minimum 1)
        N_clusters_in_field__min = INPUT.LB.N_clusters_in_field; % number of clusters in the field (minimum 1)
        lns__min = INPUT.LB.lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
        lew__min =INPUT.LB.lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
        alpha_min = INPUT.LB.alpha; % field parralelogram angle (from E-W line, +CW) [deg]
        % HEX_distance_from_field_center__min = 0; % HEX distance from center field - HEX linkage [m]
            hot_duct_R_min = INPUT.LB.HTS_hot.duct_R;
            cold_duct_R_min = INPUT.LB.HTS_cold.duct_R;
            hot_insulation_thickness_min(1) = INPUT.LB.HTS_hot.insulation_thickness(1);
            hot_insulation_thickness_min(2) = INPUT.LB.HTS_hot.insulation_thickness(2);
            hot_insulation_thickness_min(3) = INPUT.LB.HTS_hot.insulation_thickness(3);
            cold_insulation_thickness_min(1) = INPUT.LB.HTS_cold.insulation_thickness(1);
            cold_insulation_thickness_min(2) = INPUT.LB.HTS_cold.insulation_thickness(2);
            cold_insulation_thickness_min(3) = INPUT.LB.HTS_cold.insulation_thickness(3);

        field_azimuth_or__max = INPUT.UB.field_azimuth_or;
        field_elevation_or__max = INPUT.UB.field_elevation_or;
        N_dishes_per_cluster__max = INPUT.UB.N_dishes_per_cluster; % number of dishes at each cluster (minimum 1)
        N_clusters_in_field__max = INPUT.UB.N_clusters_in_field; % number of clusters in the field (minimum 1)
        lns__max = INPUT.UB.lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
        lew__max = INPUT.UB.lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
        alpha_max = INPUT.UB.alpha; % field parralelogram angle (from E-W line, +CW) [deg]
        % HEX_distance_from_field_center__max = 0; % HEX distance from center field - HEX linkage [m]
            hot_duct_R_max = INPUT.UB.HTS_hot.duct_R;
            cold_duct_R_max = INPUT.UB.HTS_cold.duct_R;
            hot_insulation_thickness_max(1) = INPUT.UB.HTS_hot.insulation_thickness(1);
            hot_insulation_thickness_max(2) = INPUT.UB.HTS_hot.insulation_thickness(2);
            hot_insulation_thickness_max(3) = INPUT.UB.HTS_hot.insulation_thickness(3);
            cold_insulation_thickness_max(1) = INPUT.UB.HTS_cold.insulation_thickness(1);
            cold_insulation_thickness_max(2) = INPUT.UB.HTS_cold.insulation_thickness(2);
            cold_insulation_thickness_max(3) = INPUT.UB.HTS_cold.insulation_thickness(3);  
        
        X_min = [field_azimuth_or__min,field_elevation_or__min,N_dishes_per_cluster__min,N_clusters_in_field__min,lns__min,lew__min,alpha_min,...
            hot_duct_R_min,cold_duct_R_min,hot_insulation_thickness_min(1),hot_insulation_thickness_min(2),...
            hot_insulation_thickness_min(3),cold_insulation_thickness_min(1),cold_insulation_thickness_min(2),...
            cold_insulation_thickness_min(3)]';
        X_min(lb>X_min) = lb(lb>X_min)
                X_max = [field_azimuth_or__max,field_elevation_or__max,N_dishes_per_cluster__max,N_clusters_in_field__max,lns__max,lew__max,alpha_max,...
            hot_duct_R_max,cold_duct_R_max,hot_insulation_thickness_max(1),hot_insulation_thickness_max(2),...
            hot_insulation_thickness_max(3),cold_insulation_thickness_max(1),cold_insulation_thickness_max(2),...
            cold_insulation_thickness_max(3)]';
        X_max(ub<X_max) = ub(ub<X_max)
        
        ub = X_max;
        lb = X_min;
        
        k = 1;
        for i = 1:length(ub)
            if ub(i)==lb(i)
                Aeq(k,:) = zeros(1,length(ub));
                Aeq(k,i) = 1;
                beq(k,1) = ub(i);
                k = k+1;
            end
        end
        
%     end
    
    Aeq
    beq

if stop_flag==1;return;end
options = psoptimset('UseParallel','always','TolX',0.0005,'TolMesh',0.0005,'TolFun',0.0005,'MeshAccelerator','on','Display','iter','CompletePoll','on','UseParallel','always','ScaleMesh','on','PlotFcns',{@psplotfuncount,@psplotbestx,@psplotbestf,@psplotmeshsize});%
[x,fval,exitflag,output] = patternsearch(OPT,X0,A,b,Aeq,beq,lb,ub,nonlcon,options);

%% OUTPUT is now INPUT for a model run:

x(1) = ceil(x(1));
x(2) = ceil(x(2));
INPUT.field_azimuth_or = x(1);
INPUT.field_elevation_or = x(2);
INPUT.N_dishes_per_cluster = ceil(x(3)); % number of dishes at each cluster (minimum 1)
INPUT.N_clusters_in_field = ceil(x(4)); % number of clusters in the field (minimum 1)
INPUT.lns = x(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
INPUT.lew =x(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
INPUT.alpha = x(7); % field parralelogram angle (from E-W line, +CW) [deg]
% HEX_distance_from_field_center = 0;%X(5); % HEX distance from center field - HEX linkage [m]
            INPUT.HTS_hot.duct_R = x(8);%INPUT.HTS_hot.duct_R;
            INPUT.HTS_cold.duct_R = x(9);%INPUT.HTS_cold.duct_R;
            INPUT.HTS_hot.insulation_thickness(1) = x(10);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_hot.insulation_thickness(2) = x(11);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_hot.insulation_thickness(3) = x(12);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_cold.insulation_thickness(1) = x(13);%INPUT.HTS_cold.insulation_thickness;
            INPUT.HTS_cold.insulation_thickness(2) = x(14);%INPUT.HTS_cold.insulation_thickness;
            INPUT.HTS_cold.insulation_thickness(3) = x(15);%INPUT.HTS_cold.insulation_thickness;

%%
message = ['Optimization Complete'];
refreshdata
drawnow expose update
end

function [Cineq,Ceq] = nonlinear_constraints(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
% Field Variables:

OPT_USER = FIELD.OPT_USER;

    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = X(3);%ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = X(4);%ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]

% Pre_or_Post:
% if Pre_or_Post = 0 : Pre-annual-calculation
% if Pre_or_Post = 1 : Post-annual-calculation

A_max = OPT_USER.A_max;
B_max = OPT_USER.B_max;
C_max = OPT_USER.C_max;
D_max = OPT_USER.D_max;
E_max = OPT_USER.E_max;
F_max = OPT_USER.F_max;
G_max = OPT_USER.G_max;
H_max = OPT_USER.H_max;

A_min = OPT_USER.A_min;
B_min = OPT_USER.B_min;
C_min = OPT_USER.C_min;
D_min = OPT_USER.D_min;
E_min = OPT_USER.E_min;
F_min = OPT_USER.F_min;
G_min = OPT_USER.G_min;
H_min = OPT_USER.H_min;


PIPES = field_structure_pipes(0);
FIELD.PIPES = PIPES;


    Ceq(1) = X(3) - ceil(X(3));
    Ceq(2) = X(4) - ceil(X(4));
    
    if X(4)==1
        j =2;
    else
        Ceq(3) = X(4)./2 - ceil(X(4)./2);
        j=3;
    end
    

% minimal distance between dishes constraint X=39m :
Cineq(1) = -sqrt((X(6)).^2+(X(5)./2).^2)+39;
i = 1;

%% Field size (area)
% Limited area or specific span of area – specified by user
% Maximizes the value: A(var)=  (net production)/Cost

% goal function:    

    if (A_max~=0 && A_min~=0) && (A_max~=A_min)
        
        i = i+1;
        Cineq(i) = field_size(ceil(X(3)),ceil(X(4)),X(5),X(6),X(7)) - A_max; % field_size(N_dishes_per_cluster,N_clusters_in_field,lns,lew)
        Cineq(i+1) = -field_size(ceil(X(3)),ceil(X(4)),X(5),X(6),X(7)) + A_min; 
        
    elseif (A_max~=0 && A_min~=0) && (A_max==A_min)
        
        j = j+1;
        Ceq(j) = field_size(ceil(X(3)),ceil(X(4)),X(5),X(6),X(7)) - A_max; % field_size(N_dishes_per_cluster,N_clusters_in_field,lns,lew)

    end

%% Project cost (CAPEX)
% Limited cost (max cost) – specified by user
% Maximizes the value: B(var)=  (net production)/Cost

% goal function:    

    if (B_max~=0 && B_min~=0) && (B_max~=B_min)
        
        i = i+1;
        Cineq(i) = project_cost(X,system_pressure) - B_max; % project_cost(X,system_pressure,HTS_type)
        Cineq(i+1) = -project_cost(X,system_pressure) + B_min;
        
    elseif (B_max~=0 && B_min~=0) && (B_max==B_min)
        
        j = j+1;
        Ceq(j) = project_cost(X,system_pressure) - B_max; % project_cost(X,system_pressure,HTS_type)
        
    end

%% Annual Energy (Enet)
% Minimal Enet – specified by user
% Maximizes the value: C(var)=  (net production)/Cost

% goal function:    
    
    if (C_max~=0 && C_min~=0) && (C_max~=C_min)

        i = i+1;
        Cineq(i) = net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - C_max; % net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
        Cineq(i+1) = net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) + C_min;
            
    elseif (C_max~=0 && C_min~=0) && (C_max==C_min)
        
        j = j+1;        
        Ceq(j) = net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - C_max; % net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)

    end

    
%% Installed Capacity in KW installed
% Limited span of Pnet – specified by user
% Maximizes the value: D(var)=  (net production)/Cost

    if (D_max~=0 && D_min~=0) && (D_max~=D_min)
        
        i = i+1;
        Cineq(i) = capacity(X,Slope_Error,reflectivity,FIELD) - D_max; % capacity(X,Slope_Error,reflectivity,FIELD)
        Cineq(i+1) = -capacity(X,Slope_Error,reflectivity,FIELD) + D_min;
    
    elseif (D_max~=0 && D_min~=0) && (D_max==D_min)
        
        j = j+1;         
        Ceq(j) = capacity(X,Slope_Error,reflectivity,FIELD) - D_max; % capacity(X,Slope_Error,reflectivity,FIELD)

    end


%% LCOE (required LCOE, meaning targeted cost of electricity) in $/KWh
% Maximum or target span LCOE ($/(KW_installed )–) - specified by user
% Maximize the value:  E(var)=   (net production)/Cost

% goal function:  

    if (E_max~=0 && E_min~=0) && (E_max~=E_min)
        
        i = i+1;
        Cineq(i) = lcoe(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - E_max; % lcoe(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
        Cineq(i+1) = -lcoe(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) + E_min;
        
    elseif (E_max~=0 && E_min~=0) && (E_max==E_min)
        
        j = j+1;  
        Ceq(j) = lcoe(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - E_max; % lcoe(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)


    end


%% Number of dishes (a specified number of dishes)
% Maximum or span of dishes in field - specified by user
% Maximize the value:  F(var)=   (net production)/Cost

% goal function:    

    if (F_max~=0 && F_min~=0) && (F_max~=F_min)
        
        i = i+1;   
        Cineq(i) = X(3).*X(4) - F_max;
        Cineq(i+1) = -X(3).*X(4) + F_min;
        
    elseif (F_max~=0 && F_min~=0) && (F_max==F_min)
        
        j = j+1;  
        Ceq(j) = X(3).*X(4) - F_max;

    end


%% Steam output in t/h pressure and temperature?
% Minimum or span of steam output, annual or design point - specified by user
% Maximize the value:  G(var)=   (net production)/Cost

% goal function:    

    if (G_max~=0 && G_min~=0) && (G_max~=G_min)
        
        i = i+1;
        Cineq(i) = steam_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - G_max; % steam_annual(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
        Cineq(i+1) = -steam_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) + G_min;
        
    elseif (G_max~=0 && G_min~=0) && (G_max==G_min)
        
        j = j+1;  
        Ceq(j) = steam_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - G_max; % steam_annual(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)


    end


%% Consumption load (the required consumption curve given by a customer. Can be with time of day and price he gets for electricity)
% Maximum acceptable consumption[KWe]/(Gross production [KWe] )  ?[%] (0<x<1) - specified by user
% Maximize the value:  H(var)=   (net production)/Cost

% goal function:    

    if (H_max~=0 && H_min~=0) && (H_max~=H_min)

        i = i+1;
        
        Cineq(i) = consumption_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - G_max; % consumption_annual(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
        Cineq(i+1) = -consumption_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) + G_min;
        
    elseif (H_max~=0 && H_min~=0) && (H_max==H_min)
        
        j = j+1;  
        Ceq(j) = consumption_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD) - G_max; % consumption_annual(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)

    end

    
end

function area = field_size(N_dishes_per_cluster,N_clusters_in_field,lns,lew,alpha)
%% for optimization:
N_dishes_per_cluster = ceil(N_dishes_per_cluster);
N_clusters_in_field = ceil(N_clusters_in_field);
if N_clusters_in_field>1 && N_clusters_in_field./2~=ceil(N_clusters_in_field./2)
    N_clusters_in_field = 2.*ceil(N_clusters_in_field./2);
end
%%

% field size:
n = N_clusters_in_field./2;
n(n<1) = 1; n(n>=1) = 2;
Xns = n.*(N_dishes_per_cluster./2).*lns+2.*(11+4);
Xew = ((N_clusters_in_field./n)).*lew+2.*(11+4);
xns = Xew.*tand(alpha);
area = Xns.*Xew - Xew.*xns;

end

function cost = project_cost(X,system_pressure)

    transient_state = 0;
    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]
    HEX_distance_from_field_center = 0;%X(5); % HEX distance from center field - HEX linkage [m]

PIPES = field_structure_pipes(0);

cost = PIPES.project_cost;
  
end

function E_net_annual = net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD)
%FIELD_OPT Summary of this function goes here
%   Detailed explanation goes here

    transient_state = 0;
    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]
    HEX_distance_from_field_center = 50;%X(5); % HEX distance from center field - HEX linkage [m]

    PIPES = field_structure_pipes(0);

            FIELD.PIPES = PIPES;

                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 time_length = size(DISH.dish_power_pr_sqr_m,1);  

                 % Model calculation:

                    for i=1:time_length

                        FIELD = steady_field_pipes(DISH.dish_power_pr_sqr_m(i),FIELD,DATA.wind_velocity(i),DATA.ambient_temp(i),DATA.humidity(i));
                        Steam(i) = (DATA.DNI_commonality(i).*(FIELD.POWER.mdot_HP+FIELD.POWER.mdot_MP+FIELD.POWER.mdot_LP));
                        Consumption(i) = DATA.DNI_commonality(i).*(PIPES.dishes_in_field.*DISH.consumption(i)+FIELD.PIPES.all_blowers_consumption+FIELD.POWER.consumption.all_pumps_consumption); 
                        NET_Electric(i) = (DATA.DNI_commonality(i).*FIELD.POWER.Electric - Consumption(i));

                    end
                    
                    Steam(isnan(Steam)==1) = 0;
                    Consumption(isnan(Consumption)==1) = 0;
                    NET_Electric(isnan(NET_Electric)==1) = 0;

                    mass_steam_annual = (sum(Steam)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Steam [t/h] production
                    Consumption_annual = (sum(Consumption)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Consumption [KWh] production
                    E_net_annual = (sum(NET_Electric)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Nete [KWh] production
                    
end

function installed_capacity = capacity(X,Slope_Error,reflectivity,FIELD)

    transient_state = 0;
    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]
    HEX_distance_from_field_center = 0;%X(5); % HEX distance from center field - HEX linkage [m]

D0 = 0.5.*0.9.*0.9.*450.*950.*N_dishes_per_cluster.*N_clusters_in_field; % maximal fantastic installed

% goal function:    
        % design_DATA.time
        design_DATA.Azimuth = 180;
        design_DATA.Elevation = 90;
        design_DATA.DNI = 950;
        design_DATA.wind_velocity = 0;
        design_DATA.wind_direction = 180;
        design_DATA.ambient_temp = 20;
        design_DATA.humidity = 30;
        % design_DATA.step_duration

design_DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,design_DATA);
design_FIELD = steady_field_pipes(design_DISH.dish_power_pr_sqr_m,FIELD,design_DATA.wind_velocity,design_DATA.ambient_temp,design_DATA.humidity);
design_E_net = design_FIELD.POWER.Electric - (FIELD.PIPES.dishes_in_field.*design_DISH.consumption+design_FIELD.PIPES.all_blowers_consumption+design_FIELD.POWER.consumption.all_pumps_consumption);

installed_capacity = design_E_net;

end

function consumption_load = consumption_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD)
%FIELD_OPT Summary of this function goes here
%   Detailed explanation goes here

    transient_state = 0;
    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]
    HEX_distance_from_field_center = 50;%X(5); % HEX distance from center field - HEX linkage [m]
    hot_r = X(8); % hot on dish pipe radius [m]
    hot_insulation1 = X(9); % hot inner insulation1 (MicroThern) thickness [m];
    hot_insulation2 = X(10); % hot inner insulation2 thickness [m];
    cold_r = X(11); % cold on dish pipe radius [m]
    cold_insulation1 = X(12); % cold inner insulation1 thickness [m];
    cold_insulation2 = X(13); % cold inner insulation2 thickness [m];


            PIPES = field_structure_pipes(0);
            FIELD.PIPES = PIPES;

                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 time_length = size(DISH.dish_power_pr_sqr_m,1);  

                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 % Model calculation:

                    for i=1:time_length

                        FIELD = steady_field_pipes(DISH.dish_power_pr_sqr_m(i),FIELD,DATA.wind_velocity(i),DATA.ambient_temp(i),DATA.humidity(i));
                        Steam(i) = (DATA.DNI_commonality(i).*(FIELD.POWER.mdot_HP+FIELD.POWER.mdot_MP+FIELD.POWER.mdot_LP));
                        Consumption(i) = DATA.DNI_commonality(i).*(PIPES.dishes_in_field.*DISH.consumption(i)+FIELD.PIPES.all_blowers_consumption+FIELD.POWER.consumption.all_pumps_consumption); 
                        NET_Electric(i) = (DATA.DNI_commonality(i).*FIELD.POWER.Electric - Consumption(i));

                    end
                    
                    Steam(isnan(Steam)==1) = 0;
                    Consumption(isnan(Consumption)==1) = 0;
                    NET_Electric(isnan(NET_Electric)==1) = 0;

                    mass_steam_annual = (sum(Steam)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Steam [t/h] production
                    Consumption_annual = (sum(Consumption)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Consumption [KWh] production
                    E_net_annual = (sum(NET_Electric)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Nete [KWh] production
                    
                    consumption_load = Consumption_annual./(E_net_annual+Consumption_annual);
                    
end

function mass_steam_annual = steam_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD)
%FIELD_OPT Summary of this function goes here
%   Detailed explanation goes here

    transient_state = 0;
    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]
    HEX_distance_from_field_center = 50;%X(5); % HEX distance from center field - HEX linkage [m]

            PIPES = field_structure_pipes(0);
            FIELD.PIPES = PIPES;

                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 time_length = size(DISH.dish_power_pr_sqr_m,1);  

                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 % Model calculation:

                    for i=1:time_length

                        FIELD = steady_field_pipes(DISH.dish_power_pr_sqr_m(i),FIELD,DATA.wind_velocity(i),DATA.ambient_temp(i),DATA.humidity(i));
                        Steam(i) = (DATA.DNI_commonality(i).*(FIELD.POWER.mdot_HP+FIELD.POWER.mdot_MP+FIELD.POWER.mdot_LP));
                        Consumption(i) = DATA.DNI_commonality(i).*(PIPES.dishes_in_field.*DISH.consumption(i)+FIELD.PIPES.all_blowers_consumption+FIELD.POWER.consumption.all_pumps_consumption); 
                        NET_Electric(i) = (DATA.DNI_commonality(i).*FIELD.POWER.Electric - Consumption(i));

                    end
                    
                    Steam(isnan(Steam)==1) = 0;
                    Consumption(isnan(Consumption)==1) = 0;
                    NET_Electric(isnan(NET_Electric)==1) = 0;

                    mass_steam_annual = (sum(Steam)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Steam [t/h] production
                    Consumption_annual = (sum(Consumption)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Consumption [KWh] production
                    E_net_annual = (sum(NET_Electric)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Nete [KWh] production
                    
end

function LCOE = lcoe(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD)

    E_net_annual = net_annual(X,DATA,Slope_Error,reflectivity,system_pressure,FIELD);
    kd = 0.06; %Real debt interest rate
    n = 20; %Lifetime of the project [years]
    crf = (kd.*(1+kd).^n)./((1+kd).^n-1);
    Capex = project_cost(X,system_pressure,HTS_type);
    Opex = FIELD.PIPES.OnM.*n;
    LCOE = (crf.*Capex+Opex)./E_net_annual;

end