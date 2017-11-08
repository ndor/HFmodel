function [] = V42_pipes()
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


optimization = 0;

TMY_file_name = INPUT.TMY_file_name;

file_name = INPUT.file_name;
[status,msginfo] = xlswrite(file_name,[]);

% Site Parameters:
site_name = INPUT.site_name;
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
HEX_distance_from_field_center = INPUT.HEX_distance_from_field_center; % HEX distance from center field - HEX linkage [m]


% Dish Parameters:
Dish_effective_Area = INPUT.Dish_effective_Area; % effective reflective area of each dish [m2]
Slope_Error = INPUT.Slope_Error;%sqrt(2.*2.5.^2); % mirror slope error (production quality) [mili.radians]
reflectivity = INPUT.reflectivity; % mirror max reflectivity [.%]

% Receiver Parameters:
FIELD.RECEIVER.concentration_ratio = Dish_effective_Area/(pi*(0.275.^2));
FIELD.RECEIVER.receiver_peak_efficiency = INPUT.receiver_peak_efficiency;

% HTS Parameters:

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
%       V2 - Model *V2.2
%============================%
transient_state = 1;

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

PIPES = field_structure_pipes(transient_state);
%         PIPES.pressure
%         PIPES.dishes_per_cluster
%         PIPES.clusters_in_field
%         PIPES.dishes_in_field
%         PIPES.dt
%         PIPES.N
            N = PIPES.N;
            
%             % GLOBALS:
%             G_sections_length(1) = 0;
%             G_sections_tag(1) = {'receiver'};
%                      for j = 1:N
%                         G_sections_length(j+1) = PIPES.SECTION{j}.length+G_sections_length(j); 
%                         G_sections_tag{j+1} = {PIPES.SECTION{j}.tag};
%                      end
         
DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);
%         DISH.Slope_Error_intercept
%         DISH.effdust
%         DISH.shading_percent
%         DISH.effstruct
%         DISH.efficiency
%         DISH.dish_power_pr_sqr_m
%         DISH.consumption
        
 time_length = size(DISH.dish_power_pr_sqr_m,1)       
 N = size(PIPES.SECTION,2);
 FIELD.PIPES = PIPES;
 x = 0;
message = ['Prompting Calculation, Please Wait....'];
status_window=waitbar(x,message);
refreshdata
drawnow expose update

tic      
        
% Model calculation:
for i = 1:time_length
        if stop_flag==1
            delete(status_window);
            break
        end
    
    FIELD = transient_field_pipes(DATA.step_duration,DISH.dish_power_pr_sqr_m(i),FIELD,DATA.wind_velocity(i),DATA.ambient_temp(i),DATA.humidity(i));
%         FIELD.RECEIVER
%                 FIELD.RECEIVER.efficiency
%                 FIELD.RECEIVER.concentration_ratio
%                 FIELD.RECEIVER.receiver_peak_efficiency
%                 FIELD.RECEIVER.mdot
%         FIELD.PIPES
%                 FIELD.PIPES.SECTION{}.tag
%                 FIELD.PIPES.SECTION{}.flow_multiplier
%                 FIELD.PIPES.SECTION{}.geometry
%                 FIELD.PIPES.SECTION{}.elbows
%                 FIELD.PIPES.SECTION{}.length
%                 FIELD.PIPES.SECTION{}.mesh
%                 FIELD.PIPES.SECTION{}.hot_temperature1
%                 FIELD.PIPES.SECTION{}.cold_temperature1
%                 FIELD.PIPES.SECTION{}.hot_temperature2
%                 FIELD.PIPES.SECTION{}.cold_temperature2
%                 FIELD.PIPES.SECTION{}.hot_pressure1
%                 FIELD.PIPES.SECTION{}.cold_pressure1
%                 FIELD.PIPES.SECTION{}.hot_pressure2
%                 FIELD.PIPES.SECTION{}.cold_pressure2
%                 FIELD.PIPES.SECTION{}.stainless_steel_cost
%                 FIELD.PIPES.SECTION{}.carbon_steel_cost
%                 FIELD.PIPES.SECTION{}.insulation1_cost
%                 FIELD.PIPES.SECTION{}.insulation2_cost
%                 FIELD.PIPES.SECTION{}.aluminium_cost
%                 FIELD.PIPES.SECTION{}.section_cost
%                 FIELD.PIPES.SECTION{}.sections_field_cost
%                 FIELD.PIPES.SECTION{}.dP_hot
%                 FIELD.PIPES.SECTION{}.dP_cold
%                 FIELD.PIPES.SECTION{}.dP_total
%             FIELD.PIPES.all_blowers_consumption
%             FIELD.PIPES.blower_consumption
%             FIELD.PIPES.dP_cluster
%             FIELD.PIPES.dP_blowers
%             FIELD.PIPES.dP_RECEIVER
%             FIELD.PIPES.dP_HEX
%             FIELD.PIPES.dt
%             FIELD.PIPES.dishes_in_field
%             FIELD.PIPES.clusters_in_field
%             FIELD.PIPES.dishes_per_cluster
%             FIELD.PIPES.pressure
%         FIELD.HEX
%                 FIELD.HEX.HP
%                     FIELD.HEX.HP.temperature
%                     FIELD.HEX.HP.pressure
%                     FIELD.HEX.HP.pinch
%                     FIELD.HEX.HP.aproach
%                 FIELD.HEX.MP
%                     FIELD.HEX.MP.temperature
%                     FIELD.HEX.MP.pressure
%                     FIELD.HEX.MP.pinch
%                     FIELD.HEX.MP.aproach
%                 FIELD.HEX.LP
%                     FIELD.HEX.LP.temperature
%                     FIELD.HEX.LP.pressure
%                     FIELD.HEX.LP.pinch
%                     FIELD.HEX.LP.aproach
%                 FIELD.HEX.WS
%                     FIELD.HEX.WS.temperature
%                     FIELD.HEX.WS.pressure
%                 FIELD.HEX.power_block_peak
%                 FIELD.HEX.design_inlet_temperature
%         FIELD.N_dishes_in_field
%                 FIELD.POWER.Electrical
%                 FIELD.POWER.PB_efficiency
%                 FIELD.POWER.Thermal_efficiency
%                 FIELD.POWER.Thermal
%                 FIELD.POWER.mdot_LP
%                 FIELD.POWER.mdot_MP
%                 FIELD.POWER.mdot_HP
%                 FIELD.POWER.used_air_energy
%                 FIELD.POWER.T_HEX_outlet
%                 FIELD.POWER.consumption
%                         FIELD.POWER.consumption.HP_pump_consumption
%                         FIELD.POWER.consumption.MP_pump_consumption
%                         FIELD.POWER.consumption.LP_pump_consumption
%                         FIELD.POWER.consumption.All_pumps_consumption

        Tin_receiver(i) = FIELD.PIPES.SECTION{1}.cold_temperature1;
        Tout_receiver(i) = FIELD.PIPES.SECTION{1}.hot_temperature1;
        Tin_HEX(i) = FIELD.PIPES.SECTION{N}.hot_temperature2;
        Tout_HEX(i) = FIELD.PIPES.SECTION{N}.cold_temperature2;
        dP_cluster(i) = FIELD.PIPES.dP_cluster;
        mdot(i) = FIELD.RECEIVER.mdot;
        receiver_efficiency(i) = FIELD.RECEIVER.efficiency;
        dP_receiver(i) = FIELD.PIPES.dP_HEX;
        dP_HEX(i) = FIELD.PIPES.dP_RECEIVER;
        
         for j = 1:N
            Tout_HOT_sections(i,j) = FIELD.PIPES.SECTION{j}.hot_temperature2;
            dTinstruments_HOT_sections(i,j) = FIELD.PIPES.SECTION{j}.hot_temperature_loss_by_components;
            Tout_COLD_sections(i,j) = FIELD.PIPES.SECTION{j}.cold_temperature2;
            dTinstruments_COLD_sections(i,j) = FIELD.PIPES.SECTION{j}.cold_temperature_loss_by_components;
            dP_HOT_sections(i,j) = FIELD.PIPES.SECTION{j}.dP_hot;
            dP_COLD_sections(i,j) = FIELD.PIPES.SECTION{j}.dP_cold;
         end

        dP_blowers(i,1:size(FIELD.PIPES.dP_blowers,2)) = FIELD.PIPES.dP_blowers;
        blowers_consumption(i,1:size(FIELD.PIPES.blowers_consumption,2)) = FIELD.PIPES.blowers_consumption;
        all_blowers_consumption(i) = FIELD.PIPES.all_blowers_consumption;

        mdot_HP(i) = FIELD.POWER.mdot_HP;
        mdot_MP(i) = FIELD.POWER.mdot_MP;
        mdot_LP(i) = FIELD.POWER.mdot_LP;
        mdot_total(i) = (FIELD.POWER.mdot_HP+FIELD.POWER.mdot_MP+FIELD.POWER.mdot_LP);
        Thermal_efficiency(i) = FIELD.POWER.Thermal_efficiency;
        Thermal(i) = FIELD.POWER.Thermal;
        PB_efficiency(i) = FIELD.POWER.PB_efficiency;
        GROSS_Electric(i) = FIELD.POWER.Electric;
        HP_pump_consumption(i) = FIELD.POWER.consumption.HP_pump_consumption;
        MP_pump_consumption(i) = FIELD.POWER.consumption.MP_pump_consumption;
        LP_pump_consumption(i) = FIELD.POWER.consumption.LP_pump_consumption;
        all_pumps_consumption(i) = FIELD.POWER.consumption.all_pumps_consumption;
        NET_Electric(i) = GROSS_Electric(i) - (PIPES.dishes_in_field.*DISH.consumption(i)+all_blowers_consumption(i)+all_pumps_consumption(i));

        singularity_count(i) = FIELD.model_run_singularity; 
        
        t = toc./60;
        DONE = 100*i/time_length;
        disp([num2str(DONE) ' % DONE'])
        x = i/time_length;
        message = ['Calculating Model, (' num2str(t) ' min) ' num2str(DONE) '% DONE...'];
        waitbar(x,status_window,message);
        refreshdata
        drawnow expose update
        
        % GLOBALS:
        G_time = DATA.time(i);
        G_sun_azimuth = DATA.Azimuth(i);
        G_sun_elevation = DATA.Elevation(i);
        G_DNI = DATA.DNI(i);
        G_wind_velocity = DATA.wind_velocity(i);
        G_wind_direction = DATA.wind_direction(i);
        G_Tambient = DATA.ambient_temp(i);
        G_humidity = DATA.humidity(i);
        
        G_dish_efficiency = DISH.efficiency(i);

        G_Thot = [Tout_receiver(i),Tout_HOT_sections(i,:)];
        G_Tcold = [Tin_receiver(i),Tout_COLD_sections(i,:)];
        G_mdot = FIELD.RECEIVER.mdot;
        G_receiver_efficiency = FIELD.RECEIVER.efficiency;
        
        G_HEX_efficiency = FIELD.POWER.Thermal_efficiency;
        G_PB_efficiency = FIELD.POWER.PB_efficiency;
        G_steam_tot = mdot_total(i);
        G_all_blowers_consumption = FIELD.PIPES.all_blowers_consumption;
        G_all_pumps_consumption = FIELD.POWER.consumption.all_pumps_consumption;
        G_all_dishes_consumption = PIPES.dishes_in_field.*DISH.consumption(i);
        G_gross_electric = FIELD.POWER.Electric;
        G_net_electric = NET_Electric(i);
        
end

total_consumption = PIPES.dishes_in_field.*DISH.consumption + all_blowers_consumption' + all_pumps_consumption';
system_efficiency = NET_Electric'./(PIPES.dishes_in_field.*Dish_effective_Area.*DATA.DNI./1000);
system_efficiency(system_efficiency<0) = 0;
system_efficiency(isnan(system_efficiency)==1) = 0;

shading_output_impact_approximation_hourly = NET_Electric(DATA.Elevation>0).*(2-DISH.non_shaded(DATA.Elevation>0)');
annual_shading_output_impact_approximation = ((sum(shading_output_impact_approximation_hourly)-sum(NET_Electric(DATA.Elevation>0)))./sum(NET_Electric(DATA.Elevation>0)));
% annual_shading_output_impact_approximation(isinf(annual_shading_output_impact_approximation) | isnan(annual_shading_output_impact_approximation)) = 0;

%% MODEL OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [file_name '.mat']
save(f);

%% INPUT:

xfile_name = {'File Name: ',file_name};

% Site Parameters:
xsite_name = {'Site name: ',site_name};
% dUTC = 0; % delta UTC (i.e. UTC1)
% longitude = 105.696716;%34.77499; % (decimal)
xlatitude = ['Latitude [Deg]: ',num2cell(latitude)];%31.024694; % (decimal)
xaltitude = ['Altitude [m]: ',num2cell(altitude)];%1081.9; % meters over (or under) sea level [m]

% Field Parameters:
xfield_azimuth_or = ['Field Azimuth Orientation [Deg]; ',num2cell(field_azimuth_or)]; % rectangular field shape azimuthial angle tilt from north [Deg]
xfield_elevation_or = ['Field Elevation Orientation [Deg]; ',num2cell(field_elevation_or)]; % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
xN_dishes_per_cluster = ['Number of Dishes per Cluster: ',num2cell(N_dishes_per_cluster)]; % number of dishes at each cluster (minimum 1)
xN_clusters_in_field = ['Number of Clusters in Field: ',num2cell(N_clusters_in_field)]; % number of clusters in the field (minimum 1)
xN_dishes_in_field = ['Number of Dishes in Field: ',num2cell(PIPES.dishes_in_field)]; % number of dishs in the field (minimum 1)
xlns = ['N-S Dishes Distance [m]: ',num2cell(lns)]; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
xlew = ['E-W Dishes Distance [m]: ',num2cell(lew)]; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
xalpha = ['Alpha Angle [deg]: ',num2cell(alpha)]; % field parralelogram angle (from E-W line, +CW) [deg]
xHEX_distance_from_field_center = ['HEX Eccentric Distance [m]: ',num2cell(HEX_distance_from_field_center)]; % HEX distance from center field - HEX linkage [m]
xsystem_pressure = ['System Pressure [BarA]: ',num2cell(system_pressure)]; % air system pressure [Bar]
xR1 = ['Inner N1 Pipe Radius [m]: ',num2cell(INPUT.HTS_hot.duct_R)]; % N1 pipe minimal radius [m]
xR2 = ['Outer N1 Pipe Radius [m]: ',num2cell(INPUT.HTS_cold.duct_R)]; % pipe proportion = R2/R1;
xHOT_ins1 = [['1st HOT Insulation Layer:' insulation_type_inv(INPUT.HTS_hot.insulation_material(1)) ', ' num2str(INPUT.HTS_hot.insulation_quality(1).*100) ' % Quality Thickness [m]: '],num2cell(INPUT.HTS_hot.insulation_thickness(1))]; % pipe cut area scalar
xHOT_ins2 = [['2nd HOT Insulation Layer:' insulation_type_inv(INPUT.HTS_hot.insulation_material(2)) ', ' num2str(INPUT.HTS_hot.insulation_quality(2).*100) ' % Quality Thickness [m]: '],num2cell(INPUT.HTS_hot.insulation_thickness(2))]; % inner insulation1 (MicroTherm) thickness / occupation of cross flow available [%]
xHOT_ins3 = [['3rd HOT Insulation Layer:' insulation_type_inv(INPUT.HTS_hot.insulation_material(3)) ', ' num2str(INPUT.HTS_hot.insulation_quality(3).*100) ' % Quality Thickness [m]: '],num2cell(INPUT.HTS_hot.insulation_thickness(3))]; % inner insulation2 thickness  / occupation of cross flow available [%]
xCOLD_ins1 = [['1st COLD Insulation Layer:' insulation_type_inv(INPUT.HTS_cold.insulation_material(1)) ', ' num2str(INPUT.HTS_cold.insulation_quality(1).*100) ' % Quality Thickness [m]: '],num2cell(INPUT.HTS_cold.insulation_thickness(1))]; % pipe cut area scalar
xCOLD_ins2 = [['2nd COLD Insulation Layer:' insulation_type_inv(INPUT.HTS_cold.insulation_material(2)) ', ' num2str(INPUT.HTS_cold.insulation_quality(2).*100) ' % Quality Thickness [m]: '],num2cell(INPUT.HTS_cold.insulation_thickness(2))]; % inner insulation1 (MicroTherm) thickness / occupation of cross flow available [%]
xCOLD_ins3 = [['3rd COLD Insulation Layer:' insulation_type_inv(INPUT.HTS_cold.insulation_material(3)) ', ' num2str(INPUT.HTS_cold.insulation_quality(3).*100) ' % Quality Thickness [m]: '],num2cell(INPUT.HTS_cold.insulation_thickness(3))]; % inner insulation2 thickness  / occupation of cross flow available [%]
if HTS_type==2
    HTS = ['After Each Dish'];
elseif HTS_type==3
    HTS = ['After Two Dishes'];
elseif HTS_type==1
    HTS = ['Separate Pipes'];
end
xHTS = {'HTS Type: ',HTS}; % altering between 1 section between dishes (sections_alteration=0) or 1 section between two dishes (sections_alteration=1)

% Dish Parameters:
xDish_effective_Area = ['Dish Effective Area [m2]: ',num2cell(Dish_effective_Area)]; % effective reflective area of each dish [m2]
xSlope_Error = ['Slope Error Deviation (RMS) [m.rad]: ',num2cell(Slope_Error)]; % mirror slope error (production quality) [mili.radians]
xreflectivity = ['Mirror Reflectivity: ',num2cell(reflectivity)]; % mirror max reflectivity [.%]

% Receiver Parameters:
xconcentration_ratio = ['Concentration Ration: ',num2cell(FIELD.RECEIVER.concentration_ratio)] ;
xreceiver_peak_efficiency = ['Receiver Peak Efficiency: ',num2cell(FIELD.RECEIVER.receiver_peak_efficiency)] ;

% HEX Parameters:
xHP_aproach = ['HP Approach [*C]: ',num2cell(FIELD.HEX.HP.approach)];
xHP_pinch = ['HP Pinch [*C]: ',num2cell(FIELD.HEX.HP.pinch)];
xHP_pressure = ['HP Pressure [BarA]: ',num2cell(FIELD.HEX.HP.pressure)];
xHP_temperature = ['HP Temperature [*C]: ',num2cell(FIELD.HEX.HP.temperature)];
xMP_aproach = ['MP Approach [*C]: ',num2cell(FIELD.HEX.MP.approach)];
xMP_pinch = ['MP Pinch [*C]: ',num2cell(FIELD.HEX.MP.pinch)];
xMP_pressure = ['MP Pressure [BarA]: ',num2cell(FIELD.HEX.MP.pressure)];
xMP_temperature = ['MP Temperature [*C]: ',num2cell(FIELD.HEX.MP.temperature)];
xLP_aproach = ['LP Approach [*C]: ',num2cell(FIELD.HEX.LP.approach)];
xLP_pinch = ['LP Pinch [*C]: ',num2cell(FIELD.HEX.LP.pinch)];
xLP_pressure = ['LP Pressure [BarA]: ',num2cell(FIELD.HEX.LP.pressure)];
xLP_temperature = ['LP Temperature [*C]: ',num2cell(FIELD.HEX.LP.temperature)];
xwater_supply_pressure = ['Water Supply Pressure [BarA]: ',num2cell(FIELD.HEX.WS.pressure)];
xwater_supply_temperature = ['Water Supply Temperature [*C]: ',num2cell(FIELD.HEX.WS.temperature)];
xpower_block_peak = ['PB Peak Efficiency: ',num2cell(FIELD.HEX.power_block_peak)];
xdesign_inlet_temperature = ['HEX Design Inlet Temperature [*C]: ',num2cell(FIELD.HEX.design_inlet_temperature)];

if FIELD.HEX.plant_type==1
    xplant_type = {'Plant Type: ','Stand-Alone'};
elseif FIELD.HEX.plant_type==2
    xplant_type = {'Plant Type: ','Fixed Efficiency'};
elseif FIELD.HEX.plant_type==3
    xplant_type = {'Plant Type: ','Combined Cycle'};
end

xINPUT = [xfile_name;{'',''};{'Site Parameters: ',''};xsite_name;xlatitude;xaltitude;{'',''};{'Field Parameters: ',''};xfield_azimuth_or;xfield_elevation_or;xN_dishes_per_cluster;xN_clusters_in_field;...
    xN_dishes_in_field;xlns;xlew;xalpha;xHEX_distance_from_field_center;xsystem_pressure;xR1;xR2;xHOT_ins1;xHOT_ins2;xHOT_ins3;xCOLD_ins1;xCOLD_ins2;xCOLD_ins3;xHTS;{'',''};{'Dish Parameters: ',''};xDish_effective_Area;xSlope_Error;xreflectivity;{'',''};...
    {'Receiver Parameters: ',''};xconcentration_ratio;xreceiver_peak_efficiency;{'',''};{'HEX Parameters: ',''};xHP_aproach;xHP_pinch;xHP_pressure;xHP_temperature;xMP_aproach;xMP_pinch;...
    xMP_pressure;xMP_temperature;xLP_aproach;xLP_pinch;xLP_pressure;xLP_temperature;xwater_supply_pressure;xwater_supply_temperature;xpower_block_peak;xdesign_inlet_temperature;xplant_type];

%% SUMMARY:

% Model Run Summary:
field_area = PIPES.field_size;
tot_DNI = sum(DATA.DNI);
tot_GROSS = sum(GROSS_Electric);
tot_consumption = sum(total_consumption);
tot_NET = sum(NET_Electric);
tot_consumption_from_GROSS = tot_consumption./tot_GROSS;
tot_steam = sum(mdot_total.*DATA.step_duration./3600);
avg_physical_shading = 1-sum(DISH.non_shaded)./numel(DISH.non_shaded(DISH.non_shaded>0));
shading_impact_approx = (annual_shading_output_impact_approximation);
avg_dish_efficiency = sum(DISH.efficiency)./numel(DISH.efficiency(DISH.efficiency>0));
avg_receiver_efficiency = sum(receiver_efficiency)./numel(receiver_efficiency(receiver_efficiency>0));
avg_HEX_efficiency = sum(Thermal_efficiency)./numel(Thermal_efficiency(Thermal_efficiency>0));
avg_PB_efficiency = mean(PB_efficiency(PB_efficiency>0));
avg_system_efficiency = sum(NET_Electric)./sum(PIPES.dishes_in_field.*Dish_effective_Area.*DATA.DNI./1000);

xfield_area = ['Field Size [m2]: ',num2cell(field_area)];
xtot_DNI = ['DNI Sum [Wh/m2]: ',num2cell(tot_DNI)];
xtot_GROSS = ['GROSS (tot) [KWh]: ',num2cell(tot_GROSS)];
xtot_consumption = ['Consumption (tot) [KWh]:',num2cell(tot_consumption)];
xtot_consumption_from_GROSS = ['Consumption from GROSS: ',num2cell(tot_consumption_from_GROSS)];
xtot_NET = ['NET (tot) [KWh]: ',num2cell(tot_NET)];
xtot_steam = ['Steam (tot) [ton]: ',num2cell(tot_steam)];
xavg_physical_shading = ['Average Physical Shading: ',num2cell(avg_physical_shading)];
xshading_impact_approx = ['Average Practical Shading: ',num2cell(shading_impact_approx)];
xavg_dish_efficiency = ['Average Dish Efficiency: ',num2cell(avg_dish_efficiency)];
xavg_receiver_efficiency = ['Average Receiver Efficiency: ',num2cell(avg_receiver_efficiency)];
xavg_HEX_efficiency = ['Average HEX Efficiency: ',num2cell(avg_HEX_efficiency)];
xavg_PB_efficiency = ['Average Power Block Efficiency: ',num2cell(avg_PB_efficiency)];
xavg_system_efficiency = ['Average HF System Efficiency: ',num2cell(avg_system_efficiency)];
xstep_duration = ['Step [sec]: ',num2cell(DATA.step_duration)];
xtitle = {'Model Run Summary: ',''};

xSUMMARY_Annual = [xtitle;xfield_area;xstep_duration;xtot_DNI;xtot_GROSS;xtot_consumption;xtot_consumption_from_GROSS;xtot_NET;xtot_steam;...
    xavg_physical_shading;xshading_impact_approx;xavg_dish_efficiency;xavg_receiver_efficiency;xavg_HEX_efficiency;xavg_PB_efficiency;xavg_system_efficiency]; % 13 x 2

% Design Point Summary:
design_DATA.Azimuth = 180;
design_DATA.Elevation = 90;
design_DATA.DNI = 950; xdesign_DNI = ['Design DNI [W/m2]: ',num2cell(design_DATA.DNI)];
design_DATA.wind_velocity = 0; xdesign_wind_velocity = ['Design Wind Velocity [m/sec]: ',num2cell(design_DATA.wind_velocity)];
design_DATA.wind_direction = 180; xdesign_wind_direction = ['Design Wind Direction [Az.deg]: ',num2cell(design_DATA.wind_direction)];
design_DATA.ambient_temp = 20; xdesign_ambient_temp = ['Design Ambient Temperature [*C]: ',num2cell(design_DATA.ambient_temp)];
design_DATA.humidity = 25; xdesign_humidity = ['Design Humidity [%]: ',num2cell(design_DATA.humidity)];
xtitle = {'Design Point: ',''};

xSUMMARY_point = [xtitle;xdesign_DNI;xdesign_wind_velocity;xdesign_wind_direction;xdesign_ambient_temp;xdesign_humidity]; % 8 x 2

design_DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,design_DATA);

xdesign_number_of_dishes = ['Number of Dishes in Field: ',num2cell(N_dishes_per_cluster.*N_clusters_in_field)];
xdesign_effstruct = ['Design Structural Intercept: ',num2cell(design_DISH.effstruct)];
xdesign_efficiency = ['Design Dish Efficiency (tot): ',num2cell(design_DISH.efficiency)];
xdesign_dish_power_pr_sqr_m = ['Design Dish Output [W/m2]: ',num2cell(design_DISH.dish_power_pr_sqr_m)];
xdesign_dish_power_tot = ['Design Dish Output (tot) [KW]: ',num2cell(design_DISH.dish_power_pr_sqr_m.*Dish_effective_Area./1000)];
xtitle = {'Design Dish: ',''};

xSUMMARY_Dish = [xtitle;xdesign_number_of_dishes;xdesign_effstruct;xdesign_efficiency;xdesign_dish_power_pr_sqr_m;xdesign_dish_power_tot]; % 6 x 2

design_FIELD = steady_field_pipes(design_DISH.dish_power_pr_sqr_m,FIELD,design_DATA.wind_velocity,design_DATA.ambient_temp,design_DATA.humidity);

        design_Tin_receiver = design_FIELD.PIPES.SECTION{1}.cold_temperature1; xdesign_Tin_receiver = ['Design Tin Receiver [*C]: ',num2cell(design_Tin_receiver)];
        design_Tout_receiver = design_FIELD.PIPES.SECTION{1}.hot_temperature1; xdesign_Tout_receiver = ['Design Tout Receiver [*C]: ',num2cell(design_Tout_receiver)];
        design_Tin_HEX = design_FIELD.PIPES.SECTION{N}.hot_temperature2; xdesign_Tin_HEX = ['Design Tin HEX [*C]: ',num2cell(design_Tin_HEX)];
        design_Tout_HEX = design_FIELD.PIPES.SECTION{N}.cold_temperature2; xdesign_Tout_HEX = ['Design Tout HEX [*C]: ',num2cell(design_Tout_HEX)];
        design_dP_cluster = 100000.*design_FIELD.PIPES.dP_cluster; xdesign_dP_cluster = ['Design Cluster dP [Pa]: ',num2cell(design_dP_cluster)];
        design_mdot = design_FIELD.RECEIVER.mdot; xdesign_mdot = ['Design m* [Kg/sec]: ',num2cell(design_mdot)];
        design_receiver_efficiency = design_FIELD.RECEIVER.efficiency; xdesign_receiver_efficiency = ['Design Receiver Efficiency: ',num2cell(design_receiver_efficiency)];
        design_dP_receiver = 100000.*design_FIELD.PIPES.dP_HEX; xdesign_dP_receiver = ['Design Receiver dP [Pa]: ',num2cell(design_dP_receiver)];
        design_dP_HEX = 100000.*design_FIELD.PIPES.dP_HEX; xdesign_dP_HEX = ['Design HEX dP [Pa]: ',num2cell(design_dP_HEX)];
        xtitle = {'PIPES Endpoints: ',''};
        
xSUMMARY_PIPES_Points = [xtitle;xdesign_mdot;xdesign_receiver_efficiency;xdesign_Tin_receiver;xdesign_Tout_receiver;xdesign_dP_receiver;xdesign_Tin_HEX;xdesign_Tout_HEX;xdesign_dP_HEX;xdesign_dP_cluster]; % 10 x 2
              
        design_Tout_Receiver = design_FIELD.PIPES.SECTION{1}.hot_temperature1;
        design_Tin_Receiver = design_FIELD.PIPES.SECTION{1}.cold_temperature1;        
        
         for j = 1:N
            design_Tout_HOT_sections(j) = design_FIELD.PIPES.SECTION{j}.hot_temperature2; 
            design_dTinstruments_HOT_sections(j) = design_FIELD.PIPES.SECTION{j}.hot_temperature_loss_by_components;
            design_Tout_COLD_sections(j) = design_FIELD.PIPES.SECTION{j}.cold_temperature2; 
            design_dTinstruments_COLD_sections(j) = design_FIELD.PIPES.SECTION{j}.cold_temperature_loss_by_components;
            design_dP_HOT_sections(j) = 100000.*design_FIELD.PIPES.SECTION{j}.dP_hot; 
            design_dP_COLD_sections(j) = 100000.*design_FIELD.PIPES.SECTION{j}.dP_cold; 
            length(j) = design_FIELD.PIPES.SECTION{j}.length; 
            hot_elbows(j) = design_FIELD.PIPES.SECTION{j}.hot_elbows; 
            cold_elbows(j) = design_FIELD.PIPES.SECTION{j}.cold_elbows; 
            tag(j) = {design_FIELD.PIPES.SECTION{j}.tag};
         end
         
         tag = ['Receiver',tag];
         xdesign_section_tag = ['Section Tag: ',tag];
         xsection_length = ['Section Length [m]: ',num2cell([0,length])];
         xdesign_section_hot_elbows = ['Number of  Hot Elbows: ',num2cell([0,hot_elbows])];
         xdesign_section_cold_elbows = ['Number of  Cold Elbows: ',num2cell([0,cold_elbows])];
         xdesign_Tout_HOT_sections = ['Design Thot out [*C]: ',num2cell([design_Tout_Receiver,design_Tout_HOT_sections])];
         xdesign_dTinstruments_HOT_sections = ['Hot Instruments Temperature Loss [*C]: ',num2cell([0,design_dTinstruments_HOT_sections])];
         xdesign_Tout_COLD_sections = ['Design Tcold out [*C]: ',num2cell([design_Tin_Receiver,design_Tout_COLD_sections])];
         xdesign_dTinstruments_COLD_sections = ['Cold Instruments Temperature Loss [*C]: ',num2cell([0,design_dTinstruments_COLD_sections])];
         xdesign_dP_HOT_sections = ['Design dP hot [Pa]: ',num2cell([1500.*(design_FIELD.RECEIVER.mdot./0.6).^2,design_dP_HOT_sections])];
         xdesign_dP_COLD_sections = ['Design dP cold [Pa]: ',num2cell([1500.*(design_FIELD.RECEIVER.mdot./0.6).^2,design_dP_COLD_sections])];
         xtitle = ['PIPES Sections: ',num2cell(nan(1,N+1))];
         
xSUMMARY_PIPES_Sections = [xtitle;xdesign_section_tag;xsection_length;xdesign_section_hot_elbows;xdesign_section_cold_elbows;xdesign_Tout_HOT_sections;xdesign_dTinstruments_HOT_sections;xdesign_Tout_COLD_sections;xdesign_dTinstruments_COLD_sections;xdesign_dP_HOT_sections;xdesign_dP_COLD_sections]; % 8xN+1

        design_dP_blowers(1,1:size(design_FIELD.PIPES.dP_blowers,2)) = 100000.*design_FIELD.PIPES.dP_blowers; 
        p1 = ['Design Blowers {'  num2str(size(design_FIELD.PIPES.dP_blowers,2)) '} dP [Pa]: '];
        xdesign_dP_blowers = [p1,num2cell(design_dP_blowers)];
        design_blowers_consumption(1,1:size(design_FIELD.PIPES.blowers_consumption,2)) = design_FIELD.PIPES.blowers_consumption;
        p2 = ['Design Blowers {' num2str(size(design_FIELD.PIPES.dP_blowers,2)) '} Consumption [KW]: '];
        xdesign_blowers_consumption = [p2,num2cell(design_blowers_consumption)];
        design_all_blowers_consumption = design_FIELD.PIPES.all_blowers_consumption;
        xdesign_all_blowers_consumption = ['Design Blowers Consumption (tot) [KW]: ',num2cell(design_all_blowers_consumption),num2cell(nan(1,size(design_FIELD.PIPES.blowers_consumption,2)-1))];
        xtitle = ['Blowers: ',num2cell(nan(1,size(design_FIELD.PIPES.dP_blowers,2)))];

xSUMMARY_Blowers = [xtitle;xdesign_dP_blowers;xdesign_blowers_consumption;xdesign_all_blowers_consumption]; % 4 x Number of clusters
    
        design_mdot_HP = design_FIELD.POWER.mdot_HP; xdesign_mdot_HP = ['Design m* HP [ton/hr]: ',num2cell(design_mdot_HP)];
        design_mdot_MP =design_FIELD.POWER.mdot_MP; xdesign_mdot_MP = ['Design m* MP [ton/hr]: ',num2cell(design_mdot_MP)];
        design_mdot_LP = design_FIELD.POWER.mdot_LP; xdesign_mdot_LP = ['Design m* LP [ton/hr]: ',num2cell(design_mdot_LP)];
        design_mdot_total = (design_FIELD.POWER.mdot_HP+design_FIELD.POWER.mdot_MP+design_FIELD.POWER.mdot_LP); xdesign_mdot_total = ['Design m* steam (tot) [tom/hr]: ',num2cell(design_mdot_total)];
        design_Thermal_efficiency = design_FIELD.POWER.Thermal_efficiency; xdesign_Thermal_efficiency = ['Design HEX Thermal Efficiency: ',num2cell(design_Thermal_efficiency)];
        design_Thermal = design_FIELD.POWER.Thermal; xdesign_Thermal = ['Design Thermal Power [KWth]: ',num2cell(design_Thermal)];
        design_PB_efficiency = design_FIELD.POWER.PB_efficiency; xdesign_PB_efficiency = ['Design PB Efficiency: ',num2cell(design_PB_efficiency)];
        design_GROSS_Electric = design_FIELD.POWER.Electric; xdesign_GROSS_Electric = ['Design GROSS Power [KWe]: ',num2cell(design_GROSS_Electric)];
        design_HP_pump_consumption = design_FIELD.POWER.consumption.HP_pump_consumption; xdesign_HP_pump_consumption = ['Design HP Pump Consumption [KW]: ',num2cell(design_HP_pump_consumption)];
        design_MP_pump_consumption = design_FIELD.POWER.consumption.MP_pump_consumption; xdesign_MP_pump_consumption = ['Design MP Pump Consumption [KW]: ',num2cell(design_MP_pump_consumption)];
        design_LP_pump_consumption = design_FIELD.POWER.consumption.LP_pump_consumption; xdesign_LP_pump_consumption = ['Design LP Pump Consumption [KW]: ',num2cell(design_LP_pump_consumption)];
        design_all_pumps_consumption = design_FIELD.POWER.consumption.all_pumps_consumption; xdesign_all_pumps_consumption = ['Design Pumps Consumption (tot) [KW]: ',num2cell(design_all_pumps_consumption)];
        design_NET_Electric = design_GROSS_Electric - (PIPES.dishes_in_field.*design_DISH.consumption+design_all_blowers_consumption+design_all_pumps_consumption); xdesign_NET_Electric = ['Installed Capacity (NET design) [KWe]: ',num2cell(design_NET_Electric)];
        design_system_efficiency = design_NET_Electric./(PIPES.dishes_in_field.*design_DATA.DNI.*Dish_effective_Area./1000); xdesign_system_efficiency = ['Design System Efficiency: ',num2cell(design_system_efficiency)];
        xtitle = {'PB Design Performance: ',''};
        
xSUMMARY_PB = [xtitle;xdesign_mdot_HP;xdesign_mdot_MP;xdesign_mdot_LP;xdesign_mdot_total;xdesign_Thermal_efficiency;xdesign_Thermal;xdesign_PB_efficiency;xdesign_GROSS_Electric;...
    xdesign_HP_pump_consumption;xdesign_MP_pump_consumption;xdesign_LP_pump_consumption;xdesign_all_pumps_consumption;'';xdesign_NET_Electric;'';xdesign_system_efficiency]; % 15 x 2

[status,msginfo] = xlswrite(file_name,xSUMMARY_Annual,'SUMMARY','A4');
[status,msginfo] = xlswrite(file_name,xSUMMARY_point,'SUMMARY','A21');
[status,msginfo] = xlswrite(file_name,xSUMMARY_PB,'SUMMARY','A28');
[status,msginfo] = xlswrite(file_name,xSUMMARY_Dish,'SUMMARY','D4');
[status,msginfo] = xlswrite(file_name,xSUMMARY_PIPES_Points,'SUMMARY','D11');
[status,msginfo] = xlswrite(file_name,xSUMMARY_Blowers,'SUMMARY','D22');
[status,msginfo] = xlswrite(file_name,xSUMMARY_PIPES_Sections,'SUMMARY','D27');

%% DATA:
xtime = ['Time';DATA.time];
xsun_azimuth = ['Sun Azimuth [Deg]';num2cell(DATA.Azimuth)];         
xsun_elevation = ['Sun Elevation [Deg]';num2cell(DATA.Elevation)];         
xDNI = ['DNI [W/m2]';num2cell(DATA.DNI)];  
xwind_velocity = ['Wind Velocity [m/sec]';num2cell(DATA.wind_velocity)];  
xwind_direction = ['Wind Direction [Deg]';num2cell(DATA.wind_direction)];  
xambient_temp = ['Ambient Temperature [*C]';num2cell(DATA.ambient_temp)];  
xhumidity = ['Humidity [%]';num2cell(DATA.humidity)];  
xsingularity_count = ['Model Singularity Conditions';num2cell(singularity_count')];

xDATA = [xtime,xsun_azimuth,xsun_elevation,xDNI,xwind_velocity,xwind_direction,xambient_temp,xhumidity,xsingularity_count];

%% DISH:
xnon_shaded = ['Non Shaded (1-shade)';num2cell(DISH.non_shaded)];
xeffstruct = ['Structiral Intercept';num2cell(DISH.effstruct)];
xeffdust = ['Reflectivity (1-dirt)';num2cell(DISH.effdust)];
xefficiency = ['Efficiency (tot)';num2cell(DISH.efficiency)];
xdish_power_pr_sqr_m = ['Dish Output [W/m2]';num2cell(DISH.dish_power_pr_sqr_m)];
xdish_power_tot = ['Dish Output (tot) [KW]';num2cell(DISH.dish_power_pr_sqr_m.*Dish_effective_Area./1000)];
xdish_consumption = ['Single Dish Consumption [KW]';num2cell(DISH.consumption)];
xall_dishes_consumption = ['All Dishes Consumption [KW]';num2cell(DISH.consumption.*PIPES.dishes_in_field)];

xDISH = [xnon_shaded,xeffstruct,xeffdust,xefficiency,xdish_power_pr_sqr_m,xdish_power_tot,xdish_consumption,xall_dishes_consumption];

%% RECEIVER:
xreceiver_efficiency = ['Receiver Efficiency';num2cell(receiver_efficiency')];
xmdot = ['m* [Kg/sec]';num2cell(mdot')];
xTin_receiver = ['Receiver Tin [*C]';num2cell(Tin_receiver')];
xTout_receiver = ['Receiver Tout [*C]';num2cell(Tout_receiver')];
xdP_receiver = ['Receiver dP [Pa]';num2cell(dP_receiver'.*100000)];

xRECEIVER = [xdish_power_tot,xreceiver_efficiency,xmdot,xTin_receiver,xTout_receiver,xdP_receiver];

%% PIPES:
for j = 1:N

        xTout_HOT_sections(:,j) =  ['Section ' FIELD.PIPES.SECTION{j}.tag ' Thot out [*C]';num2cell(Tout_HOT_sections(:,j))];
        xdTinstruments_HOT_sections(:,j) =  ['Section ' FIELD.PIPES.SECTION{j}.tag ' Hot Instruments Loss [*C]';num2cell(dTinstruments_HOT_sections(:,j))];
        xTout_COLD_sections(:,j) =  ['Section ' FIELD.PIPES.SECTION{j}.tag ' Tcold in [*C]';num2cell(Tout_COLD_sections(:,j))];
        xdTinstruments_COLD_sections(:,j) =  ['Section ' FIELD.PIPES.SECTION{j}.tag ' Cold Instruments Loss [*C]';num2cell(dTinstruments_COLD_sections(:,j))];
        xdP_HOT_sections(:,j) = ['Section ' FIELD.PIPES.SECTION{j}.tag ' hot dP [Pa]';num2cell(100000.*dP_HOT_sections(:,j))];
        xdP_COLD_sections(:,j) = ['Section ' FIELD.PIPES.SECTION{j}.tag ' cold dP [Pa]';num2cell(100000.*dP_COLD_sections(:,j))];

end

xall_blowers_consumption = ['All Blowers Consumption [KW]';num2cell(all_blowers_consumption')];
xdP_cluster = ['Cluster dP [Pa]';num2cell(100000.*dP_cluster')];
xdP_HEX = ['HEX dP [Pa]';num2cell(100000.*dP_HEX')];

for k = 1:size(FIELD.PIPES.dP_blowers,2)

        xdP_blowers(:,k) = [['Blower ' num2str(k) ' dP [Pa]'];num2cell(100000.*dP_blowers(:,k))];
        xblowers_consumption(:,k) = [['Blower ' num2str(k) ' Consumption [KW]'];num2cell(blowers_consumption(:,k))];

end

xPIPES = [xmdot,xTout_receiver,xTout_HOT_sections,xTin_receiver,xTout_COLD_sections,...
    xdTinstruments_HOT_sections,xdTinstruments_COLD_sections,xdP_receiver,xdP_HOT_sections,xdP_COLD_sections,xdP_HEX,xdP_cluster,xdP_blowers,...
    xblowers_consumption,xall_blowers_consumption];

%% HEX:
xTin_HEX = ['HEX Air Tin [*C]';num2cell(Tin_HEX')];
xTout_HEX = ['HEX Air Tout [*C]';num2cell(Tout_HEX')];
xmdot_HP = ['HP m* [ton/hr]';num2cell(mdot_HP')];
xmdot_MP = ['MP m* [ton/hr]';num2cell(mdot_MP')];
xmdot_LP = ['LP m* [ton/hr]';num2cell(mdot_LP')];
xmdot_total = ['(tot) m* [ton/hr]';num2cell(mdot_total')];
xThermal= ['Thermal Power [KWth]';num2cell(Thermal')];
xThermal_efficiency = ['Thermal Efficiency';num2cell(Thermal_efficiency')];
xPB_efficiency = ['Power Block Efficiency';num2cell(PB_efficiency')];
xGROSS_Electric= ['GROSS Production [KWe]';num2cell(GROSS_Electric')];
xNET_Electric = ['NET Production [KWe]';num2cell(NET_Electric')];
xsystem_efficiency = ['HF System Efficiency';num2cell(system_efficiency)];
xHP_pump_consumption = ['HP Pump Consumption [KW]';num2cell(HP_pump_consumption')];
xMP_pump_consumption = ['MP Pump Consumption [KW]';num2cell(MP_pump_consumption')];
xLP_pump_consumption = ['LP Pump Consumption [KW]';num2cell(LP_pump_consumption')];
xall_pumps_consumption = ['All Pumps Consumption [KW]';num2cell(all_pumps_consumption')];

xPB = [xmdot,xTin_HEX,xTout_HEX,xmdot_HP,xmdot_MP,xmdot_LP,xmdot_total,xThermal,xThermal_efficiency,...
    xPB_efficiency,xGROSS_Electric,xNET_Electric,xsystem_efficiency,xHP_pump_consumption,xMP_pump_consumption,...
    xLP_pump_consumption,xall_pumps_consumption];

%% CONSUMPTION:
% 
% xtotal_consumption = ['Field Consumption [KW]';num2cell(total_consumption)];
% 
% xCONSUMPTION = [xdish_consumption,xall_dishes_consumption,xblowers_consumption,xall_blowers_consumption,...
%     xHP_pump_consumption,xMP_pump_consumption,xLP_pump_consumption,xall_pumps_consumption];

%% COST:

xCOST_single_dish_cost = ['Single Dish Cost [$]:',num2cell(PIPES.single_dish_cost)];
xCOST_all_dishes_cost = [[num2str(PIPES.dishes_in_field) ' Dishes in the Field - Total Cost [$]:'],num2cell(PIPES.all_dishes_cost)];
xCOST_steam_cycle_cost = ['Steam Cycle Cost [$]:',num2cell(PIPES.steam_cycle_cost)];
xCOST_power_unit_cost = ['Power Unit Cost [$]:',num2cell(PIPES.power_unit_cost)];
xCOST_pumps_cost = ['Pumps Cost [$]:',num2cell(PIPES.pumps_cost)];
xCOST_water_treatment_system_cost = ['Water Treatment System Cost [$]:',num2cell(PIPES.water_treatment_system_cost)];
xCOST_air_system_cost = ['Air System Cost [$]:',num2cell(PIPES.air_system_cost)];
xCOST_blowers_cost = ['Blowers Cost [$]:',num2cell(PIPES.blowers_cost)];
xCOST_Instruments_cost = ['Instruments Cost [$]:',num2cell(PIPES.Instruments_cost)];
xCOST_control_cost = ['Control System Cost [$]:',num2cell(PIPES.control_cost)];
xCOST_other_essentials_cost = ['Other Essentials Cost [$]:',num2cell(PIPES.other_essentials_cost)];
xCOST_all_pipes_cost = ['HTS (PIPES) Cost [$]:',num2cell(PIPES.all_pipes_cost)];
xCOST_project_cost = ['Project (tot) Cost [$]:',num2cell(PIPES.project_cost)];

xCOST = [xCOST_single_dish_cost;xCOST_all_dishes_cost;xCOST_steam_cycle_cost;xCOST_power_unit_cost;...
    xCOST_pumps_cost;xCOST_water_treatment_system_cost;xCOST_air_system_cost;xCOST_blowers_cost;...
    xCOST_Instruments_cost;xCOST_control_cost;xCOST_other_essentials_cost;xCOST_all_pipes_cost;...
    xCOST_project_cost];

% PIPES:  x-number of dishes in field



for i = 1:N
    tagC(i) = {design_FIELD.PIPES.SECTION{i}.tag};
    Length_L(i) = FIELD.PIPES.SECTION{i}.length;
    hot_Radii(:,i) = FIELD.PIPES.SECTION{i}.hot_geometry';
    cold_Radii(:,i) = FIELD.PIPES.SECTION{i}.cold_geometry';
    COST_L_stainless_steel(i) = FIELD.PIPES.SECTION{i}.stainless_steel_cost;
    COST_L_carbon_steel(i) = FIELD.PIPES.SECTION{i}.carbon_steel_cost;
    COST_L_hot_insulation1(i) = FIELD.PIPES.SECTION{i}.hot_insulation1_cost;
    COST_L_hot_insulation2(i) = FIELD.PIPES.SECTION{i}.hot_insulation2_cost;
    COST_L_hot_insulation3(i) = FIELD.PIPES.SECTION{i}.hot_insulation3_cost;
    COST_L_aluminium(i) = FIELD.PIPES.SECTION{i}.aluminium_cost;
    COST_L_cold_insulation1(i) = FIELD.PIPES.SECTION{i}.cold_insulation1_cost;
    COST_L_cold_insulation2(i) = FIELD.PIPES.SECTION{i}.cold_insulation2_cost;
    COST_L_cold_insulation3(i) = FIELD.PIPES.SECTION{i}.cold_insulation3_cost;
    COST_L_section_cost(i) = FIELD.PIPES.SECTION{i}.section_cost;
    COST_L_sections_field(i) = FIELD.PIPES.SECTION{i}.sections_field_cost;
    COST_1m_stainless_steel(i) = FIELD.PIPES.SECTION{i}.stainless_steel_cost./Length_L(i);
    COST_1m_carbon_steel(i) = FIELD.PIPES.SECTION{i}.carbon_steel_cost./Length_L(i);
    COST_1m_hot_insulation1(i) = FIELD.PIPES.SECTION{i}.hot_insulation1_cost./Length_L(i);
    COST_1m_hot_insulation2(i) = FIELD.PIPES.SECTION{i}.hot_insulation2_cost./Length_L(i);
    COST_1m_hot_insulation3(i) = FIELD.PIPES.SECTION{i}.hot_insulation3_cost./Length_L(i);
    COST_1m_aluminium(i) = FIELD.PIPES.SECTION{i}.aluminium_cost./Length_L(i);
    COST_1m_cold_insulation1(i) = FIELD.PIPES.SECTION{i}.cold_insulation1_cost./Length_L(i);
    COST_1m_cold_insulation2(i) = FIELD.PIPES.SECTION{i}.cold_insulation2_cost./Length_L(i);
    COST_1m_cold_insulation3(i) = FIELD.PIPES.SECTION{i}.cold_insulation3_cost./Length_L(i);
    COST_1m_section_cost(i) = FIELD.PIPES.SECTION{i}.section_cost./Length_L(i);
    COST_1m_sections_field(i) = FIELD.PIPES.SECTION{i}.sections_field_cost./Length_L(i);
    cold_bypasses(i) = PIPES.SECTION{i}.cold_bypasses;
    cold_bypass_length(i) = PIPES.SECTION{i}.cold_bypass_length;
%     cold_bypass_thickness(i) = PIPES.SECTION{i}.cold_bypass_thickness;
end

number_of_sections = COST_L_sections_field./COST_L_section_cost;

for p = 1:2.*N
    if p./2==ceil(p./2)
        xxtag(p) = tagC(p./2);
        Length(p) = Length_L(p./2);
        hot_Radiix(:,p) = num2cell(hot_Radii(:,p./2));
        cold_Radiix(:,p) = num2cell(cold_Radii(:,p./2));
        hot_Diax(:,p) = num2cell(2.*hot_Radii(:,p./2));
        cold_Diax(:,p) = num2cell(2.*cold_Radii(:,p./2));
        Nsections(p) = num2cell(number_of_sections(p./2));
        COST_stainless_steel(p) = COST_L_stainless_steel(p./2);
        COST_carbon_steel(p) = COST_L_carbon_steel(p./2);
        COST_hot_insulation1(p) = COST_L_hot_insulation1(p./2);
        COST_hot_insulation2(p) = COST_L_hot_insulation2(p./2);
        COST_hot_insulation3(p) = COST_L_hot_insulation3(p./2);
        COST_cold_insulation1(p) = COST_L_cold_insulation1(p./2);
        COST_cold_insulation2(p) = COST_L_cold_insulation2(p./2);
        COST_cold_insulation3(p) = COST_L_cold_insulation3(p./2);
        COST_aluminium(p) = COST_L_aluminium(p./2);
        COST_section_cost(p) = COST_L_section_cost(p./2);
        COST_sections_field(p) = COST_L_sections_field(p./2);
        xcold_bypasses(p) = num2cell(cold_bypasses(p./2));
        xcold_bypass_length(p) = num2cell(cold_bypass_length(p./2));
%         xcold_bypass_thickness(p) = num2cell(cold_bypass_thickness(p./2));
    else
        xxtag(p) = {''};
        Length(p) = 1;
        cold_Radiix(:,p) = {'';'';'';'';''};
        hot_Radiix(:,p) = {'';'';'';'';''};
        cold_Diax(:,p) = {'';'';'';'';''};
        hot_Diax(:,p) = {'';'';'';'';''};
        Nsections(p) = {''};
        COST_stainless_steel(p) = COST_1m_stainless_steel(ceil(p./2));
        COST_carbon_steel(p) = COST_1m_carbon_steel(ceil(p./2));
        COST_hot_insulation1(p) = COST_1m_hot_insulation1(ceil(p./2));
        COST_hot_insulation2(p) = COST_1m_hot_insulation2(ceil(p./2));
        COST_hot_insulation3(p) = COST_1m_hot_insulation3(ceil(p./2));
        COST_cold_insulation1(p) = COST_1m_cold_insulation1(ceil(p./2));
        COST_cold_insulation2(p) = COST_1m_cold_insulation2(ceil(p./2));
        COST_cold_insulation3(p) = COST_1m_cold_insulation3(ceil(p./2));
        COST_aluminium(p) = COST_1m_aluminium(ceil(p./2));
        COST_section_cost(p) = COST_1m_section_cost(ceil(p./2));
        COST_sections_field(p) = COST_1m_sections_field(ceil(p./2));
        xcold_bypasses(p) = {''};
        xcold_bypass_length(p) = {''};
%         xcold_bypass_thickness(p) = {''};
    end

end

xtag = ['Section Tag [m]: ',xxtag];
xLength = ['Section Length [m]: ',num2cell(Length)];
xRadiix = [[{'R hot pipe in [m]: ';'R hot pipe out [m]: ';'R hot insulation 1 [m]: ';'R hot insulation 2 [m]: ';'R hot insulation 3 [m]: '},hot_Radiix];[{'R cold pipe in [m]: ';'R cold pipe out [m]: ';'R cold insulation 1 [m]: ';'R cold insulation 2 [m]: ';'R cold insulation 3 [m]: '},cold_Radiix]];
xDiax = [[{'D hot pipe in [m]: ';'D hot pipe out [m]: ';'D hot insulation 1 [m]: ';'D hot insulation 2 [m]: ';'D hot insulation 3 [m]: '},cold_Diax];[{'D cold pipe in [m]: ';'D cold pipe out [m]: ';'D cold insulation 1 [m]: ';'D cold insulation 2 [m]: ';'D cold insulation 3 [m]: '},cold_Diax]];
xNelements = ['Number in Field: ',Nsections];
xCOST_stainless_steel = ['Stainless Steel [$]: ',num2cell(COST_stainless_steel)];
xCOST_carbon_steel = ['Carbon Steel [$]: ',num2cell(COST_carbon_steel)];
xCOST_aluminium = ['Aluminium [$]: ',num2cell(COST_aluminium)];

for j = 1:2
    for i = 1:3
            if j==1
                hot_insulation{i} = ['HOT Layer ' num2str(i) ' Insulation Matereial: ' insulation_type_inv(INPUT.HTS_hot.insulation_material(i)) 'Cost [$]: ']
            else
                cold_insulation{i} = ['COLD Layer ' num2str(i) ' Insulation Matereial: ' insulation_type_inv(INPUT.HTS_hot.insulation_material(i)) 'Cost [$]: ']
            end
    end
end

xCOST_hot_insulation1 = [hot_insulation{1},num2cell(COST_hot_insulation1)];
xCOST_hot_insulation2 = [hot_insulation{2},num2cell(COST_hot_insulation2)];
xCOST_hot_insulation3 = [hot_insulation{3},num2cell(COST_hot_insulation3)];

xCOST_cold_insulation1 = [cold_insulation{1},num2cell(COST_cold_insulation1)];
xCOST_cold_insulation2 = [cold_insulation{2},num2cell(COST_cold_insulation2)];
xCOST_cold_insulation3 = [cold_insulation{3},num2cell(COST_cold_insulation3)];

xCOST_section_cost = ['Section Cost [$]: ',num2cell(COST_section_cost)];
xCOST_sections_field = ['Sections (tot in field) Cost [$]: ',num2cell(COST_sections_field)];
xxcold_bypasses = ['Number of Cold Bypasses: ',xcold_bypasses];
xxcold_bypass_length = ['Cold Bypasses Length [m]: ',xcold_bypass_length];
% xxcold_bypass_thickness = ['Cold Bypasses Thickness [m]:',xcold_bypass_thickness];

xCOST_Pipes = [xtag;xNelements;xLength;xRadiix;xDiax;xxcold_bypasses;xxcold_bypass_length;xCOST_stainless_steel;...
    xCOST_carbon_steel;xCOST_hot_insulation1;xCOST_hot_insulation2;xCOST_hot_insulation3;xCOST_cold_insulation1;xCOST_cold_insulation2;xCOST_cold_insulation3;...
    xCOST_aluminium;xCOST_section_cost;xCOST_sections_field];


%% MS Excell WRITING:

%% MS Excell WRITING:

[status,msginfo] = xlswrite(file_name,xINPUT,'INPUT','A3');
[status,msginfo] = xlswrite(file_name,[xDATA,xDISH],'DISH','A3');
[status,msginfo] = xlswrite(file_name,[xDATA,xRECEIVER],'RECEIVER','A3');
[status,msginfo] = xlswrite(file_name,[xDATA,xPIPES],'PIPES','A3');
[status,msginfo] = xlswrite(file_name,[xDATA,xPB],'HEX','A3');
[status,msginfo] = xlswrite(file_name,xCOST,'COST','A3');
[status,msginfo] = xlswrite(file_name,xCOST_Pipes,'COST','A17');
t = toc./60

% HF title:
HF_title = {'Helio Focus. LTD',' * Steady State Model * ','Site: ' site_name,'Date: ' datestr(now),'Calculation Duration : ' datestr(toc/(24.*3600),'HH:MM:SS')};
[status,msginfo] = xlswrite(file_name,HF_title,'INPUT','A1');
[status,msginfo] = xlswrite(file_name,HF_title,'SUMMARY','A1');
[status,msginfo] = xlswrite(file_name,HF_title,'DISH','A1');
[status,msginfo] = xlswrite(file_name,HF_title,'RECEIVER','A1');
[status,msginfo] = xlswrite(file_name,HF_title,'PIPES','A1');
[status,msginfo] = xlswrite(file_name,HF_title,'HEX','A1');
[status,msginfo] = xlswrite(file_name,HF_title,'COST','A1');

OPT_title = {'Project Cost [$]: ' num2str(PIPES.project_cost),'','Net Power [KWh]: ' num2str(tot_NET),'','Optimization Criterion [Net(KWh.e)/Cost($)] : ' num2str(tot_NET/PIPES.project_cost)};
[status,msginfo] = xlswrite(file_name,OPT_title,'INPUT','A2');
[status,msginfo] = xlswrite(file_name,OPT_title,'SUMMARY','A2');
[status,msginfo] = xlswrite(file_name,OPT_title,'DISH','A1');
[status,msginfo] = xlswrite(file_name,OPT_title,'RECEIVER','A2');
[status,msginfo] = xlswrite(file_name,OPT_title,'PIPES','A2');
[status,msginfo] = xlswrite(file_name,OPT_title,'HEX','A2');
[status,msginfo] = xlswrite(file_name,OPT_title,'COST','A2');

f = [file_name '.mat']
save(f);
message = ['Model Run Finished (in ' num2str(t) 'min)'];
waitbar(x,status_window,message);
refreshdata
drawnow expose update
