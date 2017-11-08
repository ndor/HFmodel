function OPT = field_opt_annulus(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
%FIELD_OPT Summary of this function goes here
%   Detailed explanation goes here
global INPUT...
    stop_flag... 
    time_length...
    G_sections_length... 
    G_sections_tag...
    pus...
        toc...
    G_time...
    G_sun_azimuth...
    G_sun_elevation...
    G_DNI...
    G_wind_velocity...
    G_wind_direction...
    G_Tambient...
    G_humidity...
    G_dish_efficiency...
    G_T_hot...
    G_Tcold...
    G_mdot...
    G_receiver_efficiency...
    G_HEX_efficiency...
    G_PB_efficiency...
    G_steam_tot...
    G_all_blowers_consumption...
    G_all_pumps_consumption....
    G_all_dishes_consumption...
    G_gross_electric...
    G_net_electric

if stop_flag<=1
    
    transient_state = 0;

%     X(isnan(X)==1) = X0(isnan(X)==1);    
    field_azimuth_or = X(1);
    field_elevation_or = X(2);
    N_dishes_per_cluster = ceil(X(3)); % number of dishes at each cluster (minimum 1)
    N_clusters_in_field = ceil(X(4)); % number of clusters in the field (minimum 1)
    lns = X(5); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
    lew = X(6); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
    alpha = X(7); % field parralelogram angle (from E-W line, +CW) [deg]
    HEX_distance_from_field_center = 0;%X(5); % HEX distance from center field - HEX linkage [m]

    % HTS Parameters:
        if HTS_type>1
            INPUT.HTS_hot.duct_R = X(8);%INPUT.HTS_hot.duct_R;
            INPUT.HTS_cold.duct_R = X(9);%INPUT.HTS_cold.duct_R;
            
            INPUT.HTS_hot.insulation_occupation(1) = X(10);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_hot.insulation_occupation(2) = X(11);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_hot.insulation_occupation(3) = X(12);%INPUT.HTS_hot.insulation_occupation;
%             hLm = INPUT.HTS_hot.insulation_material;
%             hLq = INPUT.HTS_hot.insulation_quality;

            INPUT.HTS_cold.insulation_thickness(1) = X(13);%INPUT.HTS_cold.insulation_thickness;
            INPUT.HTS_cold.insulation_thickness(2) = X(14);%INPUT.HTS_cold.insulation_thickness;
            INPUT.HTS_cold.insulation_thickness(3) = X(15);%INPUT.HTS_cold.insulation_thickness;
%             cLm = INPUT.HTS_cold.insulation_material;
%             cLq = INPUT.HTS_cold.insulation_quality;
        else
            INPUT.HTS_hot.duct_R = X(8);%INPUT.HTS_hot.duct_R;
            INPUT.HTS_cold.duct_R = X(9);%INPUT.HTS_cold.duct_R;
            
            INPUT.HTS_hot.insulation_thickness(1) = X(10);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_hot.insulation_thickness(2) = X(11);%INPUT.HTS_hot.insulation_occupation;
            INPUT.HTS_hot.insulation_thickness(3) = X(12);%INPUT.HTS_hot.insulation_occupation;
%             hLm = INPUT.HTS_hot.insulation_material;
%             hLq = INPUT.HTS_hot.insulation_quality;

            INPUT.HTS_cold.insulation_thickness(1) = X(13);%INPUT.HTS_cold.insulation_thickness;
            INPUT.HTS_cold.insulation_thickness(2) = X(14);%INPUT.HTS_cold.insulation_thickness;
            INPUT.HTS_cold.insulation_thickness(3) = X(15);%INPUT.HTS_cold.insulation_thickness;
%             cLm = INPUT.HTS_cold.insulation_material;
%             cLq = INPUT.HTS_cold.insulation_quality;
        end

        
        
        
            ANNULUS = field_structure_annulus(transient_state);
%                     ANNULUS.pressure
%                     ANNULUS.dishes_per_cluster
%                     ANNULUS.clusters_in_field
%                     ANNULUS.dishes_in_field
%                     ANNULUS.dt
%                     ANNULUS.N
            N = ANNULUS.N;
            FIELD.ANNULUS = ANNULUS;

            
%             go_optimize = GONOGO(0,FIELD,OPT_USER,Slope_Error,reflectivity,field_azimuth_or,field_elevation_or,0,0,0);
%             
%             
%             if go_optimize==0
%                 
%                 OPT = 0;
%                 
%             else
            
                % GLOBALS:
                G_sections_length(1) = 0;
                G_sections_tag{1} = {'receiver'};
                
                         for j = 1:N
                             
                             
                            if stop_flag==1
                                delete(status_window);
                                break
                            end
                            
                            G_sections_length(j+1) = ANNULUS.SECTION{j}.length+G_sections_length(j);
                            G_sections_tag{j+1} = {ANNULUS.SECTION{j}.tag};
   
                         end


                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 time_length = size(DISH.dish_power_pr_sqr_m,1);  
                 N = size(ANNULUS.SECTION,2);
                 FIELD.ANNULUS = ANNULUS;

                 % Model calculation:

% %                     F = @ (dish_power_pr_sqr_m,FIELD,wind_velocity,ambient_temp,humidity) steady_field(dish_power_pr_sqr_m,FIELD{:},wind_velocity,ambient_temp,humidity);
% %                     NEF=cell(time_length,1);
% %                     for i=1:time_length; NEF{i} = FIELD; end
% % 
% %                     FIELD = arrayfun(F,DISH.dish_power_pr_sqr_m,NEF,DATA.wind_velocity,DATA.ambient_temp,DATA.humidity,'UniformOutput',false);
% %                     for i=1:time_length
% % 
% %                         Steam(i) = real(FIELD{i}.POWER.mdot_HP+FIELD{i}.POWER.mdot_MP+FIELD{i}.POWER.mdot_LP);
% %                         Consumption(i) = (ANNULUS.dishes_in_field.*DISH.consumption(i)+FIELD{i}.ANNULUS.all_blowers_consumption+FIELD{i}.POWER.consumption.all_pumps_consumption); 
% %                         NET_Electric(i) = real(DATA.DNI_commonality(i).*(FIELD{i}.POWER.Electric - (ANNULUS.dishes_in_field.*DISH.consumption(i)+FIELD{i}.ANNULUS.all_blowers_consumption+FIELD{i}.POWER.consumption.all_pumps_consumption)));
% % 
% %                     end
                    
                    for i=1:time_length

                        FIELD = steady_field_annulus(DISH.dish_power_pr_sqr_m(i),FIELD,DATA.wind_velocity(i),DATA.ambient_temp(i),DATA.humidity(i));
                        Steam(i) = (DATA.DNI_commonality(i).*(FIELD.POWER.mdot_HP+FIELD.POWER.mdot_MP+FIELD.POWER.mdot_LP));
                        Consumption(i) = DATA.DNI_commonality(i).*(ANNULUS.dishes_in_field.*DISH.consumption(i)+FIELD.ANNULUS.all_blowers_consumption+FIELD.POWER.consumption.all_pumps_consumption); 
                        NET_Electric(i) = (DATA.DNI_commonality(i).*FIELD.POWER.Electric - Consumption(i));

                            refreshdata
                            drawnow expose update   
                    end
                    
%                     NET_Electric(NET_Electric<0) = 0;
                    
                    Steam(isnan(Steam)==1) = 0;
                    Consumption(isnan(Consumption)==1) = 0;
                    NET_Electric(isnan(NET_Electric)==1) = 0;

                    mass_steam_annual = (sum(Steam)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Steam [t/h] production
                    Consumption_annual = (sum(Consumption)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Consumption [KWh] production
                    E_net_annual = (sum(NET_Electric)./sum(DATA.DNI_commonality)).*(0.5.*8762./time_length); % crude estimation of an annual Nete [KWh] production
                    
                     OPT = -sum(NET_Electric)./FIELD.ANNULUS.all_annulus_cost;


%                     go_optimize = GONOGO(1,FIELD,OPT_USER,Slope_Error,reflectivity,field_azimuth_or,field_elevation_or,E_net_annual,mass_steam_annual,Consumption_annual);
%                     
%                     if go_optimize==0
% 
%                         OPT = 0;
% 
%                     else
%                 
%                         OPT = -sum(NET_Electric)./FIELD.ANNULUS.all_annulus_cost;
%                         
%                     end
            
            
%         end

end
    
end

function GO = GONOGO(Pre_or_Post,FIELD,OPT_USER,Slope_Error,reflectivity,field_azimuth_or,field_elevation_or,E_net_annual,mass_steam_annual,Consumption_annual)
%

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

%% Field size (area)
% Limited area or specific span of area – specified by user
% Maximizes the value: A(var)=  (net production)/Cost

% goal function:    

A = FIELD.ANNULUS.field_size;

if (A>A_max || A<A_min) && (A_max~=0 && A_min~=0)
    
    go_optimize(1) = 0;
    
else
    
    go_optimize(1) = 1;
    
end

%% Project cost (CAPEX)
% Limited cost (max cost) – specified by user
% Maximizes the value: B(var)=  (net production)/Cost

% goal function:    

B = FIELD.ANNULUS.project_cost;

if (B>B_max || B<B_min) && (B_max~=0 && B_min~=0)
    
    go_optimize(2) = 0;
    
else
    
    go_optimize(2) = 1;
    
end

%% Annual Energy (Enet)
% Minimal Enet – specified by user
% Maximizes the value: C(var)=  (net production)/Cost

% goal function:    

C = E_net_annual;

if Pre_or_Post==1
    
    if (C>C_max || C<C_min) && (C_max~=0 && C_min~=0)

        go_optimize(3) = 0;

    else

        go_optimize(3) = 1;

    end

else
    
    go_optimize(3) = 1;
    
end

%% Installed Capacity in KW installed
% Limited span of Pnet – specified by user
% Maximizes the value: D(var)=  (net production)/Cost

D0 = 0.5.*0.9.*0.9.*450.*950.*FIELD.ANNULUS.dishes_in_field; % maximal fantastic installed

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

design_DISH = dish(FIELD.ANNULUS.dishes_per_cluster,FIELD.ANNULUS.clusters_in_field,Slope_Error,reflectivity,FIELD.ANNULUS.Xns,FIELD.ANNULUS.Xew,FIELD.ANNULUS.alpha_angle,field_azimuth_or,field_elevation_or,design_DATA);
design_FIELD = steady_field_annulus(design_DISH.dish_power_pr_sqr_m,FIELD,design_DATA.wind_velocity,design_DATA.ambient_temp,design_DATA.humidity);
D = design_FIELD.POWER.Electric - (FIELD.ANNULUS.dishes_in_field.*design_DISH.consumption+design_FIELD.ANNULUS.all_blowers_consumption+design_FIELD.POWER.consumption.all_pumps_consumption);

if ((D>D_max || D<D_min) || (D0/D>2 || D0/D<0.5)) && (D_max~=0 && D_min~=0)%
    
    go_optimize(4) = 0;
    
else
    
    go_optimize(4) = 1;
    
end

%% LCOE (required LCOE, meaning targeted cost of electricity) in $/KWh
% Maximum or target span LCOE ($/(KW_installed )–) - specified by user
% Maximize the value:  E(var)=   (net production)/Cost


% goal function:  

if Pre_or_Post==1
    
    kd = 0.06; %Real debt interest rate
    n = 20; %Lifetime of the project [years]
    crf = (kd.*(1+kd ).^n)./((1+kd ).^n-1);
    Capex = FIELD.ANNULUS.project_cost;
    Opex = FIELD.ANNULUS.OnM.*n;
    LCOE = (crf.*Capex+Opex)./E_net_annual;

    E = LCOE;

    if (E>E_max || E<E_min) && (E_max~=0 && E_min~=0)

        go_optimize(5) = 0;

    else

        go_optimize(5) = 1;

    end

else
    
    go_optimize(5) = 1;
    
end

%% Number of dishes (a specified number of dishes)
% Maximum or span of dishes in field - specified by user
% Maximize the value:  F(var)=   (net production)/Cost

% goal function:    

F = FIELD.ANNULUS.dishes_in_field;

if (F>F_max || F<F_min) && (F_max~=0 && F_min~=0)
    
    go_optimize(6) = 0;
    
else
    
    go_optimize(6) = 1;
    
end


%% Steam output in t/h pressure and temperature?
% Minimum or span of steam output, annual or design point - specified by user
% Maximize the value:  G(var)=   (net production)/Cost

% goal function:    

G = mass_steam_annual;

if Pre_or_Post==1
    
    if (G>G_max || G<G_min) && (G_max~=0 && G_min~=0)

        go_optimize(7) = 0;

    else

        go_optimize(7) = 1;

    end

else
    
    go_optimize(7) = 1;
    
end

%% Consumption load (the required consumption curve given by a customer. Can be with time of day and price he gets for electricity)
% Maximum acceptable consumption[KWe]/(Gross production [KWe] )  ?[%] (0<x<1) - specified by user
% Maximize the value:  H(var)=   (net production)/Cost

% goal function:    

H = Consumption_annual./(E_net_annual+Consumption_annual);

if Pre_or_Post==1
    
    if (H>H_max || H<H_min) && (H_max~=0 && H_min~=0)

        go_optimize(8) = 0;

    else

        go_optimize(8) = 1;

    end

else
    
    go_optimize(8) = 1;
    
end


%% final verdict:

if sum(go_optimize)==8
    GO = 1;
else
    GO = 0; % NO GO!
end

end

