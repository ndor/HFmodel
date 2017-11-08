function OPT = field_opt_pipes(X,DATA,Slope_Error,reflectivity,system_pressure,HTS_type,FIELD)
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


            PIPES = field_structure_pipes(transient_state);
%                     PIPES.pressure
%                     PIPES.dishes_per_cluster
%                     PIPES.clusters_in_field
%                     PIPES.dishes_in_field
%                     PIPES.dt
%                     PIPES.N
            N = PIPES.N;
            FIELD.PIPES = PIPES;

            
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
                            
                            G_sections_length(j+1) = PIPES.SECTION{j}.length+G_sections_length(j);
                            G_sections_tag{j+1} = {PIPES.SECTION{j}.tag};
   
                         end


                DISH = dish(N_dishes_per_cluster,N_clusters_in_field,Slope_Error,reflectivity,lns,lew,alpha,field_azimuth_or,field_elevation_or,DATA);

                 time_length = size(DISH.dish_power_pr_sqr_m,1);  
                 N = size(PIPES.SECTION,2);
                 FIELD.PIPES = PIPES;

                 % Model calculation:

% %                     F = @ (dish_power_pr_sqr_m,FIELD,wind_velocity,ambient_temp,humidity) steady_field(dish_power_pr_sqr_m,FIELD{:},wind_velocity,ambient_temp,humidity);
% %                     NEF=cell(time_length,1);
% %                     for i=1:time_length; NEF{i} = FIELD; end
% % 
% %                     FIELD = arrayfun(F,DISH.dish_power_pr_sqr_m,NEF,DATA.wind_velocity,DATA.ambient_temp,DATA.humidity,'UniformOutput',false);
% %                     for i=1:time_length
% % 
% %                         Steam(i) = real(FIELD{i}.POWER.mdot_HP+FIELD{i}.POWER.mdot_MP+FIELD{i}.POWER.mdot_LP);
% %                         Consumption(i) = (PIPES.dishes_in_field.*DISH.consumption(i)+FIELD{i}.PIPES.all_blowers_consumption+FIELD{i}.POWER.consumption.all_pumps_consumption); 
% %                         NET_Electric(i) = real(DATA.DNI_commonality(i).*(FIELD{i}.POWER.Electric - (PIPES.dishes_in_field.*DISH.consumption(i)+FIELD{i}.PIPES.all_blowers_consumption+FIELD{i}.POWER.consumption.all_pumps_consumption)));
% % 
% %                     end
                    
                    for i=1:time_length

                        FIELD = steady_field_pipes(DISH.dish_power_pr_sqr_m(i),FIELD,DATA.wind_velocity(i),DATA.ambient_temp(i),DATA.humidity(i));
                        Steam(i) = (DATA.DNI_commonality(i).*(FIELD.POWER.mdot_HP+FIELD.POWER.mdot_MP+FIELD.POWER.mdot_LP));
                        Consumption(i) = DATA.DNI_commonality(i).*(PIPES.dishes_in_field.*DISH.consumption(i)+FIELD.PIPES.all_blowers_consumption+FIELD.POWER.consumption.all_pumps_consumption); 
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
                    
                     OPT = -sum(NET_Electric)./FIELD.PIPES.all_PIPES_cost;


%                     go_optimize = GONOGO(1,FIELD,OPT_USER,Slope_Error,reflectivity,field_azimuth_or,field_elevation_or,E_net_annual,mass_steam_annual,Consumption_annual);
%                     
%                     if go_optimize==0
% 
%                         OPT = 0;
% 
%                     else
%                 
%                         OPT = -sum(NET_Electric)./FIELD.PIPES.all_PIPES_cost;
%                         
%                     end
            
            
%         end

end
    
end

