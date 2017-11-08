function FIELD = steady_field_annulus(dish_power_per_sqr_m,FIELD,wind_velocity,T_ambient,humidity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        % initiation:
        FIELD_0 = FIELD;
        N_dishes_in_field = FIELD.ANNULUS.dishes_in_field;
        T_in_HEX_design = FIELD.HEX.design_inlet_temperature;
        N = FIELD_0.ANNULUS.N;%size(FIELD_0.ANNULUS.SECTION,2); % amount of sections to be calculated
        T_out_receiver_design = 650;
        T01 = FIELD_0.ANNULUS.SECTION{1}.hot_temperature1;
        T02 = FIELD_0.ANNULUS.SECTION{1}.cold_temperature1;
        % initial estimation:
        mdot = FIELD_0.RECEIVER.receiver_peak_efficiency.*0.6.*(dish_power_per_sqr_m./950);        
        FIELD.model_run_singularity = 0;
        dT_in_HEX = 1;
        dT_out_HEX = 1;
        n = 0;
        m = 0;
        
        if  mdot>=0.2%dish_power_per_sqr_m>350
            
            FIELD_0.ANNULUS.SECTION{1}.hot_temperature1 = T_out_receiver_design;
            T_in_receiver = FIELD_0.ANNULUS.SECTION{1}.cold_temperature1;
            
                while abs(dT_out_HEX)>0.01 || abs(dmdot)>0.0001

                        mdot1 = mdot;
                        T_in_receiver1 = T_in_receiver;
                    FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver_design,T_ambient,wind_velocity,dish_power_per_sqr_m,FIELD_0.RECEIVER);
                        mdot = FIELD.RECEIVER.mdot;
                    FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,mdot,wind_velocity,T_ambient);
                        FIELD.ANNULUS.SECTION{N}.cold_temperature2(FIELD.ANNULUS.SECTION{N}.cold_temperature2<FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch) = FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch;
                        T_HEX_in = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
                        T_out_HEX_annulus = FIELD.ANNULUS.SECTION{N}.cold_temperature2;
                        FIELD.ANNULUS.SECTION{N}.cold_temperature2 = T_out_HEX_annulus;
                    FIELD.POWER = PB(N_dishes_in_field,mdot,T_HEX_in,T_ambient,humidity,FIELD_0.HEX);
                        T_out_HEX = FIELD.POWER.T_HEX_outlet;

                    % fit to HEX outlet, by altering T inlet receiver:
                        dT_out_HEX = (T_out_HEX - T_out_HEX_annulus).*mdot;
                        T_in_receiver = T_in_receiver + dT_out_HEX;
                        T_in_receiver(T_in_receiver<T_ambient) = T_ambient;
                        T_in_receiver(T_in_receiver>T_HEX_in) = T_HEX_in;
                        FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T_in_receiver;
                        
                        dmdot = mdot - mdot1;
                        n=n+1;
                        if n>50 && T_in_receiver>T_ambient || (abs(T_in_receiver-T_in_receiver1)<0.0001 && n>10 && abs(dmdot)<0.0001) || isnan(T_HEX_in) || isinf(T_HEX_in)% || T_out_HEX_annulus>FIELD.ANNULUS.SECTION{1}.hot_temperature1
                            
                            break
                            
                        elseif n>5 && mdot<0.2
                            
                            	FIELD_0.ANNULUS.SECTION{1}.hot_temperature1 = T01;
                                FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T02;
                            FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,0,wind_velocity,T_ambient);
                                T_out_receiver = FIELD.ANNULUS.SECTION{1}.hot_temperature1;
                                T_in_receiver = FIELD.ANNULUS.SECTION{1}.cold_temperature1;
                            FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,0,FIELD_0.RECEIVER);
                            FIELD.POWER = PB(N_dishes_in_field,0,0,T_ambient,humidity,FIELD_0.HEX);
                            FIELD.model_run_singularity = 0;
                            break

%                         elseif n>50 && T_in_receiver<T_ambient %&& mdot>0.2
%                     
%                              F_dT_out_HEX = @ (T_in_receiver) steady_dT_out_HEX(T_in_receiver,step_duration,dish_power_per_sqr_m,FIELD_0,wind_velocity,T_ambient,humidity);
%                              options = optimset('Display','iter','TolFun',1e-1);
%                              [T_in_receiver,fval,exitflag,output] = fzero(F_dT_out_HEX,T_ambient,options)
% 
%                             FIELD.model_run_singularity = 1;
%                             FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T_in_receiver;
% 
%                              while abs(dmdot)>0.0001
% 
%                                         mdot1 = mdot;
%                                     FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver_design,T_ambient,dish_power_per_sqr_m,FIELD_0.RECEIVER);
%                                         mdot = FIELD.RECEIVER.mdot;
%                                     FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,mdot,wind_velocity,T_ambient);
%                                         T_HEX_in = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
%                                         T_out_HEX_annulus = FIELD.ANNULUS.SECTION{N}.cold_temperature2;
%                                         T_out_HEX_annulus(T_out_HEX_annulus<FIELD_0.HEX.WS.temperature) = FIELD_0.HEX.WS.temperature;
%                                         FIELD.ANNULUS.SECTION{N}.cold_temperature2 = T_out_HEX_annulus;
%                                     FIELD.POWER = PB(N_dishes_in_field,mdot,T_HEX_in,T_ambient,humidity,FIELD_0.HEX);
%                                         T_out_HEX = FIELD.POWER.T_HEX_outlet;
%                                         dT_out_HEX = (T_out_HEX - T_out_HEX_annulus);%.*abs(T_in_receiver./T_out_HEX)
%                                         dmdot = mdot - mdot1;
%                                         break
% 
%                              end
%                              
%                         elseif n>100
%                             break
                            
                        end
                
        end
%                      
%                     
%                 if n>50 && mdot<0.2
%                     
%                             	FIELD_0.ANNULUS.SECTION{1}.hot_temperature1 = T01;
%                                 FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T02;
%                                 mdot = 0 ;
%                             FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,0,wind_velocity,T_ambient);
%                                 T_out_receiver = FIELD.ANNULUS.SECTION{1}.hot_temperature1;
%                                 T_in_receiver = FIELD.ANNULUS.SECTION{1}.cold_temperature1;
%                             FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,0,FIELD_0.RECEIVER);
%                             FIELD.POWER = PB(N_dishes_in_field,0,0,T_ambient,humidity,FIELD_0.HEX);
%                                 FIELD.model_run_singularity = 0;
%                             
%                 end
%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 T_in_HEX = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
%                 if n>100 && mdot>=0.2 && T_in_HEX < T_in_HEX_design
%                     v=0;
%                     
% 
%                         while abs(dT_out_HEX)>0.01
%                             dmdot=1;
%                             p = 0;
%                             while abs(dmdot)>0.0001
%                                     mdot1 = mdot;
%                                 FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver_design,T_ambient,dish_power_per_sqr_m,FIELD_0.RECEIVER);
%                                     mdot = FIELD.RECEIVER.mdot;
%                                 FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,mdot,wind_velocity,T_ambient);
%                                     T_HEX_in = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
%                                     T_out_HEX_annulus = FIELD.ANNULUS.SECTION{N}.cold_temperature2;
%                                     T_out_HEX_annulus(T_out_HEX_annulus<T_ambient) = FIELD_0.HEX.WS.temperature;
%                                     FIELD.ANNULUS.SECTION{N}.cold_temperature2 = T_out_HEX_annulus;
%                                     FIELD.POWER = PB(N_dishes_in_field,mdot,T_HEX_in,T_ambient,humidity,FIELD_0.HEX);
%                                     T_out_HEX = FIELD.POWER.T_HEX_outlet;
%                                     dmdot = mdot - mdot1;
%                                 p=p+1
%                                 if p>50; break; end
%                             end
%                             
%                         % fit to HEX outlet, by altering T inlet receiver:
%                             dT_out_HEX = (T_out_HEX - T_out_HEX_annulus);
%                             T_in_receiver = T_in_receiver + dT_out_HEX
%                             T_in_receiver(T_in_receiver<T_ambient) = T_ambient;
%                             FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T_in_receiver;
% 
%                             
%                             v=v+1
%                             if v>100; break; end
%                         end
% 
%                 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                T_in_HEX = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
                T_out_receiver = T_out_receiver_design;
            
            if T_in_HEX > T_in_HEX_design && mdot>=0.2
                
                while abs(dT_in_HEX)>0.01 || abs(dT_out_HEX)>0.01 || abs(dmdot)>0.0001 || T_out_receiver>T_out_receiver_design
                    
                        mdot1 = mdot;
                    FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,dish_power_per_sqr_m,FIELD_0.RECEIVER);
                        mdot = FIELD.RECEIVER.mdot;
                    FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,mdot,wind_velocity,T_ambient);
                        FIELD.ANNULUS.SECTION{N}.cold_temperature2(FIELD.ANNULUS.SECTION{N}.cold_temperature2<FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch) = FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch;
                        T_HEX_in = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
                        T_out_HEX_annulus = FIELD.ANNULUS.SECTION{N}.cold_temperature2;
                        FIELD.ANNULUS.SECTION{N}.cold_temperature2 = T_out_HEX_annulus;
                    FIELD.POWER = PB(N_dishes_in_field,mdot,T_HEX_in,T_ambient,humidity,FIELD_0.HEX);
                        T_out_HEX = FIELD.POWER.T_HEX_outlet;

                    % fit to HEX outlet, by altering T inlet receiver:
                        dT_out_HEX = (T_out_HEX - T_out_HEX_annulus);%.*abs(1-T_HEX_in./T_in_HEX_design);
                        T_in_receiver = T_in_receiver + dT_out_HEX;
                        T_in_receiver(T_in_receiver<T_ambient) = T_ambient;
                        T_in_receiver(T_in_receiver>T_HEX_in) = T_HEX_in;
                        FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T_in_receiver;
                    
                    % fit to HEX inlet, by altering T outlet receiver:
                        T_in_HEX = FIELD.ANNULUS.SECTION{N}.hot_temperature2;
                        dT_in_HEX = (T_in_HEX_design - T_in_HEX);
                        T_out_receiver = T_out_receiver + dT_in_HEX;
                        T_out_receiver(T_out_receiver>T_out_receiver_design) = T_out_receiver_design;
                        FIELD_0.ANNULUS.SECTION{1}.hot_temperature1 = T_out_receiver;
                        
                        dmdot = mdot - mdot1;
                        m=m+1;

                        if m>150; break; end

                end
                
            end

            
        else
            
            FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,0,wind_velocity,T_ambient);
                T_out_receiver = FIELD.ANNULUS.SECTION{1}.hot_temperature1;
                T_in_receiver = FIELD.ANNULUS.SECTION{1}.cold_temperature1;
            FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,0,FIELD_0.RECEIVER);
            FIELD.POWER = PB(N_dishes_in_field,0,0,T_ambient,humidity,FIELD_0.HEX);
            FIELD.model_run_singularity = 0;
            
        end
        
        if mdot<0.2 && (n>0 || m>0)
                            
                            	FIELD_0.ANNULUS.SECTION{1}.hot_temperature1 = T01;
                                FIELD_0.ANNULUS.SECTION{1}.cold_temperature1 = T02;
                            FIELD.ANNULUS  = annulus(FIELD_0.ANNULUS,0,wind_velocity,T_ambient);
                                T_out_receiver = FIELD.ANNULUS.SECTION{1}.hot_temperature1;
                                T_in_receiver = FIELD.ANNULUS.SECTION{1}.cold_temperature1;
                            FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,0,FIELD_0.RECEIVER);
                            FIELD.POWER = PB(N_dishes_in_field,0,0,T_ambient,humidity,FIELD_0.HEX);
                            FIELD.model_run_singularity = 0;
        end
        
end

function ANNULUS  = annulus(ANNULUS,mdot,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here        
        
        ANNULUS_0 = ANNULUS;
        system_pressure = ANNULUS.pressure;
        N_dishes_per_cluster = ANNULUS.dishes_per_cluster;
        N_clusters_in_field = ANNULUS.clusters_in_field;
        N = ANNULUS_0.N;%size(ANNULUS_0.SECTION,2); % amount of sections to be calculated
        mdot1 = mdot; % mass flow of one dish/receiver {Kg/sec]
               
        if N_dishes_per_cluster==1
            N_Cluster_Header = 2;
        else
            N_Cluster_Header = 2 + N_dishes_per_cluster; % N max in each cluster
        end
        
        if N_clusters_in_field/4==round(N_clusters_in_field./4)
            center_blower_exist = 0;
        else
            center_blower_exist = 1;
        end
        
        cluster_multi = ANNULUS.SECTION{N_Cluster_Header}.flow_multiplier;
              
        ANNULUS.SECTION{1}.hot_temperature1 = ANNULUS_0.SECTION{1}.hot_temperature1;
        ANNULUS.SECTION{1}.cold_temperature1 = ANNULUS_0.SECTION{1}.cold_temperature1;
        
         
            for i = 1:N
                
                    if mdot > 0.2
                                                
                            multi = ANNULUS.SECTION{i}.flow_multiplier;
                            mdot = mdot1.*multi; % mass flow amount in each section:
                            ANNULUS.SECTION{i} = section(ANNULUS.SECTION{i},mdot,system_pressure,wind_velocity,T_ambient);
                            
                                L = ANNULUS.SECTION{i}.length;
                                if i==1
                                    ANNULUS.SECTION{i}.hot_temperature_loss_by_components = -((30.7./75).*L).*(0.6./mdot1)./1.2;
                                    ANNULUS.SECTION{i}.cold_temperature_loss_by_components = -((22.3./100).*L).*(0.6./mdot1);                                    
                                elseif i>1 && i<=N_Cluster_Header
                                    ANNULUS.SECTION{i}.hot_temperature_loss_by_components = -((27./75).*L).*(0.6./mdot1)./2.3;
                                    ANNULUS.SECTION{i}.cold_temperature_loss_by_components = -((20.7./100).*L).*(0.6./mdot1);         
                                elseif i>N_Cluster_Header
                                    ANNULUS.SECTION{i}.hot_temperature_loss_by_components = -((20.3./75).*L).*(0.6./mdot1)./3.4;
                                    ANNULUS.SECTION{i}.cold_temperature_loss_by_components = -((18.3./100).*L).*(0.6./mdot1);         
                                end  
                                ANNULUS.SECTION{i}.hot_temperature2 = ANNULUS.SECTION{i}.hot_temperature2 + ANNULUS.SECTION{i}.hot_temperature_loss_by_components;
                                ANNULUS.SECTION{i}.cold_temperature1 = ANNULUS.SECTION{i}.cold_temperature1 + ANNULUS.SECTION{i}.cold_temperature_loss_by_components;
                                
                            if i>3 && i<N_Cluster_Header
                                ANNULUS.SECTION{i}.hot_temperature2 = (ANNULUS.SECTION{i}.hot_temperature2.*multi.*air_Cp(ANNULUS.SECTION{i}.hot_temperature2)+Thot_dish.*air_Cp(Thot_dish))./(multi.*air_Cp(ANNULUS.SECTION{i}.hot_temperature2)+1.*air_Cp(Thot_dish)); % new dish heat contribution
                            elseif i<N && i-N_Cluster_Header>1
                                ANNULUS.SECTION{i}.hot_temperature2 = (ANNULUS.SECTION{i}.hot_temperature2.*multi.*air_Cp(ANNULUS.SECTION{i}.hot_temperature2)+(2.*cluster_multi).*Thot_cluster.*air_Cp(Thot_cluster))./(multi.*air_Cp(ANNULUS.SECTION{i}.hot_temperature2)+2.*cluster_multi.*air_Cp(Thot_cluster)); % new cluster heat contribution
                            end
                            
                            if i<N
                                ANNULUS.SECTION{i+1}.hot_temperature1 = ANNULUS.SECTION{i}.hot_temperature2;
                                ANNULUS.SECTION{i+1}.cold_temperature1 = ANNULUS.SECTION{i}.cold_temperature2;
                                ANNULUS.SECTION{i+1}.hot_pressure1 = ANNULUS.SECTION{i}.hot_pressure2;
                                ANNULUS.SECTION{i+1}.cold_pressure1 = ANNULUS.SECTION{i}.cold_pressure2;
                            end                

                            if i==N_Cluster_Header && i>2
                                Thot_cluster = ANNULUS.SECTION{i}.hot_temperature2;
                                Tcold_cluster = ANNULUS.SECTION{i}.cold_temperature2; 
                                Pcold_cluster = ANNULUS.SECTION{i}.cold_pressure2;
                                mdot_cluster = mdot; 
                                blower_air_density = (Pcold_cluster.*100000./0.27805./(Tcold_cluster+273));
                            elseif i==2
                                Thot_dish = ANNULUS.SECTION{i}.hot_temperature2;
                                Thot_cluster = ANNULUS.SECTION{i}.hot_temperature2;
                                Tcold_cluster = ANNULUS.SECTION{i}.cold_temperature2; 
                                Pcold_cluster = ANNULUS.SECTION{i}.cold_pressure2;
                                mdot_cluster = mdot; 
                                blower_air_density = (Pcold_cluster.*100000./0.27805./(Tcold_cluster+273));
                            end
                        
                    else

                            ANNULUS.SECTION{i} = section(ANNULUS.SECTION{i},0,system_pressure,wind_velocity,T_ambient);
                            ANNULUS.SECTION{i}.hot_temperature = 0;
                            ANNULUS.SECTION{i}.cold_temperature = 0;
                            ANNULUS.SECTION{i}.hot_temperature_loss_by_components = 0;
                            ANNULUS.SECTION{i}.cold_temperature_loss_by_components = 0;   
                            dP_section(i) = 0;
                            blower_air_density = 0;
                            mdot_cluster = 0;
                        
                    end
                    
                            % dP output:         
                            dP_hot(i) = abs(ANNULUS.SECTION{i}.hot_pressure1 - ANNULUS.SECTION{i}.hot_pressure2);
                            ANNULUS.SECTION{i}.dP_hot = dP_hot(i);
                            dP_cold(i) = abs(ANNULUS.SECTION{i}.cold_pressure2 - ANNULUS.SECTION{i}.cold_pressure1);
                            ANNULUS.SECTION{i}.dP_cold = dP_cold(i);
                            dP_section(i) = dP_hot(i) + dP_cold(i);
                            ANNULUS.SECTION{i}.dP_total = dP_section(i);                    
                    
            end

            
                        % dP & Blowers consumption calculation:         
                        ANNULUS.dP_HEX = 2000.*(mdot1./0.6).^2./100000;
                        ANNULUS.dP_RECEIVER = 2000.*(mdot1./0.6).^2./100000;                        
                        dP_cluster = (sum(dP_section(1:N_Cluster_Header))+ ANNULUS.dP_RECEIVER).*(2.*N_Cluster_Header./N_dishes_per_cluster);

                        for j = 1:(N-N_Cluster_Header - (1-center_blower_exist))
                            c = ['dP_cluster_' num2str(j)];
                            b = ['cluster_' num2str(j)  '_blower_consumption'];
                            dP = sum(dP_section((N_Cluster_Header+1):(N_Cluster_Header+j))) + ANNULUS.dP_HEX + dP_cluster;
%                             ANNULUS.(genvarname(c)) = dP;
                            blowers_consumption(j) = (dP.*100000).*mdot_cluster./(blower_air_density.*(0.75.*((1+mdot1)./1.85).^2)); % [KW]
                            blowers_consumption = abs(blowers_consumption);
                            blowers_consumption(isnan(blowers_consumption)) = 0;
%                             ANNULUS.(genvarname(b)) = blower_consumption(j);
                            deltaP(j) = dP;
                        end

                    if N_clusters_in_field>1
                        if N_clusters_in_field/4~=round(N_clusters_in_field./4)
                            all_blowers_consumption = 2.*(blowers_consumption(1)+2.*sum(blowers_consumption(2:size(blowers_consumption,2)),2)); % [KW]
                        elseif N_clusters_in_field/4==round(N_clusters_in_field./4)
                            all_blowers_consumption = 4.*sum(blowers_consumption(1:floor(N_clusters_in_field)./4),2); % [KW]
                        else
                            all_blowers_consumption = N_clusters_in_field.*sum(blowers_consumption,2);
                        end
                    else
                        all_blowers_consumption = blowers_consumption(1); % [KW]
                    end
                    

                             ANNULUS.dP_cluster = dP_cluster;
                             ANNULUS.dP_blowers = deltaP;
                             ANNULUS.blowers_consumption = blowers_consumption;
                             ANNULUS.all_blowers_consumption = all_blowers_consumption;
                             
end

function RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,dish_power_per_sqr_m,RECEIVER)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

concentration_ratio = RECEIVER.concentration_ratio;
receiver_peak_efficiency = RECEIVER.receiver_peak_efficiency;

% receiver_efficiency_factor =receiver_peak_efficiency./0.96;%0.8/0.918; % default is 1 (i.e. 100% for max receiver efficiency of 92%)
% dish_power_per_sqr_m = dish_power_per_sqr_m./1000; % [W/m2] to [KW/m2]
% dish_power_pr_sqr_m0 = dish_power_per_sqr_m;
% K0 = 273;
% Tout = T_out_receiver+K0;
% Tin = T_in_receiver+K0;
% Tlow = T_ambient+K0;
% e = 1; % effective receiver emissivity [%]
% sigma = 5.67*10^-11; % Stefan-Boltzman constant [kW/(m^2-K^4)]
% c1= 405.157; % design constant or istantenious = 398.6
% c2 = 200;
% c3 = 76;
% n1 = -0.4;
% n2 = 0.5;
% 
% S = concentration_ratio.*(pi.*0.275.^2); % [m2]
% 
%  if dish_power_pr_sqr_m0>0
%     
%             if dish_power_pr_sqr_m0<0.3
%                 dish_power_per_sqr_m=0.3;
%             end
% 
%             Cp_air_inlet = air_Cp(T_in_receiver);
%             Cp_air_outlet = air_Cp(T_out_receiver);
%             dh = Cp_air_outlet.*Tout - Cp_air_inlet.*Tin;
% 
%         f = 1;
%         m0dot = dish_power_per_sqr_m;
%         m2 = m0dot; m1 = 0.17;
% 
%         while abs(m1-m2)>=0.0000001
%             
%             f = m0dot.*dh./(receiver_efficiency_factor.*(1-sigma.*e.*((Tout+c2.*m0dot.^n1+c3.*dish_power_per_sqr_m.^n2).^4-Tlow.^4)./(dish_power_per_sqr_m.*concentration_ratio)))-dish_power_per_sqr_m.*S;
%             if f>0
%                 m2= m0dot;
%             elseif f<0
%                 m1 = m0dot;
%             elseif f==0 || abs(m0dot)==inf || isnan(m0dot)==1 || isnan(f)==1 || abs(f)==inf
%                 break
%             end
%             m0dot = (m1+m2)./2;
% 
%         end
% 
%             mdot = m0dot;
%             efficiency = receiver_efficiency_factor.*(1-sigma.*e.*((Tout+c2.*mdot.^n1+c3.*dish_power_per_sqr_m.^n2).^4-Tlow.^4)./(dish_power_per_sqr_m.*concentration_ratio)); % receiver efficiency
% 
%             if dish_power_pr_sqr_m0<0.3
%                 mdot = mdot.*dish_power_pr_sqr_m0./dish_power_per_sqr_m;
%             end


if dish_power_per_sqr_m>0

    K0 = 273.2;
    SBC = 5.6704e-8;
    h = natural_convection_coefficient(T_in_receiver,T_ambient,0.394,wind_velocity,1);%[W/m2k]
    emissivity = 0.5;
    Acirc_receiver = 3.18;%[m2]
    Ap_receiver = (pi.*0.275.^2);%[m2]
    Ap_collector = concentration_ratio.*Ap_receiver;%[m2]

    dTwindow = 200; %[C]
    dTambient = 0; %[C]
    dTface = 0; %[C]
    C_rerad = 0.1;

    Qcollected = dish_power_per_sqr_m.*Ap_collector;
    Qrerad  = C_rerad.*Ap_receiver.*emissivity.*SBC.*((T_out_receiver+dTwindow+K0).^4-(T_ambient+dTambient+K0).^4); Qrerad(Qrerad<0) = 0;
    Qconv = h.*Acirc_receiver.*(T_in_receiver+dTface-T_ambient); Qconv(Qconv<0) = 0;

    mdot = receiver_peak_efficiency.*(Qcollected - Qrerad - Qconv)./(1000.*((T_out_receiver - T_in_receiver).*air_Cp((T_in_receiver+T_out_receiver)./2))); %0.8948 (correction factor to Rec 2.1 performance by experiment

    efficiency =  mdot.*(1000.*((T_out_receiver - T_in_receiver).*air_Cp((T_in_receiver+T_out_receiver)./2)))./Qcollected;

else
      
    mdot=0;
    efficiency = 0;
    
end

efficiency(isnan(efficiency)) = 0;
mdot(isnan(mdot)) = 0;
efficiency(isinf(efficiency)) = 0;
mdot(isinf(mdot)) = 0; 


RECEIVER.mdot = mdot;
RECEIVER.efficiency = efficiency;
 
 
end

function SECTION = section(SECTION,mdot,P_system,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        Geometry = SECTION.geometry;
        Phot1 = SECTION.hot_pressure1;
        Pcold1 = SECTION.cold_pressure1;
        Thot1 = SECTION.hot_temperature1;
        Tcold1 = SECTION.cold_temperature1;
        L = SECTION.length;
        
    if mdot > 0.2

                T1 = steady_state_step(Geometry,Thot1,Tcold1,mdot,Phot1,Pcold1,wind_velocity,T_ambient); 
                Thot_dm = T1.hot(1);
                Tcold_dm = T1.cold(1);
                
%___________________________________________________________________________________________________________________________                
% HOT dT____________________________________________________________________________________________________________________                                
                                %||||| Tx = Tamb+exp(x*B+A) |||||%
                                Ahot = log(Thot1-Tcold1);
                                Bhot = log(Thot_dm-Tcold1) - Ahot;
                                Thot2_dL = Tcold1+real(exp(Ahot + Bhot.*(L-2)));
                                
                            if isnan(Thot2_dL) || Thot2_dL>=Thot1
                                dThotm = Thot1 - Thot_dm;%abs(Thot_d2m - Thot2); dT2hotm(dT2hotm>10) = 10;
                                Thot2_dL = Thot1-(L-1).*dThotm;%_not_e
                            end
                            
                                Thot2_dL = real(Thot2_dL);
                                Thot2_dL(isnan(Thot2_dL) || isinf(Thot2_dL || Thot2_dL<mesh_map(1,index1+2)) || Thot2_dL<T_ambient) = Tcold1;                  
%___________________________________________________________________________________________________________________________                
% COLD dT__________________________________________________________________________________________________________________
                                %||||| Tx = Tamb+exp(x*B+A) |||||%
                                
                                if Tcold1-Tcold_dm>0
                                    
                                    Acold = log(Tcold1-T_ambient);
                                    Bcold = log(Tcold_dm-T_ambient) - Acold;
                                    Tcold2_dL = T_ambient+real(exp(Acold - Bcold.*(L-2)));
                                
                                else
                                    
                                    Tcold2_dL = Tcold1 + (L-1).*(Tcold1-Tcold_dm);
                                    
                                end
                                
%                                 dT1coldm = Tcold1 - Tcold_dm;
%                                 Tcold2_dL = Tcold1 + (L-1).*dT1coldm;
                                Tcold2_dL(isnan(Tcold2_dL) || Tcold2_dL<T_ambient) = T_ambient;
%___________________________________________________________________________________________________________________________ 

                                dPhot = section_dP(Thot1,Phot1,L,0,SECTION.geometry(1),SECTION.hot_elbows,mdot);
                                dPcold = section_dP(Tcold1,Pcold1,L,SECTION.geometry(5)+0.0005,SECTION.geometry(6),SECTION.cold_elbows,mdot);                                
                                Phot2 = Phot1 - dPhot;
                                Pcold2 = Pcold1 - dPcold;                 

                T2 = steady_state_step(Geometry,Thot2_dL,Tcold2_dL,mdot,Phot2,Pcold2,wind_velocity,T_ambient); 

                SECTION.hot_pressure2 = Phot2;
                SECTION.cold_pressure2 = Pcold2;
                     
    else
        
                Phot = P_system;
                Pcold = P_system;
                T2 = steady_state_step(Geometry,Thot1,Tcold1,0,Phot,Pcold,wind_velocity,T_ambient);
                SECTION.hot_pressure1 = P_system;
                SECTION.cold_pressure1 = P_system;
                SECTION.hot_pressure2 = P_system;
                SECTION.cold_pressure2 = P_system;
                SECTION.hot_temperature1 = 0;
                SECTION.cold_temperature1 = 0;
                SECTION.hot_temperature2 = 0;
                SECTION.cold_temperature2 = 0;

    end
    
%                 T_mesh1(1,isnan(T_mesh1(1,:))) = T_ambient;
%                 T_mesh2(5,isnan(T_mesh2(5,:))) = T_ambient;
%                 T_mesh1(1,(abs(T_mesh1(1,:))==inf)) = T_ambient;
%                 T_mesh2(5,(abs(T_mesh1(5,:))==inf)) = T_ambient;
%                 SECTION.hot_temperature2 = T_mesh2(5,1);
%                 SECTION.cold_temperature2 = T_mesh2(5,index1);
                SECTION.hot_temperature2 = T2.hot;
                SECTION.cold_temperature2 = T2.cold;  

end

function T = steady_state_step(Geometry,Tinner_air,Touter_air,mdot,Phot,Pcold,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
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

% HTS Parameters:
        if HTS_type>1
            R1 = INPUT.HTS_hot.duct_R;
            R2 = INPUT.HTS_cold.duct_R;
            
            hLo = INPUT.HTS_hot.insulation_occupation;
            hLm = INPUT.HTS_hot.insulation_material;
            hLq = INPUT.HTS_hot.insulation_quality;

            cLt = INPUT.HTS_cold.insulation_thickness;
            cLm = INPUT.HTS_cold.insulation_material;
            cLq = INPUT.HTS_cold.insulation_quality;
        else
            R1 = INPUT.HTS_hot.duct_R;
            R2 = INPUT.HTS_cold.duct_R;
            
            hLt = INPUT.HTS_hot.insulation_thickness;
            hLm = INPUT.HTS_hot.insulation_material;
            hLq = INPUT.HTS_hot.insulation_quality;

            cLt = INPUT.HTS_cold.insulation_thickness;
            cLm = INPUT.HTS_cold.insulation_material;
            cLq = INPUT.HTS_cold.insulation_quality;
        end
        
if mdot>0.2

%% parameters:
    
    K0 = 273.2; % [Deg.K]
    emissivity = 0.3;
    SBC = 5.67E-8; %Stefan Boltzman coefficient [W/(m.^2-K.^4)]  
    
    imperfection_factor = 0.8;
     % Steel (pipe):
    Cp_steel = 500; % [J/Kg.*K]
    k_steel = @(x) 0.015.*x+15.458333333333327; % [W/m.*C]
    rho_steel = 8000; % [Kg/m3]

    % Aluminium (liner):
    Cp_aluminium = 800; % [J/Kg.*K]
    k_aluminium = @(x) 0.00125.*x.^2-0.0875.*x+206.40625; % [W/m.*C]
    rho_aluminium = 2800; % [Kg/m3]

    % Geometry: annulus profile geometry
    r_ii = Geometry(1);
    r_io = Geometry(2);
    r_is1 = Geometry(3);
    r_is2 = Geometry(4);
    r_is3 = Geometry(5);
    r_ic = r_is3+0.0005;
    r_oi = Geometry(6);
    r_oo = Geometry(7);
    r_os1 = Geometry(8);
    r_os2 = Geometry(9);
    r_os3 = Geometry(10);
    r_oc = r_os3+0.0005;

    % Temperature IC: annulus profile

    % [Tinner_air,Tinner_tube_outer_wall,Tinner_1st_insulation_outer_wall,Tinner_2nd_insulation_outer_wall,Tinner_liner_outer_wall,Touter_air,Touter_tube_inner_wall,Touter_tube_outer_wall,Touter_insulation_outer_wall,Touter_liner_outer_wall]

    %% calculation:

                    Tin1 = Tinner_air+K0;
                    Tin2 = Touter_air+K0;
                    Tin3 = T_ambient+K0;
                    
                    Pamb = 1; %[Bar]
                    
                    dT4K1 = (Tin1.*r_io+Tin2.*r_is3)./(2.*(-r_io+r_is3))-K0;
                    dT4K2 = (Tin2.*r_oo+Tin3.*r_os3)./(2.*(-r_oo+r_os3))-K0;
                    
                    h1 = convection_coefficient(Tinner_air,r_ii,0,r_ii,mdot,Phot);
                    h2r_ic = convection_coefficient(Touter_air,r_ic,r_ic,r_oi,mdot,Pcold);
                    h2r_oi = convection_coefficient(Touter_air,r_oi,r_ic,r_oi,mdot,Pcold);
                    h3 = natural_convection_coefficient(T_ambient,T_ambient,r_oc,wind_velocity,Pamb);
                    
                    Rii = 1./(h1.*2.*pi.*r_ii); %(inner tube gas-wall)   
                    Rio = log(r_io./r_ii)./(k_steel(Tinner_air).*2.*pi); % (inner pipe wall)    
                    Ris1 = log(r_is1./r_io)./(insulation_k(hLm(1),hLq(1),dT4K1).*2.*pi); % (inner insulation, high temp)    
                    Ris2 = log(r_is2./r_is1)./(insulation_k(hLm(2),hLq(2),dT4K1).*2.*pi); % (inner insulation, high temp)    
                    Ris3 = log(r_is3./r_is2)./(insulation_k(hLm(3),hLq(3),dT4K1).*2.*pi); % (inner insulation, high temp)    
                    Ric = 1./(h2r_ic.*2.*pi.*r_ic); % (outer tube: r_ic - gas)    
                    Roi = 1./(h2r_oi.*2.*pi.*r_oi); % (outer tube: gas-r_oi)    
                    Roo = log(r_oo./r_oi)./(k_steel(Touter_air).*2.*pi); % (outer pipe wall)    
                    Ros1 = log(r_os1./r_oo)./(insulation_k(cLm(1),cLq(1),dT4K2).*2.*pi); % (inner insulation, high temp)    
                    Ros2 = log(r_os2./r_os1)./(insulation_k(cLm(2),cLq(2),dT4K2).*2.*pi); % (inner insulation, high temp)    
                    Ros3 = log(r_os3./r_os2)./(insulation_k(cLm(3),cLq(3),dT4K2).*2.*pi); % (inner insulation, high temp)    
                    Ril = log(r_ic./r_is3)./(k_aluminium(dT4K1).*2.*pi); % (inner liner)    
                    Roc = log(r_oc./r_os3)./(k_aluminium(dT4K2).*2.*pi); % (outer liner) 
                    Rambient =  1./(h3.*2.*pi.*r_oc);%(outer wal-ambience)   
                    
                    Q1to2 = (Tin1-Tin2)./sum([Rii,Rio,Ris1,Ris2,Ris3,Ril,Ric]); %Conduction from inner to outer tube [W/m]
                    Tr_ii = Tin1-Q1to2.*Rii; % r_ii wall temperature  [degC/m]
                    Tr_io = Tr_ii-Q1to2.*Rio; % r_io wall temperature  [degC/m]
                    Tr_is1 = Tr_io-Q1to2.*Ris1; % r_is1 wall temperature  [degC/m]
                    Tr_is2 = Tr_is1-Q1to2.*Ris2; % r_is2 wall temperature  [degC/m]
                    Tr_is3 = Tr_is2-Q1to2.*Ris3; % r_is2 wall temperature  [degC/m]
                    Tr_ic = Tr_is3-Q1to2.*Ril; % r_ic wall temperature  [degC/m]
                    Tr_oi = Tin2;
                    Tr_oi = T_annular_radiation((Tr_ic-K0),Touter_air,(Tr_oi-K0),r_ic,r_oi,r_oo,k_steel(Tr_oi-K0),h2r_oi)+K0; % dT on outer pipe inner surface due to radiation from inner liner
                    Q2to3 = (Tr_oi-Tin3)./sum([Roi,Roo,Ros1,Ros2,Ros3,Roc,Rambient]); % Conduction from inner to outer tube [W/m]
                    
                    Qconvgas_to3 = 2.*pi.*r_oi.*h2r_oi.*(Tin2-Tr_oi); %Convection from gas to outer tube  [W/m]
                    Q1to2tot = Q1to2-Qconvgas_to3; %Total Heat Transfer from inner to outer gas [W/m]
                    Tr_oo = Tr_oi-Q2to3.*Roo; %r_oo wall temperature  [degC/m]
                    Tr_os1 = Tr_oo-Q2to3.*Ros1;%r_os wall temperature  [degC/m]
                    Tr_os2 = Tr_os1-Q2to3.*Ros2;%r_os wall temperature  [degC/m]
                    Tr_os3 = Tr_os2-Q2to3.*Ros3;%r_os wall temperature  [degC/m]
                    Tr_oc = Tr_os3-Q2to3.*Roc;%r_oc wall temperature  [degC/m]
                    
                    Cp1 = air_Cp(Tinner_air);
                    Cp2 = air_Cp(Touter_air);
                    
                    dTgas1 = (Q1to2tot./(mdot.*Cp1))./1000; %per [m] 
                    dTgas2 = ((Q1to2tot-Q2to3)./(mdot.*Cp2))./1000; %per [m] 
                    %components losses:
                    dTgas1_by_components = (8.5.*0.0176.*(Tin2-Tr_ic)./Ric)./(mdot.*Cp1)./1000;
                    dTgas2_by_components = (3.1.*0.0147.*(Tin3-Tr_oo)./Rambient)./(mdot.*Cp2)./1000;
                    
                    Tin1 = Tin1-dTgas1+dTgas1_by_components; %per [m] 
                    Tin2 = Tin2+dTgas2+dTgas2_by_components; %per [m]
                    T.hot = Tin1-K0; 
                    T.cold = Tin2-K0;
                    T.cold(T.cold<T_ambient) = T_ambient;
                    T.hot(T.hot<T.cold) = T.cold;
                    
else
                    T.hot = 0;
                    T.cold = 0;
end
                    
end

function T = T_annular_radiation(T_inner,T_fluid,T_outer,R_radiative,R_inner,R_outer,k,h)

if sum(isnan([T_inner,T_fluid,T_outer]))>0 || sum(isinf([T_inner,T_fluid,T_outer]))
    
    T = T_outer;
    
else
    
    K0 = 273.2;
    SBC = 5.7603e-8;
    emiss = 0.6;
    r1 = R_inner;
    r2 = R_outer;
    Tf = T_fluid+K0;
    Tb = T_outer+K0;
    Tem = T_inner+K0;

    A = 2.*pi.*r2;

    Qk = @ (T) 2.*pi.*k.*(T-Tb)./log(r2./r1);
    Qh = @ (T) A.*h.*(T-Tf);
    Qrad = @ (T) 2.*pi.*R_radiative.*SBC.*emiss.*((Tem.^4)-(T.^4))+eps;
    F  = @ (T) (Qk(T)+Qh(T))./Qrad(T)-1;

    C1 = 2.*pi.*R_radiative.*SBC.*emiss;
    C2 = A.*h;
    C3 = 2.*pi.*k./log(r2./r1);

    Ca = C1;
    Cb = C2 + C3;
    Cc = -(C3.*Tb + C2.*Tf + C1.*Tem.^4);

    P = [Ca,0,0,Cb,Cc];
    
    if sum(isnan(P))>0 || sum(isinf(P))>0
        
        T = T_fluid;
        
    else
        
        T = roots(P); 
        T = T(imag(T)==0)-K0;

        T = T(T>=Tb);
        T = min(T);

        if numel(T) == 0 || T>=Tem
            T = Tf;
        end

        T(isnan(T) | isinf(abs(T)) | T>Tem | T<Tb) = Tf;

        T = T - K0;
        
    end

end


    if isnan(T) || isinf(T) || T>T_inner || T<T_outer

        T = T_outer;

    end
   
end

%% PARAMETERS FUNCTIONS:

function h = convection_coefficient(T,R_to_h,R_inner,R_outer,mdot,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    P = P.*100000;
    K0 = 273.2; 
    Cp = 1000.*air_Cp(T); % Specific Heat - Cp  [J/(Kg.*K)]
    T = T+K0;
    A = pi.*(R_outer.^2-R_inner.^2);
    Dh = 2.*(R_outer-R_inner); % Hydraulic Diameter [m]
    Rh = R_outer-R_inner;
    epsilon = 5e-5; % pipe wall roughness
    
        % convection coeficient calculation:
        air_conductivity = 0.000000000015207.*T.^3-0.000000048574.*T.^2+0.00010184.*T-0.00039333; % Gas Conductivity - air [W/(m.*K)]
        dens = P./287.05./T; % Density - r [Kg/m3]
        velm = mdot./(A.*dens); % air flow velocity (inner pipe) [m/sec]
        visco = 0.00001827.*(291.15+120)./(T+120).*(T./291.15).^1.5; % Viscosity - m [N.*sec/m2]
        kinevisco = visco./dens; % Kinematic Viscosity - n [m2/sec]
        Re = velm.*2.*Rh./kinevisco; % Reinolds number (inner pipe)
        f = (1./(-1.8.*log10(6.9./Re+(epsilon./(Dh.*3.7)).^1.11))).^2;%0.25./((log10(epsilon./(3.7.*Dh)+5.74./(Re.^0.9))).^2);
        Pr = visco.*Cp./air_conductivity; % Prandtl No. - Pr
        Nuss = ((f./8).*(Re-1000).*Pr)./(1+12.7.*((f./8).^0.5).*(Pr.^(2./3)-1)); % 0.022.*Pr.^0.5.*Re.^0.8; % Nusselt number (inner pipe)
        h = Nuss.*air_conductivity./R_to_h; % h (inner pipe)
        
        h = abs(h);
        
        if h<10 || isnan(h) || isinf(h)
            h=abs(mdot.*10);
        end
        
end

function h = natural_convection_coefficient(T_surface,T_ambient,R_outer,wind_velocity,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    P = P.*100000;
    K0 = 273.2;
    T = T_ambient + K0;


                    air_density = P./287.05./T; % Density - r [Kg/m3]
                    air_conductivity = 0.000000000015207.*T.^3-0.000000048574.*T.^2+0.00010184.*T-0.00039333; % Gas Conductivity - air [W/(m.*K)]
                    if air_conductivity<0.05; air_conductivity = 0.05; end 
                    air_viscosity = 0.00001827.*(291.15+120)./(T+120).*(T./291.15).^1.5; % Viscosity - m [N.*sec/m2]
                    air_kinematic_viscosity = air_viscosity./air_density; % Kinematic Viscosity - n [m2/sec]
                    Cp = 1000.*air_Cp(T_ambient); % Specific Heat - Cp  [KJ/(Kg.*K)]
                    air_diffusivity = air_conductivity./(air_density.*Cp); % Diffusivity - a [m2/sec]
                    Pr = air_viscosity.*Cp./air_conductivity; % Prandtl No. - Pr

                    dTr_oc2amb = T_surface-T_ambient; %Delta T (T r_os - T ambience) - to be changet to a function%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    Ra = 9.81./T.*(dTr_oc2amb.*(2.*R_outer).^3)./(air_kinematic_viscosity.*air_diffusivity); %Rayleigh number: Ra = gbDTL.^3/an
                    Nusselt_number_cc = 0.68+0.67.*Ra.^0.25./(1+(0.492./Pr).^(9./16)).^(4./9); %Nu   (Nat. Conv. Ref: Churchill & Chu): Nu = 0.68+0.67Ra.^0.25/[1+(0.492/Pr).^(9/16)].^(4/9)
                    Nusselt_number_e = 0.683.*Ra.^0.25.*(Pr./(0.861+Pr)).^0.25; %Nu   (Nat. Conv. Ref: Eckert): Nu = 0.683Ra.^0.25[Pr/(0.861+Pr)].^0.25
                    Nusselt_number = abs(max(Nusselt_number_cc,Nusselt_number_e));
                    h_natural = Nusselt_number.*air_conductivity./(2.*R_outer); %natural  convection coefficient
                    Re_forced = wind_velocity.*R_outer.*2./air_kinematic_viscosity; %Reynolds No. based on outer tube diam.
                    Nusselt_number_forced = (0.4.*Re_forced.^0.5+0.06.*Re_forced.^0.667).*Pr.^0.4; %Nu number - forced
                    h_forced = Nusselt_number_forced.*air_conductivity./(2.*R_outer); %forced convection coefficient

                     if wind_velocity>0
                       h=h_forced;
                    else
                       h=h_natural;
                     end
                     
                     h = abs(h);
                     h(h<5) = 5;
                     h(isnan(h)==1) = 10;
                     h(inf==abs(h)) = 10;

end

function dP = section_dP(T,P,L,R_inner,R_outer,elbows,mdot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

epsilon = 5e-5; % pipe roughness [m]

    if mdot<=0
        dP = 0;
    else
        T = T+273; %[Deg.C] to [k]
        P = P.*100000; % [Bar] to [Pa]
        A = pi.*(R_outer.^2-R_inner.^2); % cross section area [m2]
        Dh = 2.*(R_outer-R_inner); % Hydraulic Diameter [m]
        dens = P./287.05./T; % Density - r [Kg/m3]
        velm = mdot./(A.*dens); % air flow velocity (inner pipe) [m/sec]
        visco = 0.00001827.*(291.15+120)./(T+120).*(T./291.15).^1.5; % Viscosity - m [N.*sec/m2]
        kinevisco = visco./dens; % Kinematic Viscosity - n [m2/sec]
        Re = velm.*Dh./kinevisco; % Reinolds number (inner pipe)
        Head = 0.5.*dens.*velm.^2; %pressure head        
        friction_factor = (1./(-1.8.*log10(6.9./Re+(epsilon./(Dh.*3.7)).^1.11))).^2;%0.25./((log10(epsilon./(3.7.*Dh)+5.74./(Re.^0.9))).^2);
        dP_on_Branch = L.*Head.*friction_factor./Dh; %dP per meter [N/m3]
        dP_elbow = elbows.*Head.*0.5; %total elbows P drop per branch [N/m2]

        dP = abs(dP_on_Branch+dP_elbow)./100000; % [Pa] to [Bar]
    end
    
    if dP<0 || isnan(dP)==1 || inf==abs(dP)
        dP=0;
    end
    
end

function Cp = air_Cp(T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T = T+273.2;
% Cp = 0.00000000000022.*T.^4 -0.000000000891606.*T.^3 +0.000001233974292.*T.^2 -0.000480684463337.*T +1.059513120474746;
Cp = (28.09 + 0.1965.*10.^(-2).*T + 0.4799.*10.^(-5).*T.^2 - 1.965.*10.^(-9).*T.^3)./28.97;

if Cp<1 || isnan(Cp)==1 || inf==abs(Cp)
    
    Cp=1.0035;
    
elseif Cp>1.2
    
    Cp = 1.2;

end

end

function a = insulation_alpha(material,quality,T)
% T in [*C]
K0 = 273.2; %[*K]
    if material<1 || quality<0.01 || quality >1
        a = 0;
    else
        a = insulation_k(material,quality,T)./(insulation_Cp(material,quality,T).*insulation_rho(material,quality,T));
    end

end

function k = insulation_k(material,quality,T)
% T in [*C]
K0 = 273.2; %[*K]

    switch material
        
        case 0
            k = inf;
            quality = 1;
        case -1
            k = 0.00125.*T.^2-0.0875.*T+206.40625; % [W/m.*C]
            quality = 1;
        case -2
            k = 0.015.*T+15.458333333333327; % [W/m.*C]
            quality = 1;
        case 1 %'MPS (MicroTherm)'
            k = (-0.0000000233*T.^3+0.0000666.*T.^2-0.000256.*T+26.4)./1000;
        case 2 %'Blanket (MicroTherm)'
            k = (0.00005.*T.^2-0.0063.*T + 30.3)./1000;
        case 3 % 'PyroJell (XT-E)'
            k = (0.000000278.*T.^3-0.0000321.*T.^2+0.0345.*T + 19.9)./1000;
        case 4 % 'Majus (MicroTherm)'
            k = (1.0000e-5).*T+0.0055;
        case 5 %'Rock Wool'
            k = (0.000000059239299.*T.^2+0.000146330816514.*T+0.032362884481313);
        case 6 %'Glass Wool'
            k = (0.000000059239299.*T.^2+0.000146330816514.*T+0.032362884481313)./0.8;
        case 7 % 'Air Cavity (Low Pressure)'
            T = T+K0;
            k = 0.000000000015207.*T.^3-0.000000048574.*T.^2+0.00010184.*T-0.00039333; % Gas Conductivity - air [W/(m.*K)]
            if k>0.05; k = 0.05; end 
    end
    
    k = k./quality;

end

function rho = insulation_rho(material,quality,T)

K0 = 273.2; %[*K]

    switch material
        
        case 0 %
            rho = 0;
        case 1 %'MPS (MicroTherm)'
            rho = 330;
        case 2 %'Blanket (MicroTherm)'
            rho = 330;
        case 3 % 'PyroJell (XT-E)'
            rho = 200;
        case 4 % 'Majus (MicroTherm)'
            rho = 280;
        case 5 %'Rock Wool'
            rho = 120;
        case 6 %'Glass Wool'
            rho = 100;
        case 7 % 'Air Cavity (Low Pressure)'
%             T = T+K0;
%             Pcavity = 0.1; %[Bar]
%             rho = Pcavity.*100000./278.05./T;
            rho = 100;
    end

end

function Cp = insulation_Cp(material,quality,T)

K0 = 273.2; %[*K]

        switch material
        
            case 0 %
                Cp = 0;
            case 1 %'MPS (MicroTherm)'
                Cp = 1080;
            case 2 %'Blanket (MicroTherm)'
                Cp = 1080;
            case 3 % 'PyroJell (XT-E)'
                Cp = 840;
            case 4 % 'Majus (MicroTherm)'
                Cp = 950;
            case 5 %'Rock Wool'
                Cp = 1380;
            case 6 %'Glass Wool'
                Cp = 700;
            case 7 % 'Air Cavity (Low Pressure)'
                Cp = 1000.*air_Cp(T);
        end

end