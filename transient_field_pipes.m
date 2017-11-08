function FIELD = transient_field_pipes(step_duration,dish_power_per_sqr_m,FIELD,wind_velocity,T_ambient,humidity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        % initiation:
        FIELD_00 = FIELD;
        FIELD_0 = FIELD;
        N_dishes_in_field = FIELD.PIPES.dishes_in_field;
        T_in_HEX_design = FIELD.HEX.design_inlet_temperature;
        N = FIELD_0.PIPES.N;%size(FIELD_0.PIPES.SECTION,2); % amount of sections to be calculated
        T_out_receiver_design = 650;
        T01 = FIELD_0.PIPES.SECTION{1}.hot_temperature1;
        T02 = FIELD_0.PIPES.SECTION{1}.cold_temperature1;
        % initial estimation:
        mdot = FIELD_0.RECEIVER.receiver_peak_efficiency.*0.6.*(dish_power_per_sqr_m./950);        
        FIELD.model_run_singularity = 0;
        dT_in_HEX = 1;
        n = 0;
        m = 0;
        
        if  mdot > 0.2%dish_power_per_sqr_m>350
            
            FIELD_0.PIPES.SECTION{1}.hot_temperature1 = T_out_receiver_design;
            T_in_receiver = FIELD_0.PIPES.SECTION{1}.cold_temperature1;
            
            dmdot = 1;
                while abs(dmdot)>0.0001

                        mdot1 = mdot;
                    FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver_design,T_ambient,wind_velocity,dish_power_per_sqr_m,FIELD_0.RECEIVER);
                        mdot = FIELD.RECEIVER.mdot;
                    FIELD.PIPES  = pipes(step_duration,FIELD_0.PIPES,mdot,wind_velocity,T_ambient);
                        FIELD.PIPES.SECTION{N}.cold_temperature2(FIELD.PIPES.SECTION{N}.cold_temperature2<FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch) = FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch;
                        T_HEX_in = FIELD.PIPES.SECTION{N}.hot_temperature2
                    FIELD.POWER = PB(N_dishes_in_field,mdot,T_HEX_in,T_ambient,humidity,FIELD_0.HEX);
                        FIELD_0.PIPES.SECTION{N}.cold_temperature2 = FIELD.POWER.T_HEX_outlet;
                        
                        dmdot = mdot - mdot1;
                        n=n+1;
                
                end
                
                if mdot<=0.2; FIELD = FIELD_00; return; end

                T_in_HEX = FIELD.PIPES.SECTION{N}.hot_temperature2;
                T_out_receiver = T_out_receiver_design;
            
            if T_in_HEX > T_in_HEX_design && mdot>=0.2
                dmdot = 1;
                while abs(dT_in_HEX)>0.01 || abs(dmdot)>0.0001 || T_out_receiver>T_out_receiver_design
                    
                        mdot1 = mdot;
                    FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver_design,T_ambient,wind_velocity,dish_power_per_sqr_m,FIELD_0.RECEIVER);
                        mdot = FIELD.RECEIVER.mdot;
                    FIELD.PIPES  = pipes(step_duration,FIELD_0.PIPES,mdot,wind_velocity,T_ambient);
                        FIELD.PIPES.SECTION{N}.cold_temperature2(FIELD.PIPES.SECTION{N}.cold_temperature2<FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch) = FIELD_0.HEX.WS.temperature+FIELD_0.HEX.HP.pinch;
                        T_HEX_in = FIELD.PIPES.SECTION{N}.hot_temperature2;
                    FIELD.POWER = PB(N_dishes_in_field,mdot,T_HEX_in,T_ambient,humidity,FIELD_0.HEX);
                        FIELD_0.PIPES.SECTION{N}.cold_temperature2 = FIELD.POWER.T_HEX_outlet;
                    
                    % fit to HEX inlet, by altering T outlet receiver:
                        T_in_HEX = FIELD.PIPES.SECTION{N}.hot_temperature2
                        dT_in_HEX = (T_in_HEX_design - T_in_HEX);
                        T_out_receiver = T_out_receiver + dT_in_HEX;
                        T_out_receiver(T_out_receiver>T_out_receiver_design) = T_out_receiver_design;
                        FIELD_0.PIPES.SECTION{1}.hot_temperature1 = T_out_receiver;
                        
                        dmdot = mdot - mdot1;
                        m=m+1;
                        if m>100; break; end

                end
                
            end

            
        else
            
            FIELD.PIPES  = pipes(step_duration,FIELD_0.PIPES,0,wind_velocity,T_ambient);
                T_out_receiver = FIELD.PIPES.SECTION{1}.hot_temperature1;
                T_in_receiver = FIELD.PIPES.SECTION{1}.cold_temperature1;
            FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,0,FIELD_0.RECEIVER);
            FIELD.POWER = PB(N_dishes_in_field,0,0,T_ambient,humidity,FIELD_0.HEX);
            FIELD.model_run_singularity = 0;
            
        end
        
        if mdot <= 0.2 %&& (n>0 || m>0)
                            
                            	FIELD_0.PIPES.SECTION{1}.hot_temperature1 = T01;
                                FIELD_0.PIPES.SECTION{1}.cold_temperature1 = T02;
                            FIELD.PIPES  = pipes(step_duration,FIELD_0.PIPES,0,wind_velocity,T_ambient);
                                T_out_receiver = FIELD.PIPES.SECTION{1}.hot_temperature1;
                                T_in_receiver = FIELD.PIPES.SECTION{1}.cold_temperature1;
                            FIELD.RECEIVER = receiver(T_in_receiver,T_out_receiver,T_ambient,wind_velocity,0,FIELD_0.RECEIVER);
                            FIELD.POWER = PB(N_dishes_in_field,0,0,T_ambient,humidity,FIELD_0.HEX);
                            FIELD.model_run_singularity = 0;
        end
        
end

function PIPES  = pipes(step_duration,PIPES,mdot,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here        
        
        PIPES_0 = PIPES;
        dt = PIPES.dt;
        system_pressure = PIPES.pressure;
        N_dishes_per_cluster = PIPES.dishes_per_cluster;
        N_clusters_in_field = PIPES.clusters_in_field;
        N = PIPES_0.N;%size(PIPES_0.SECTION,2); % amount of sections to be calculated
        mdot1 = mdot; % mass flow of one dish/receiver {Kg/sec]
        steps = round(step_duration./dt);
        
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
        
        cluster_multi = PIPES.SECTION{N_Cluster_Header}.flow_multiplier;

    for s = 1:steps
   
        
 %% HOT        
            for i = 1:1:N 
                
                    if mdot > 0.2
                        
                            PIPES.SECTION{1}.hot_temperature1 = PIPES_0.SECTION{1}.hot_temperature1;
                            multi = PIPES.SECTION{i}.flow_multiplier;
                            mdot = mdot1.*multi; % mass flow amount in each section:
                            PIPES.SECTION{i} = section(PIPES.SECTION{i},1,mdot,system_pressure,wind_velocity,T_ambient);
                            PIPES.SECTION{i}.hot_temperature2(PIPES.SECTION{i}.hot_temperature2<T_ambient) = T_ambient;
                            
                                L = PIPES.SECTION{i}.length;
                                if i==1
                                    PIPES.SECTION{i}.hot_temperature_loss_by_components = -((30.7./75).*L).*(0.6./mdot1);                       
                                elseif i>1 && i<=N_Cluster_Header
                                    PIPES.SECTION{i}.hot_temperature_loss_by_components = -((27./75).*L).*(0.6./mdot1);
                                elseif i>N_Cluster_Header
                                    PIPES.SECTION{i}.hot_temperature_loss_by_components = -((20.3./75).*L).*(0.6./mdot1);
                                end  
                                PIPES.SECTION{i}.hot_temperature2 = PIPES.SECTION{i}.hot_temperature2 + PIPES.SECTION{i}.hot_temperature_loss_by_components;
                                PIPES.SECTION{i}.hot_temperature2(PIPES.SECTION{i}.hot_temperature2<T_ambient) = T_ambient;

                            if i>3 && i<N_Cluster_Header
                                PIPES.SECTION{i}.hot_temperature2 = (PIPES.SECTION{i}.hot_temperature2.*multi.*air_Cp(PIPES.SECTION{i}.hot_temperature2)+Thot_dish.*air_Cp(Thot_dish))./(multi.*air_Cp(PIPES.SECTION{i}.hot_temperature2)+1.*air_Cp(Thot_dish)); % new dish heat contribution
                            elseif i<N && i-N_Cluster_Header>1
                                PIPES.SECTION{i}.hot_temperature2 = (PIPES.SECTION{i}.hot_temperature2.*multi.*air_Cp(PIPES.SECTION{i}.hot_temperature2)+(2.*cluster_multi).*Thot_cluster.*air_Cp(Thot_cluster))./(multi.*air_Cp(PIPES.SECTION{i}.hot_temperature2)+2.*cluster_multi.*air_Cp(Thot_cluster)); % new cluster heat contribution
                            end
                            
                            if i<N
                                PIPES.SECTION{i+1}.hot_temperature1 = PIPES.SECTION{i}.hot_temperature2;
                                PIPES.SECTION{i+1}.hot_pressure1 = PIPES.SECTION{i}.hot_pressure2;
                            end                

                            if i==N_Cluster_Header && i>2
                                Thot_cluster = PIPES.SECTION{i}.hot_temperature2;
                                mdot_cluster = mdot; 
                            elseif i==2
                                Thot_dish = PIPES.SECTION{i}.hot_temperature2;
                                Thot_cluster = PIPES.SECTION{i}.hot_temperature2;
                                mdot_cluster = mdot; 
                            end
                                
                    else

                            PIPES.SECTION{i} = section(PIPES.SECTION{i},1,0,system_pressure,wind_velocity,T_ambient);
                            PIPES.SECTION{i}.hot_temperature_loss_by_components = 0;
                            PIPES.SECTION{i}.cold_temperature_loss_by_components = 0;   
                            dP_section(i) = 0;
                            mdot_cluster = 0;
                            blower_air_density = 0;
                            mdot_cluster = 0;
                        
                    end
                                     
            end
%%
%% COLD
            for i = N:-1:1 
                
                    if mdot > 0.2
                            
                            PIPES.SECTION{N}.cold_temperature2 = PIPES_0.SECTION{N}.cold_temperature2; 
                            multi = PIPES.SECTION{i}.flow_multiplier;
                            mdot = mdot1.*multi; % mass flow amount in each section:
                            PIPES.SECTION{i} = section(PIPES.SECTION{i},2,mdot,system_pressure,wind_velocity,T_ambient);
                            PIPES.SECTION{i}.cold_temperature1(PIPES.SECTION{i}.cold_temperature1<T_ambient) = T_ambient;
                            
                                L = PIPES.SECTION{i}.length;
                                if i==1
                                    PIPES.SECTION{i}.cold_temperature_loss_by_components = -((22.3./100).*L).*(0.6./mdot1);                                    
                                elseif i>1 && i<=N_Cluster_Header
                                    PIPES.SECTION{i}.cold_temperature_loss_by_components = -((20.7./100).*L).*(0.6./mdot1);         
                                elseif i>N_Cluster_Header
                                    PIPES.SECTION{i}.cold_temperature_loss_by_components = -((18.3./100).*L).*(0.6./mdot1);         
                                end  
                                PIPES.SECTION{i}.cold_temperature1 = PIPES.SECTION{i}.cold_temperature1 + PIPES.SECTION{i}.cold_temperature_loss_by_components;
                                PIPES.SECTION{i}.cold_temperature1(PIPES.SECTION{i}.cold_temperature1<T_ambient) = T_ambient;
                                
                            if i>1
                                PIPES.SECTION{i-1}.cold_temperature2 = PIPES.SECTION{i}.cold_temperature1;
                                PIPES.SECTION{i-1}.cold_pressure2 = PIPES.SECTION{i}.cold_pressure1;
                            end                

                            if i==N_Cluster_Header && i>2
                                Tcold_cluster = PIPES.SECTION{i}.cold_temperature2; 
                                Pcold_cluster = PIPES.SECTION{i}.cold_pressure2;
                                blower_air_density = (Pcold_cluster.*100000./0.27805./(Tcold_cluster+273));
                            elseif i==2
                                Tcold_cluster = PIPES.SECTION{i}.cold_temperature2; 
                                Pcold_cluster = PIPES.SECTION{i}.cold_pressure2;
                                blower_air_density = (Pcold_cluster.*100000./0.27805./(Tcold_cluster+273));
                            end
                                
                    else

                            PIPES.SECTION{i} = section(PIPES.SECTION{i},2,0,system_pressure,wind_velocity,T_ambient);
                            PIPES.SECTION{i}.cold_temperature2(PIPES.SECTION{i}.cold_temperature2<T_ambient) = T_ambient;
                            PIPES.SECTION{i}.cold_temperature1(PIPES.SECTION{i}.cold_temperature1<T_ambient) = T_ambient;
                            PIPES.SECTION{i}.hot_temperature2(PIPES.SECTION{i}.hot_temperature2<T_ambient) = T_ambient;
                            PIPES.SECTION{i}.hot_temperature1(PIPES.SECTION{i}.hot_temperature1<T_ambient) = T_ambient;
                            dP_section(i) = 0;
                            blower_air_density = 0;
                            mdot_cluster = 0;
                        
                    end         
                    
            end
%%
                        for i = 1:N
                            % dP output:         
                            dP_hot(i) = abs(PIPES.SECTION{i}.hot_pressure1 - PIPES.SECTION{i}.hot_pressure2);
                            PIPES.SECTION{i}.dP_hot = dP_hot(i);
                            dP_cold(i) = abs(PIPES.SECTION{i}.cold_pressure2 - PIPES.SECTION{i}.cold_pressure1);
                            PIPES.SECTION{i}.dP_cold = dP_cold(i);
                            dP_section(i) = dP_hot(i) + dP_cold(i);
                            PIPES.SECTION{i}.dP_total = dP_section(i);
                        end
                        
                        
    end
                            
%%
                        % dP & Blowers consumption calculation:         
                        PIPES.dP_HEX = 2000.*(mdot1./0.6).^2./100000;
                        PIPES.dP_RECEIVER = 2000.*(mdot1./0.6).^2./100000;                        
                        dP_cluster = (sum(dP_section(1:N_Cluster_Header))+ PIPES.dP_RECEIVER).*(2.*N_Cluster_Header./N_dishes_per_cluster);

                        for j = 1:(N-N_Cluster_Header - (1-center_blower_exist))
                            c = ['dP_cluster_' num2str(j)];
                            b = ['cluster_' num2str(j)  '_blower_consumption'];
                            dP = sum(dP_section((N_Cluster_Header+1):(N_Cluster_Header+j))) + PIPES.dP_HEX + dP_cluster;
%                             PIPES.(genvarname(c)) = dP;
                            blowers_consumption(j) = (dP.*100000).*mdot_cluster./(blower_air_density.*(0.75.*((1+mdot1)./1.85).^2)); % [KW]
%                             PIPES.(genvarname(b)) = blower_consumption(j);
                            deltaP(j) = dP;
                        end
                        
                            blowers_consumption = abs(blowers_consumption);
                            blowers_consumption(isnan(blowers_consumption)) = 0;
                            
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
                    
                    
                             PIPES.dP_cluster = dP_cluster;
                             PIPES.dP_blowers = deltaP;
                             PIPES.blowers_consumption = blowers_consumption;
                             PIPES.all_blowers_consumption = all_blowers_consumption;
                             
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
    emissivity = 1;
    Acirc_receiver = 4.18;%[m2]
    Ap_receiver = (pi.*0.275.^2);%[m2]
    Ap_collector = concentration_ratio.*Ap_receiver;%[m2]

    dTwindow = 200; %[C]
    dTambient = 0; %[C]
    dTface = 0; %[C]
    C_rerad = 0.4;

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

function SECTION = section(SECTION,Tindex,mdot,P_system,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

SECTION.hot_mesh(1,1) = SECTION.hot_temperature1;
SECTION.cold_mesh(1,1) = SECTION.cold_temperature1;
SECTION.hot_mesh(6,1) = SECTION.hot_temperature2;
SECTION.cold_mesh(6,1) = SECTION.cold_temperature2;

T_mesh1 = SECTION.hot_mesh;
T_mesh2 = SECTION.cold_mesh;
L = SECTION.length;
        
    if mdot > 0.2

                if Tindex ==1 % HOT

                        Phot1 = SECTION.hot_pressure1;        
                        Thot1 = SECTION.hot_temperature1;        

                T_mesh1(1,1) = Thot1;
                T_mesh1 = real(Tt(1,T_mesh1,mdot,Phot1,wind_velocity,T_ambient)); 
                Thot_dm = T_mesh1(1,1);
                                
                %___________________________________________________________________________________________________________________________                
                % HOT dT____________________________________________________________________________________________________________________                                
                                                %||||| Tx = Tamb+exp(x*B+A) |||||%
                                                Ahot = log(Thot1-T_ambient);
                                                Bhot = log(Thot_dm-T_ambient) - Ahot;
                                                Thot2_dL = T_ambient + real(exp(Ahot + Bhot.*(L-2)));

                                            if isnan(Thot2_dL) || isinf(Thot2_dL) || Thot2_dL>=Thot1
                                                dThotm = -abs(Thot_dm - Thot1);%abs(Thot_d2m - Thot2); dT2hotm(dT2hotm>10) = 10;
                                                Thot2_dL = Thot1 + (L-2).*dThotm;%_not_e
                                            end

                                                Thot2_dL = real(Thot2_dL);
                                                Thot2_dL(isnan(Thot2_dL) || isinf(Thot2_dL) || Thot2_dL<T_ambient) = T_ambient;               
                                                dPhot = section_dP(Thot1,Phot1,L,0,SECTION.hot_geometry(1),SECTION.hot_elbows,mdot);
                                                Phot2 = Phot1 - dPhot;
                                                T_mesh1(6,1) = Thot2_dL;
                                                T_mesh1 = real(Tt(2,T_mesh1,mdot,Phot2,wind_velocity,T_ambient));
                                                SECTION.hot_pressure2 = Phot2;
                                                SECTION.hot_temperature2 = T_mesh1(6,1);
                                                SECTION.hot_mesh = T_mesh1;
                                                T_mesh2 = SECTION.cold_mesh;
                                                
                else % COLD

                        Pcold2 = SECTION.cold_pressure2;
                        Tcold2 = SECTION.cold_temperature2;

                T_mesh2(6,1) = Tcold2;
                T_mesh2 = real(Tt(2,T_mesh2,mdot,Pcold2,wind_velocity,T_ambient)); 
                Tcold_dm = T_mesh2(6,1);
                
                %___________________________________________________________________________________________________________________________                
                % COLD dT__________________________________________________________________________________________________________________
                                                %||||| Tx = Tamb+exp(x*B+A) |||||%
                                                Acold = log(Tcold2-T_ambient);
                                                Bcold = log(Tcold_dm-T_ambient) - Acold;
                                                Tcold1_dL = T_ambient + real(exp(Acold + Bcold.*(L-2)));

                                            if isnan(Tcold1_dL) || isinf(Tcold1_dL) || Tcold1_dL>=Tcold2
                                                dTcoldm = -abs(Tcold_dm - Tcold2);%abs(Tcold_d2m - Tcold2); dT2coldm(dT2coldm>10) = 10;
                                                Tcold1_dL = Tcold2 + (L-2).*dTcoldm;%_not_e
                                            end

                                                Tcold1_dL = real(Tcold1_dL);
                                                Tcold1_dL(isnan(Tcold1_dL) || isinf(Tcold1_dL) || Tcold1_dL<T_ambient) = T_ambient;            
                                                dPcold = section_dP(Tcold2,Pcold2,L,0,SECTION.cold_geometry(1),SECTION.cold_elbows,mdot);                                
                                                Pcold1 = Pcold2 - dPcold;                 
                                                T_mesh2(1,1) = Tcold1_dL;
                                                T_mesh2 = real(Tt(1,T_mesh2,mdot,Pcold1,wind_velocity,T_ambient)); 
                                                SECTION.cold_pressure1 = Pcold1;                                                
                                                SECTION.cold_temperature1 = T_mesh2(1,1);
                                                SECTION.cold_mesh = T_mesh2;
                                                T_mesh1 = SECTION.hot_mesh;

                end
 %___________________________________________________________________________________________________________________________ 
     
    else
        
                Phot = P_system;
                Pcold = P_system;
                T_mesh1 = Tt(1,T_mesh1,0,P_system,wind_velocity,T_ambient);
                T_mesh1 = Tt(2,T_mesh1,0,P_system,wind_velocity,T_ambient);
                T_mesh2 = Tt(1,T_mesh2,0,P_system,wind_velocity,T_ambient);
                T_mesh2 = Tt(2,T_mesh2,0,P_system,wind_velocity,T_ambient);
                SECTION.hot_pressure1 = P_system;
                SECTION.cold_pressure1 = P_system;
                SECTION.hot_pressure2 = P_system;
                SECTION.cold_pressure2 = P_system;    

    end
    
	SECTION.hot_temperature1 = T_mesh1(1,1);
    SECTION.cold_temperature1 = T_mesh2(1,1);
    SECTION.hot_temperature2 = T_mesh1(6,1);
    SECTION.cold_temperature2 = T_mesh2(6,1);  
    SECTION.hot_mesh = T_mesh1;
    SECTION.cold_mesh = T_mesh2;

end

function mesh_map = Tt(side,mesh_map,mdot,Psys,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if side==1
        Ts = mesh_map(1,:);
    else
        Ts = mesh_map(6,:);
    end
    
%     Ts01 = Ts(1);
%     Ts02 = Ts(index1);
    Rs = mesh_map(2,:);
    drs = mesh_map(3,:);
    material_tag = mesh_map(4,:);
    material_quality = mesh_map(5,:);

%         if mdot > 0
% 
%             Ts0 = dTdt(Ts,Rs,drs,material_tag,dt,FoC0,mdot,Phot,Pcold,wind_velocity,T_ambient);
%             Ti = steady_state_step(Geometry,Ts01,Ts02,mdot,Phot,Pcold,wind_velocity,T_ambient);
% 
%                 if Ts0(1)>Ti.hot 
%                     Ts(1)=Ti.hot;
%                 end
%                 if Ts0(index1)>Ti.cold 
%                     Ts(index1)=Ti.cold;
%                 end
            
    Ts = dTdt(Ts,Rs,drs,material_tag,material_quality,mdot,Psys,wind_velocity,T_ambient);
            
%         else
%             
%             Ts = dTdt(Ts,Rs,drs,material_tag,dt,FoC0,Geometry,0,Phot,Pcold,wind_velocity,T_ambient);
%             
%         end
        
    if side==1
        mesh_map(1,:) = Ts(:);
    else
        mesh_map(6,:) = Ts(:);
    end
    
end

function Ts = dTdt(Ts,Rs,drs,material_tag,material_quality,mdot,Psys,wind_velocity,T_ambient)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    index = find(material_tag==-1);
    FoC0 = Rs(1);
    dt = drs(1);
    
    P_ambient = 1;
    K0 = 273.2; % [Deg.K]
    T_ambient = T_ambient+K0;
    T_dt_old = Ts+K0;
    T_air = Ts(1) +K0;
        
     % Steel (pipe):
    k_steel = @(x) 0.015.*x+15.458333333333327; % [W/m.*C]

    % Aluminium (liner):
    k_aluminium = @(x) 0.00125.*x.^2-0.0875.*x+206.40625; % [W/m.*C]
   
%%%%%%%%%%%%%%%%%%%% Boundary Conditions Set %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    % Inner Duct Boundary:
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     if Rs(i-1)==Rs(i)
%                                         T_dt_old(i) = T_dt_old(i-1);
%                                     else
                                        FoC = FoX(Ts(3),material_tag(4),material_quality(4),drs(4),dt,FoC0);

                                        if mdot <= 0.2
                                            air_dens = Psys.*100000./278.05./T_air;
                                            volume = pi.*Rs(3).^2;
                                            mass = air_dens.*volume;
                                            h = natural_convection_coefficient(Ts(3),Ts(1),Rs(3),0,Psys);
                                        else
                                            h = convection_coefficient(Ts(1),Rs(3),0,Rs(3),mdot,Psys); 
                                        end

                                        k = insulation_k(material_tag(4),material_quality(4),(Ts(3)+Ts(4))./2);
%                                         if material_tag(4)==2; k = k_1st_insulation((Ts(3)+Ts(4))./2); end
%                                         if material_tag(4)==3; k = k_2nd_insulation((Ts(3)+Ts(4))./2); end
%                                         if material_tag(4)==4; k = k_aluminium((Ts(3)+Ts(4))./2); end 
                                        Rk = log((Rs(3)+drs(3))./Rs(3))./(2.*pi.*k_steel(Ts(3)));
                                        Rh = 1./(2.*pi.*Rs(3).*h);
                                        U = 1./(Rk+Rh);
                                        Bi = drs(4).*U./k;
                                        a = Rs(3);
                                        am = Rs(3)-drs(4)./2;
                                        A = am./a;
                                        if am<=0; a = R(3)+drs(4)./2; am = R(3); A = am./a; end
                                        T_dt_new(3) = T_dt_old(3).*(1-2.*FoC.*(1+Bi.*A))+2.*FoC.*(T_dt_old(4)+Bi.*A.*T_air);%T_dt_old(3) + 2.*FoC.*(T_dt_old(4) - T_dt_old(3) -Bi.*A.*(T_dt_old(3) - T_inner_air));
                                        dQ_air = 2.*pi.*Rs(3).*h.*(T_dt_new(3)-T_air);  

                                        if mdot <= 0.2
                                            dT_air = dQ_air./(mass.*air_Cp(T_air-K0).*1000);
                                        else
                                            dT_air = dQ_air./(mdot.*air_Cp(T_air-K0).*1000);
                                        end

                                        T_air_new = T_air + dT_air;
    %                                     T_air_new(T_air_new>Tsteady.hot+K0) = Tsteady.hot+K0;
                                        T_dt_new(1) = T_air_new;
                                        T_dt_old(3) = T_dt_new(3);                        
%                                     end
                                    % Outer Surface - To Environement:
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     if Rs(i-1)==Rs(i)
%                                         T_dt_old(i) = T_dt_old(i-1);
%                                     else

                                        FoC = FoX(Ts(index),material_tag(index-1),material_quality(index-1),drs(index-1),dt,FoC0);

                                        h = natural_convection_coefficient(Ts(index),T_ambient-K0,Rs(index),wind_velocity,P_ambient);
                                        h3 = h;
    %                                     dQ_to_amb = -2.*pi.*Rs(i).*h.*(T_dt_old(i)-T_ambient);
                                        k = insulation_k(material_tag(index-1),material_quality(index-1),(Ts(index)+Ts(index-1))./2);
%                                         if material_tag(index-1)==2;k = k_1st_insulation((Ts(3)+Ts(4))./2); end
%                                         if material_tag(index-1)==3; k = k_2nd_insulation((Ts(3)+Ts(4))./2); end
%                                         if material_tag(index-1)==4; k = k_aluminium((Ts(3)+Ts(4))./2); end
                                        Rk = log((Rs(index)+drs(index))./Rs(index))./(2.*pi.*k_aluminium(Ts(index)));
                                        Rh = 1./(2.*Rs(index).*pi.*h);
                                        U = 1./(Rk+Rh);
                                        Bi = drs(index-1).*U./k;
                                        a = Rs(index);
                                        ap = Rs(index)+drs(index-1)./2;
                                        A = ap./a;
                                        T_dt_new(index) = T_dt_old(index).*(1-2.*FoC.*(1+Bi.*A))+2.*FoC.*(T_dt_old(index-1)+Bi.*A.*T_ambient);
                                        if T_dt_new(index)<T_ambient; T_dt_new(index)=T_ambient; end
                                        T_dt_old(index) = T_dt_new(index);
%                                     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
                                % Heat Conduction in Solids:
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 3:length(Ts)-2

                                    if (i>3 && i<index)
                                        % heat conduction - solids:
                                        if Rs(i-1)==Rs(i) || drs(i)<=0
                                            T_dt_old(i) = T_dt_old(i-1);
                                        else
                                            FoC = FoX(Ts(i),material_tag(i),material_quality(i),drs(i),dt,FoC0);
                                            a1 = (Rs(i)-drs(i)./2);
                                            a2 = (Rs(i)+drs(i)./2);                         
                                            T_dt_new(i) = T_dt_old(i) + (FoC./Rs(i)).*(a1.*T_dt_old(i-1) - 2.*Rs(i).*T_dt_old(i) + a2.*T_dt_old(i+1));
                                            T_dt_old(i) = T_dt_new(i);
                                        end
                                    end
                                                                                  
                                end
                                           
    B=[1,3:index];
        
        T_dt_new(B) = T_dt_new(B)-K0;
        Ts = [T_dt_new,0,0];

end

function T = steady_state_step(Ttag,Geometry,Tin,mdot,Psys,wind_velocity,T_ambient)
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

    % Geometry: PIPES profile geometry
    r_ii = Geometry(1);
    r_io = Geometry(2);
    r_is1 = Geometry(3);
    r_is2 = Geometry(4);
    r_is3 = Geometry(5);
    r_ic = r_is2+0.0005;

    % Temperature IC: PIPES profile

    % [Tinner_air,Tinner_tube_outer_wall,Tinner_1st_insulation_outer_wall,Tinner_2nd_insulation_outer_wall,Tinner_liner_outer_wall,Touter_air,Touter_tube_inner_wall,Touter_tube_outer_wall,Touter_insulation_outer_wall,Touter_liner_outer_wall]

    %% calculation:

    if Ttag<2
        material1 = hLm(1); quality1 = hLq(1);
        material2 = hLm(2); quality2 = hLq(2);
        material3 = hLm(3); quality3 = hLq(3);
    else
        material1 = cLm(1); quality1 = cLq(1);
        material2 = cLm(2); quality2 = cLq(2);
        material3 = cLm(3); quality3 = cLq(3);
    end
    
                    Tin = Tin+K0;
                    T_ambient = T_ambient+K0;
                    
                    Pamb = 1; %[Bar]
                    
                    dT4K1 = (Tin.*r_io+T_ambient.*r_is2)./(2.*(-r_io+r_is2))-K0;
                    
                    h1 = convection_coefficient(Tin-K0,r_ii,0,r_ii,mdot,Psys);
                    h_ambient = natural_convection_coefficient(T_ambient-K0,T_ambient-K0,r_ic,wind_velocity,Pamb);
                    
                    Rii = 1./(h1.*2.*pi.*r_ii); %(inner tube gas-wall)   
                    Rio = log(r_io./r_ii)./(k_steel(Tin-K0).*2.*pi); % (inner pipe wall)    
                    Ris1 = log(r_is1./r_io)./(insulation_k(material1,quality1,dT4K1).*2.*pi); % (inner insulation, high temp)    
                    Ris2 = log(r_is2./r_is1)./(insulation_k(material2,quality2,dT4K1).*2.*pi); % (inner insulation, low temp)    
                    Ris3 = log(r_is3./r_is2)./(insulation_k(material3,quality3,dT4K1).*2.*pi); % (inner insulation, low temp)    
                    Ril = log(r_ic./r_is2)./(k_aluminium(dT4K1).*2.*pi); % (inner liner)    
                    Rambient = 1./(h_ambient.*2.*pi.*r_ic); % (outer tube: r_ic - gas)    
                    
                    dQ = (Tin-T_ambient)./sum([Rii,Rio,Ris1,Ris2,Ris3,Ril,Rambient]); %Conduction from inner to outer tube [W/m]
                    Tr_ii = Tin-dQ.*Rii; % r_ii wall temperature  [degC/m]
                    Tr_io = Tr_ii-dQ.*Rio; % r_io wall temperature  [degC/m]
                    Tr_is1 = Tr_io-dQ.*Ris1; % r_is1 wall temperature  [degC/m]
                    Tr_is2 = Tr_is1-dQ.*Ris2; % r_is2 wall temperature  [degC/m]
                    Tr_is3 = Tr_is2-dQ.*Ris3; % r_is2 wall temperature  [degC/m]
                    Tr_ic = Tr_is3-dQ.*Ril; % r_ic wall temperature  [degC/m]
                    
                    Cp = air_Cp(Tin-K0);
                    
                    dTair = (dQ./(mdot.*Cp))./1000; %per [m] 

                    T = Tin-dTair-K0;
                    
else
                    T = 0;
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

function FoC = FoX(T,tag,quality,dr,dt,FoC0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Outer insulation:
%     Cp_outer_insulation = 840; % [J/Kg.*K]
%     k_outer_insulation = @(x) (0.000000059239299.*x.^2+0.000146330816514.*x+0.032362884481313)./0.7; % [W/m.*C]
%     rho_outer_insulation = 100; % [Kg/m3]

    %a_outer_insulation = @ (T) k_outer_insulation(T)./(Cp_outer_insulation.*rho_outer_insulation);
    
                FoC = FoC0;

                                   switch  tag
                                       
                                       case -2
                                           % Steel (pipe):
                                           Cp_steel = 500; % [J/Kg.*K]
                                           k_steel = @(x) 0.015.*x+15.458333333333327; % [W/m.*C]
                                           rho_steel = 8000; % [Kg/m3]
                                           a_steel = @ (T) k_steel(T)./(Cp_steel.*rho_steel);
                                           Fo = @ (T) a_steel(T).*dt./(dr.^2); 
                                           FoC = Fo(T);
                                       case -1
                                           % Aluminium (liner):
                                           Cp_aluminium = 800; % [J/Kg.*K]
                                           k_aluminium = @(x) 0.00125.*x.^2-0.0875.*x+206.40625; % [W/m.*C]
                                           rho_aluminium = 2800; % [Kg/m3]
                                           a_aluminium = @ (T) k_aluminium(T)./(Cp_aluminium.*rho_aluminium);
                                           Fo = @ (T) a_aluminium(T).*dt./(dr.^2); 
                                           FoC = Fo(T);
                                       case 0
                                           FoC = FoC0.*0.99;
                                       case{1,2,3,4,5,6,7}
                                           FoC = insulation_alpha(tag,quality,T).*dt./(dr.^2); 
                                   end

                                   if FoC>=FoC0 || isnan(FoC) || FoC<=0
                                       FoC = FoC0.*0.99;
                                   end
                                   
end

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
        Nuss = ((f./8).*(Re-1000).*Pr)./(1+12.7.*((f./8).^0.5).*(Pr.^(2/3)-1)); % 0.022.*Pr.^0.5.*Re.^0.8; % Nusselt number (inner pipe)
        h = Nuss.*air_conductivity./R_to_h; % h (inner pipe)
        
        h = abs(h);
        
        if h<0 || isnan(h)==1 || inf==abs(h)
            h=mdot.*100;
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