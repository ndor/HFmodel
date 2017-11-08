function PIPES = field_structure_pipes(transient_state)

global  DATA...
    INPUT...
    field_table...
    dishNreceiver_table...
    ducts_size_table...
    PIPES_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    PIPES_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    HTS_type...
    file_name...
    

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
PIPES.field_size = Xns.*Xew - Xew.*xns;

% all_PIPES_hot_cost = 0;
% all_PIPES_cold_cost = 0;
all_PIPES_cost = 0;             

stainless_steel_price_m3 = @ (x) 35100; %[$/m3]
carbon_steel_price_m3 = @ (x) 4301; %[$/m3]
insulation1_price_m3 = @ (x) 4416; %[$/m3]
insulation2_price_m3 = @ (x) 4; %[$/m3]
aluminium_price_m3 = @ (x) 5994; %[$/m3]

%% minimum dish spacing alert:
    num_field_table = cell2mat({field_table{:,2}})';

            alpha = num_field_table(5);
            minimum_Lns = num_field_table(6);
            P1 = [0.5,(-sind(alpha).*minimum_Lns./2),((minimum_Lns./2).^2-39.^2)];
            r1 = roots(P1);
            r1 = r1(r1>0);
%             P2 = [(1+tand(alpha)),(39.*tand(alpha)),((minimum_Lns./2).^2-minimum_Lns.^2)];
%             r2 = roots(P2);
%             r2 = r2(r2>0);
%             minimum_Lew = 2.*max([r1,r2]);
            minimum_Lew = r1;

lew(lew<minimum_Lew) = minimum_Lew;
lns(lns<minimum_Lns) = minimum_Lns;

%% Field Geometry Determination:
if N_dishes_per_cluster==1 && N_clusters_in_field==1
   N_Cluster_Header = 0;
   N_Field_Header = 0;
   N = 2+1;%+1 is for HEX_distance_from_field_center
elseif N_dishes_per_cluster==1 && N_clusters_in_field>1
   N_Cluster_Header = 2;
   N_Field_Header = ceil((N_clusters_in_field./2 - 1)./2);
   N = N_Cluster_Header + N_Field_Header+1; %+1 is for HEX_distance_from_field_center
else
   N_Cluster_Header = 2 + N_dishes_per_cluster;
   N_Field_Header = ceil((N_clusters_in_field./2 - 1)./2);
   N = N_Cluster_Header + N_Field_Header+1; %+1 is for HEX_distance_from_field_center
end

%% Sections Assignement:
% PIPES{i,j} => PIPES { [N1,N2,...Nn] , [Geometry,Elbows,T_map] }
total_number_of_elenerts = 0;
    
     % Steel (pipe):
    Cp_steel = 500; % [J/Kg.*K]
    k_steel = @(x) 0.015.*x+15.458333333333327; % [W/m.*C]
    rho_steel = 8000; % [Kg/m3]

    % Aluminium (liner):
    Cp_aluminium = 800; % [J/Kg.*K]
    k_aluminium = @(x) 0.00125.*x.^2-0.0875.*x+206.40625; % [W/m.*C]
    rho_aluminium = 2800; % [Kg/m3]
    
    % alphas:
    a_steel = @ (T) k_steel(T)./(Cp_steel.*rho_steel);
    a_aluminium = @ (T) k_aluminium(T)./(Cp_aluminium.*rho_aluminium);
    
    
    
    T = 0;
    
    pdr(1) = R1;
    for i = 1:3
        if HTS_type>1
            pdr(i+1) = pdr(i)+(R2-R1).*hLo;
        else
            pdr(i+1) = pdr(i)+hLt(i);
        end
    end

%     % Fo post definition: 
        mesh_rank = 2;
    for i = 1:3
        a_insulation(i) = insulation_alpha(hLm(i),hLq(i),T);
        k_insulation(i) = insulation_k(hLm(i),hLq(i),T);
        dr = pdr(i+1)-pdr(i);
        dr1 = dr./mesh_rank;
        am = R1+0.003-dr1./2;
        a = R1+0.003;
        aa = am./a;
        Rk = log((R1+0.003)./R1)/(2.*pi*k_steel(0));
        Rh = 1./(2.*pi.*1);
        U = 1./(Rh+Rk);
        Bi = U.*dr1./k_insulation(i);
        FoC = 0.1./(2.*(1+Bi.*aa));  
%         FoC(FoC>0.45) = 0.45;
        dt(i) = floor(FoC.*dr1.^2./a_insulation(i));
    end

    dt = dt(~isnan(dt));
    dt = dt(~isinf(dt));
    dt = dt(dt>0)
    dt = min(dt);

%     FoC = FoC(~isnan(FoC));
%     FoC = FoC(~isinf(FoC));
%     FoC = FoC(FoC>0);
%     FoC = min(FoC);
%     FoC(FoC>0.5) = 0.1;
    
%         FoC = a_1st_inner_insulation(0).*dt./(dr1.^2);
        if DATA.step_duration>1 && dt>DATA.step_duration
            dt = DATA.step_duration;
        elseif DATA.step_duration>1
            dt  = ceil(DATA.step_duration./ceil(DATA.step_duration./dt));
        else
            dt = 1;
        end
        

%     % HOT pipe:
%         pipe_thickness = ceil(1000.*((hot_r.*(1+(system_pressure.*14.6./(15500.*0.8+system_pressure.*14.6.*0.4)))-0.00098)))./1000;
%         dr = (min(hot_insulation1,hot_insulation2)-pipe_thickness)./mesh_rank;
%         a = hot_r+pipe_thickness;
%         am = a-dr./2;
%         aa = am./a;
%         Rk = log((a)./hot_r)/(2.*pi*k_steel(0));
%         Rh = 1./(2.*pi.*1);
%         U = 1./(Rh+Rk);
%         Bi = U.*dr./k_1st_inner_insulation(0);
%         FoC = 0.1./(2.*(1+Bi.*aa));      
%         dt = floor(FoC.*dr.^2./a_1st_inner_insulation(0));
%         dt_hot(dt>DATA.step_duration) = DATA.step_duration
% 
% 
%     % COLD pipe:
%         pipe_thickness = ceil(1000.*((cold_r.*(1+(system_pressure.*14.6./(15500.*0.8+system_pressure.*14.6.*0.4)))-0.00098)))./1000;
%         dr = (min(cold_insulation1,cold_insulation2)-pipe_thickness)./mesh_rank;
%         a = cold_r+pipe_thickness;
%         am = a-dr./2;
%         aa = am./a;
%         Rk = log((a)./cold_r)/(2.*pi*k_steel(0));
%         Rh = 1./(2.*pi.*1);
%         U = 1./(Rh+Rk);
%         Bi = U.*dr./k_1st_inner_insulation(0);
%         FoC = 0.1./(2.*(1+Bi.*aa));      
%         dt = floor(FoC.*dr.^2./a_1st_inner_insulation(0));
%         dt_cold(dt>DATA.step_duration) = DATA.step_duration
%         
%         dt = min([dt_hot,dt_cold])
        

         A0 = 1;       
        
        for i = 1:N
         

           if i==1
                    L = 42.5;
          elseif i==2
                    L = lew./4-7.7;
          elseif i<=N_Cluster_Header && i>2 
              
              if i./2~=round(i./2)
                  L = (lns./2).*(1-tand(alpha));
              else
                  L = (lns./2).*(1+tand(alpha));
              end

          elseif i>N_Cluster_Header && i<N
                    L = lew./cosd(alpha);
                    if N_clusters_in_field==4 %ceil((N_clusters_in_field./2 - 1)./2)~=(N_clusters_in_field./2 - 1)./2
                        L=L./2;
                    end
           elseif i==N
               L = HEX_distance_from_field_center;
           end 
           L(L<2) = 2;
           
        %% Geometry:

            A = 1;
        
            if i<=2 % N1=N2=N3/2
                
                flow_multiplier = 1;
                tag = ['N' num2str(i)];
                
            else % N4....N_field
                
                if i>2 && i<=N_Cluster_Header %N4....Ncluster
                    flow_multiplier = i-2;              
                    A = sqrt(flow_multiplier);
                    tag = ['N' num2str(i)];
                    
                else %Ncluster+1.....Nfield
                
                    if i==N % tag_number
                        flow_multiplier = N_dishes_per_cluster.*N_clusters_in_field;
                        tag = ['X'];
                        A = sqrt(flow_multiplier);
                    else
                        flow_multiplier = 2.*N_dishes_per_cluster.*(i-N_Cluster_Header);
                        tag = ['H' num2str(i-N_Cluster_Header)];
                        A = sqrt(flow_multiplier);
                    end        

                end
                
            end

            % HOT pipe Geometry:
            r1 = ceil(10000.*(A.*R1))./10000;
            r2 = ceil(10000.*((r1.*(1+(system_pressure.*14.6./(1550.*0.5+system_pressure.*14.6.*0.7))-0.00098))))./10000;
            r3 = r2+hLt(1);
            r4 = r3+hLt(2);
            r5 = r4+hLt(3);
            
            hot_Geometry_matrix = [r1,r2,r3,r4,r5];
            
            N_bypasses = floor(L./9); N_bypasses(N_bypasses<1) = 1;
            hot_bypass_length = N_bypasses.*3.5; % cold bypass length [m]
            hot_bypass_volume = hot_bypass_length.*pi.*(r2.^2-r1.^2); % [m3]

            % COLD pipe Geometry:
            r1 = ceil(10000.*(A.*R2))./10000;
            r2 = ceil(10000.*((r1.*(1+(system_pressure.*14.6./(1550.*0.5+system_pressure.*14.6.*0.7))-0.00098))))./10000;
            r3 = r2+cLt(1);
            r4 = r3+cLt(2);
            r5 = r4+cLt(3);
            
            cold_Geometry_matrix = [r1,r2,r3,r4,r5];
            
            cold_bypass_length = N_bypasses.*3.5; % cold bypass length [m]
            cold_bypass_volume = cold_bypass_length.*pi.*(r2.^2-r1.^2); % [m3]
            
            bypass_length = max(hot_bypass_length,cold_bypass_length); % [m]
            
                    
        %% Elbows:
            if i==1
                
                hot_elbows = 15.2;
                cold_elbows = 138.96;
                
            elseif i==2
                
                hot_elbows = 10.1;
                cold_elbows = 32.4;
                
            elseif (i>2 && i<N_Cluster_Header) || (i<N && i>N_Cluster_Header)
                
                hot_elbows = 3;
                cold_elbows = 3;
                
            elseif i==N_Cluster_Header
                
                hot_elbows = 8;
                cold_elbows = 13;
            
            elseif i==N
                
                hot_elbows = 5.5;
                cold_elbows = 8.5;
                
            end                  
                
            % Cost:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                stainless_steel_cost = stainless_steel_price_m3(N_clusters_in_field.*N_dishes_per_cluster).*(L.*pi.*(hot_Geometry_matrix(2).^2 - hot_Geometry_matrix(1).^2 + cold_Geometry_matrix(2).^2 - cold_Geometry_matrix(1).^2) + cold_bypass_volume);
                carbon_steel_cost = 0;%carbon_steel_price_m3(N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(Geometry_matrix(7).^2 - Geometry_matrix(6).^2);
                hot_insulation1_cost = insulation_price_m3(hLm(1),N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(hot_Geometry_matrix(3).^2 - hot_Geometry_matrix(2).^2);
                hot_insulation2_cost = insulation_price_m3(hLm(2),N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(hot_Geometry_matrix(4).^2 - hot_Geometry_matrix(3).^2);
                hot_insulation3_cost = insulation_price_m3(hLm(3),N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(hot_Geometry_matrix(5).^2 - hot_Geometry_matrix(4).^2);
                cold_insulation1_cost = insulation_price_m3(cLm(1),N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(cold_Geometry_matrix(3).^2 - cold_Geometry_matrix(2).^2);
                cold_insulation2_cost = insulation_price_m3(cLm(2),N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(cold_Geometry_matrix(4).^2 - cold_Geometry_matrix(3).^2);
                cold_insulation3_cost = insulation_price_m3(cLm(3),N_clusters_in_field.*N_dishes_per_cluster).*L.*pi.*(cold_Geometry_matrix(5).^2 - cold_Geometry_matrix(4).^2);
                aluminium_cost = aluminium_price_m3(N_clusters_in_field.*N_dishes_per_cluster).*2.*L.*pi.*((hot_Geometry_matrix(5)+0.0005).^2 - hot_Geometry_matrix(5).^2 + (cold_Geometry_matrix(5)+0.0005).^2 - cold_Geometry_matrix(5).^2);
                
                PIPES.SECTION{i}.tag = tag;
                PIPES.SECTION{i}.flow_multiplier = flow_multiplier;
                PIPES.SECTION{i}.hot_geometry = hot_Geometry_matrix;
                PIPES.SECTION{i}.cold_geometry = cold_Geometry_matrix;
                PIPES.SECTION{i}.hot_elbows = hot_elbows;
                PIPES.SECTION{i}.cold_elbows = cold_elbows;
                PIPES.SECTION{i}.length = L;
                PIPES.SECTION{i}.cold_bypasses = N_bypasses;
                PIPES.SECTION{i}.cold_bypass_length = cold_bypass_length;
%                 PIPES.SECTION{i}.cold_bypass_thickness = cold_bypass_thickness;
                T_ambient = DATA.ambient_temp(1);
                if transient_state == 1
                    mesh = meshing(dt,mesh_rank,hot_Geometry_matrix,cold_Geometry_matrix,T_ambient);
                    PIPES.SECTION{i}.hot_mesh = mesh.hot_mesh_map;
                    PIPES.SECTION{i}.cold_mesh = mesh.cold_mesh_map;
                    PIPES.SECTION{i}.number_of_elenerts = 2.*(numel(mesh.hot_mesh_map(2,:))+numel(mesh.cold_mesh_map(2,:)) - 7);
                    total_number_of_elenerts = total_number_of_elenerts + PIPES.SECTION{i}.number_of_elenerts;
                end
                
                PIPES.SECTION{i}.hot_temperature1 = T_ambient;
                PIPES.SECTION{i}.cold_temperature1 = T_ambient;
                PIPES.SECTION{i}.hot_temperature2 = T_ambient;
                PIPES.SECTION{i}.cold_temperature2 = T_ambient;
                PIPES.SECTION{i}.hot_pressure1 = system_pressure;
                PIPES.SECTION{i}.cold_pressure1 = system_pressure;
                PIPES.SECTION{i}.hot_pressure2 = system_pressure;
                PIPES.SECTION{i}.cold_pressure2 = system_pressure;
                PIPES.SECTION{i}.stainless_steel_cost = stainless_steel_cost;
                PIPES.SECTION{i}.carbon_steel_cost = carbon_steel_cost;
                PIPES.SECTION{i}.hot_insulation1_cost = hot_insulation1_cost;
                PIPES.SECTION{i}.hot_insulation2_cost = hot_insulation2_cost;
                PIPES.SECTION{i}.hot_insulation3_cost = hot_insulation3_cost;
                PIPES.SECTION{i}.aluminium_cost = aluminium_cost;
                PIPES.SECTION{i}.cold_insulation1_cost = cold_insulation1_cost;
                PIPES.SECTION{i}.cold_insulation2_cost = cold_insulation2_cost;
                PIPES.SECTION{i}.cold_insulation3_cost = cold_insulation3_cost;
                PIPES.SECTION{i}.section_cost = stainless_steel_cost + carbon_steel_cost + hot_insulation1_cost + hot_insulation2_cost + hot_insulation3_cost + cold_insulation1_cost + cold_insulation2_cost + cold_insulation3_cost + aluminium_cost;
                
                
                
                if i<=2
                    PIPES.SECTION{i}.sections_field_cost = N_dishes_per_cluster.*N_clusters_in_field.*PIPES.SECTION{i}.section_cost;
                elseif i>2 && i<=N_Cluster_Header
                    PIPES.SECTION{i}.sections_field_cost = N_clusters_in_field.*PIPES.SECTION{i}.section_cost;
                elseif i>N_Cluster_Header && i<N
                    if N_Cluster_Header==1 && N_clusters_in_field==4
                        times = 4;
                    else
                        times = 2;
                    end
                    PIPES.SECTION{i}.sections_field_cost = times.*PIPES.SECTION{i}.section_cost;
                elseif i==N
                    PIPES.SECTION{i}.sections_field_cost = PIPES.SECTION{i}.section_cost;
                end
                all_PIPES_cost = all_PIPES_cost + PIPES.SECTION{i}.sections_field_cost;
        end
        
        PIPES.single_dish_cost = dish_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.all_dishes_cost = PIPES.single_dish_cost.*N_clusters_in_field.*N_dishes_per_cluster;
        PIPES.steam_cycle_cost = steam_cycle_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.power_unit_cost = power_unit_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.pumps_cost = pumps_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.water_treatment_system_cost = water_treatment_system_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.air_system_cost = air_system_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.blowers_cost = blowers_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.Instruments_cost = instruments_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.control_cost = control_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.other_essentials_cost = other_essentials_cost(N_clusters_in_field.*N_dishes_per_cluster);
        PIPES.all_pipes_cost = all_PIPES_cost;
        PIPES.project_cost = PIPES.all_pipes_cost + PIPES.other_essentials_cost + PIPES.control_cost + PIPES.Instruments_cost ...
            + PIPES.blowers_cost + PIPES.air_system_cost + PIPES.water_treatment_system_cost + PIPES.pumps_cost + PIPES.power_unit_cost ...
            + PIPES.steam_cycle_cost + PIPES.all_dishes_cost;
        
        PIPES.pressure = system_pressure;
        PIPES.Xns = lns;
        PIPES.Xew = lew;
        PIPES.alpha_angle = alpha;
        PIPES.dishes_per_cluster = N_dishes_per_cluster;
        PIPES.clusters_in_field = N_clusters_in_field;
        PIPES.dishes_in_field = N_clusters_in_field.*N_dishes_per_cluster;
        if transient_state == 1
            %PIPES.SECTION{:}.FoC = FoC;
            PIPES.dt = dt;%mesh.hot_mesh_map(3,1);
            PIPES.total_number_of_elenerts = total_number_of_elenerts;
        end
        PIPES.N = N;
        
        %%
        OnM = 333; % annual
        PIPES.OnM = OnM;  % annual
        
end

function mesh_map = meshing(dt,mesh_rank,hot_Geometry,cold_Geometry,T_ambient)
%MESH Summary of this function goes here
%   Detailed explanation goes here
global  DATA...
    INPUT...
    field_table...
    dishNreceiver_table...
    ducts_size_table...
    PIPES_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    PIPES_Xb_table... 
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
        
%% parameters:

    
    imperfection_factor = 0.1;
    
     % Steel (pipe):
    Cp_steel = 500; % [J/Kg.*K]
    k_steel = @(x) 0.015.*x+15.458333333333327; % [W/m.*C]
    rho_steel = 8000; % [Kg/m3]

    % Aluminium (liner):
    Cp_aluminium = 800; % [J/Kg.*K]
    k_aluminium = @(x) 0.00125.*x.^2-0.0875.*x+206.40625; % [W/m.*C]
    rho_aluminium = 2800; % [Kg/m3]
   
    % alphas:
    a_steel = @ (T) k_steel(T)./(Cp_steel.*rho_steel);
    a_aluminium = @ (T) k_aluminium(T)./(Cp_aluminium.*rho_aluminium);
        
    hot_Geometry = [hot_Geometry(1:5),hot_Geometry(5)+0.0005];
    cold_Geometry = [cold_Geometry(1:5),cold_Geometry(5)+0.0005];

    pdr(1) = R1;
    for i = 1:3
        if HTS_type>1
            pdr(i+1) = pdr(i)+(R2-R1).*hLo;
        else
            pdr(i+1) = pdr(i)+hLt(i);
        end
    end
    
%     % Fo post definition: 
        dr1 = (hot_Geometry(3) - hot_Geometry(2))./mesh_rank;
        dr2 = (hot_Geometry(4) - hot_Geometry(3))./mesh_rank;
        dr3 = (hot_Geometry(5) - hot_Geometry(4))./mesh_rank;
        dr4 = (cold_Geometry(3) - cold_Geometry(2))./mesh_rank;
        dr5 = (cold_Geometry(4) - cold_Geometry(3))./mesh_rank;
        dr6 = (cold_Geometry(5) - cold_Geometry(4))./mesh_rank;
        drHOT = [dr1,dr2,dr3];
        drCOLD = [dr4,dr5,dr6];
        dr = drHOT;%[drHOT,drCOLD];
        dr(dr<=0) = inf;
        dr1 = min(dr);
        
%     % Fo post definition: 
        index = find(dr==dr1);
        FoC = insulation_alpha(hLm(index),hLq(index),0).*dt./(dr1.^2);
%         FoC(FoC>0.5) = 0.1;

    % node size (mesh resolution):
    
    % HOT:
    hot_dr(1) = hot_Geometry(2) - hot_Geometry(1);%sqrt(a_steel(0).*dt./(FoC)); % dr_steel1;
    hot_dr(2) = sqrt(insulation_alpha(hLm(1),hLq(1),0).*dt./(FoC)); 
    hot_dr(3) = sqrt(insulation_alpha(hLm(2),hLq(2),0).*dt./(FoC)); 
    hot_dr(4) = sqrt(insulation_alpha(hLm(3),hLq(3),0).*dt./(FoC)); 
%     hot_dr(5) = 0.0005;%sqrt(a_aluminium(0).*dt./(FoC)); % dr_aluminium1;

    % COLD:
    cold_dr(1) = cold_Geometry(2) - cold_Geometry(1);%sqrt(a_steel(0).*dt./(FoC)); % dr_steel1;
    cold_dr(2) = sqrt(insulation_alpha(cLm(1),cLq(1),0).*dt./(FoC)); 
    cold_dr(3) = sqrt(insulation_alpha(cLm(2),cLq(2),0).*dt./(FoC));
    cold_dr(4) = sqrt(insulation_alpha(cLm(3),cLq(3),0).*dt./(FoC)); 
%     cold_dr(5) = 0.0005;%sqrt(a_aluminium(0).*dt./(FoC)); % dr_aluminium1;
    
    
hot_cyl(1) = hot_Geometry(2) - hot_Geometry(1);
hot_cyl(2) = hot_Geometry(3) - hot_Geometry(2); 
hot_cyl(3) = hot_Geometry(4) - hot_Geometry(3);
hot_cyl(4) = hot_Geometry(5) - hot_Geometry(4);

cold_cyl(1) = cold_Geometry(2) - cold_Geometry(1);
cold_cyl(2) = cold_Geometry(3) - cold_Geometry(2); 
cold_cyl(3) = cold_Geometry(4) - cold_Geometry(3);
cold_cyl(4) = cold_Geometry(5) - cold_Geometry(4);


hot_Nrs(1) = 1;%ceil(cyl(1)./dr(1));
hot_Nrs(2) = ceil(hot_cyl(2)./hot_dr(2));
hot_Nrs(3) = ceil(hot_cyl(3)./hot_dr(3));
hot_Nrs(4) = ceil(hot_cyl(4)./hot_dr(4));
hot_Nrs(5) = 1;

cold_Nrs(1) = 1;%ceil(cyl(1)./dr(1));
cold_Nrs(2) = ceil(cold_cyl(2)./cold_dr(2));
cold_Nrs(3) = ceil(cold_cyl(3)./cold_dr(3));
cold_Nrs(4) = ceil(cold_cyl(4)./cold_dr(4));
cold_Nrs(5) = 1;

hot_drs(1) = hot_cyl(1)./hot_Nrs(1);
hot_material_tag(1) = -2;
hot_drs(2) = hot_cyl(2)./hot_Nrs(2); 
hot_material_tag(2) = hLm(1);
hot_drs(3) = hot_cyl(3)./hot_Nrs(3);
hot_material_tag(3) = hLm(2);
hot_drs(4) = hot_cyl(4)./hot_Nrs(4);
hot_material_tag(4) = hLm(3);
% hot_drs(5) = hot_dr(5); 
hot_material_tag(5) = -1;

cold_drs(1) = cold_cyl(1)./cold_Nrs(1); 
cold_material_tag(1) = -2;
cold_drs(2) = cold_cyl(2)./cold_Nrs(2); 
cold_material_tag(2) = cLm(1);
cold_drs(3) = cold_cyl(3)./cold_Nrs(3);
cold_material_tag(3) = cLm(2);
cold_drs(4) = cold_cyl(4)./cold_Nrs(4);
cold_material_tag(4) = cLm(3);
% cold_drs(5) = cold_dr(5); 
cold_material_tag(5) = -1;

hot_drs(isnan(hot_drs)) = 0;
hot_drs(isinf(hot_drs)) = 0;
cold_drs(isnan(cold_drs)) = 0;
cold_drs(isinf(cold_drs)) = 0;
hot_Nrs(isnan(hot_Nrs)) = 0;
hot_Nrs(isinf(hot_Nrs)) = 0;
cold_Nrs(isnan(cold_Nrs)) = 0;
cold_Nrs(isinf(cold_Nrs)) = 0;

% HOT:
hot_Rs(1) = hot_Geometry(1);
hot_dr(1) = hot_drs(1);
hot_material_tag(1) = -2;
hot_material_quality(1) = 1;
if hot_Nrs(2)~=0
    for i = 2:hot_Nrs(2)+1
        hot_Rs(i) = hot_Rs(i-1)+hot_drs(2);
        hot_material_tag(i) = hLm(1);
        hot_dr(i) = hot_drs(2);
        hot_material_quality(i) = hLq(1);
    end
else
    i = 1;
end
if hot_Nrs(3)~=0
    for k = 1:hot_Nrs(3)
        hot_Rs(k+i) = hot_Rs(k+i-1)+hot_drs(3);
        hot_material_tag(k+i) = hLm(2);
        hot_dr(k+i) = hot_drs(3);
        hot_material_quality(k+i) = hLq(2);
    end
else
    k = 0;
end
if hot_Nrs(4)~=0
    for m = 1:hot_Nrs(4)
        hot_Rs(m+k+i) = hot_Rs(m+k+i-1)+hot_drs(4);
        hot_material_tag(m+k+i) = hLm(3);
        hot_dr(m+k+i) = hot_drs(4);
        hot_material_quality(m+k+i) = hLq(3);
    end
else
    m = 0;
end
hot_Rs(m+i+k+1) = hot_Geometry(6);
hot_dr(m+i+k+1) = 0.0005;
hot_material_tag(m+i+k+1) = -1;
hot_material_quality(m+i+k+1) = 1;
if numel(hot_material_tag==-1)>1
    hot_material_tag = hot_material_tag(1:(m+i+k+1));
    hot_material_quality = hot_material_quality(1:(m+i+k+1));
end
i = 0; k = 0;m=0;

% COLD:
cold_Rs(1) = cold_Geometry(1);
cold_dr(1) = cold_drs(1);
cold_material_tag(1) = -2;
cold_material_quality(1) = 1;
if cold_Nrs(2)~=0
    for i = 2:cold_Nrs(2)+1
        cold_Rs(i) = cold_Rs(i-1)+cold_drs(2);
        cold_material_tag(i) = cLm(1);
        cold_dr(i) = cold_drs(2);
        cold_material_quality(i) = cLq(1);
    end
else
    i = 1;
end
if cold_Nrs(3)~=0
    for k = 1:cold_Nrs(3)
        cold_Rs(k+i) = cold_Rs(k+i-1)+cold_drs(3);
        cold_material_tag(k+i) = cLm(2);
        cold_dr(k+i) = cold_drs(3);
        cold_material_quality(k+i) = cLq(2);
    end
else
    k = 0;
end
if cold_Nrs(4)~=0
    for m = 1:cold_Nrs(4)
        cold_Rs(m+k+i) = cold_Rs(m+k+i-1)+cold_drs(4);
        cold_material_tag(m+k+i) = cLm(3);
        cold_dr(m+k+i) = cold_drs(4);
        cold_material_quality(m+k+i) = cLq(3);
    end
else
    m = 0;
end
cold_Rs(m+i+k+1) = cold_Geometry(6);
cold_dr(m+i+k+1) = 0.0005;
cold_material_tag(m+i+k+1) = -1;
cold_material_quality(m+i+k+1) = 1;
if numel(cold_material_tag==-1)>1
    cold_material_tag = cold_material_tag(1:(m+i+k+1));
    cold_material_quality = cold_material_quality(1:(m+i+k+1));
end


% HOT mesh:
    hot_Rs = [0,0,hot_Rs,0,0]
    hot_drs = [0,0,hot_dr,0,0,]
    hot_material_tag = [0,0,hot_material_tag,0,0]
    hot_material_quality = [0,0,hot_material_quality,0,0]
    s = zeros(1,length(hot_Rs));
    s(hot_Rs~=0) = 1;
    s(1) = 1; 
    hot_Ts = T_ambient.*s;
    hot_Rs(1) = FoC;
    hot_drs(1) = dt;
    hot_mesh_map = [hot_Ts;hot_Rs;hot_drs;hot_material_tag;hot_material_quality;hot_Ts]
    
 % COLD mesh:
    cold_Rs = [0,0,cold_Rs,0,0];
    cold_drs = [0,0,cold_dr,0,0];
    cold_material_tag = [0,0,cold_material_tag,0,0];
    cold_material_quality = [0,0,cold_material_quality,0,0];
    s = zeros(1,length(cold_Rs)); 
    s(cold_Rs~=0) = 1;
    s(1) = 1; 
    cold_Ts = T_ambient.*s;
    cold_Rs(1) = FoC;
    cold_drs(1) = dt;
    cold_mesh_map = [cold_Ts;cold_Rs;cold_drs;cold_material_tag;cold_material_quality;cold_Ts]
    
    mesh_map.cold_mesh_map = cold_mesh_map;
    mesh_map.hot_mesh_map = hot_mesh_map;
    
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

%% Cost Functions:
%%%%%%%%%%%%%

function stainless_steel_cost = stainless_steel_price_m3(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [81950,74500,55875,52150,50660,48425];
    stainless_steel_cost = interp1(N,V,N_dishes_in_field,'spline');

end

function carbon_steel_cost = carbon_steel_price_m3(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [32780,29800,22350,20860,19370,17880];
    carbon_steel_cost = interp1(N,V,N_dishes_in_field,'spline');

end

function aluminium_cost = aluminium_price_m3(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [12953,11775,8831,8243,7654,7065];
    aluminium_cost = interp1(N,V,N_dishes_in_field,'spline');

end

function insulation_cost = insulation_price_m3(material,N_dishes_in_field)
% cost as [$/m3]

        switch material
        
            case {0,-inf,inf} % null
                    V = [0,0,0,0,0,0];
            case [] % null
                    V = [0,0,0,0,0,0];
            case 1 %'MPS (MicroTherm)'
                    V = [30517.5,23010,23010,31463.25,31463.25,31463.25];
            case 2 %'Blanket (MicroTherm)'
                    V = [30517.5,23010,23010,31463.25,31463.25,31463.25];
            case 3 % 'PyroJell (XT-E)'
                    V = [31322,31322,31322,31322,31322,31322];
            case 4 % 'Majus (MicroTherm)'
                    V = [205632,205632,205632,205632,205632,205632];
            case 5 %'Rock Wool'
                    V = [10070.775,7705.9125,7705.9125,7072.065,7072.065,7072.065];
            case 6 %'Glass Wool'
                    V = [10070.775,7705.9125,7705.9125,7072.065,7072.065,7072.065];
            case 7
                    V = [205632,205632,205632,205632,205632,205632]./10;
        end

    N = [1,8,80,100,120,160];   
    insulation_cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = dish_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [804260.85,699357,615434,559486,503537,447589];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = steam_cycle_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [1045333,1194667,2730667,3157333,3584000,5461333];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = power_unit_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [4875000,9750000,48000000,60000000,66000000,85200000];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = pumps_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [16145.83,19375,124000,155000,170500,209250];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = water_treatment_system_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [13541.67,16250,52000,65000,78000,86450];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = air_system_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [28303.18,33964,108684,135855,163026,180688];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = blowers_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [20500,24600,236160,295200,354240,392616];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = instruments_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [4683,18731,179819,224773,269728,269728];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = control_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [5304,21217,203680,254600,305520,269728];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end

function cost = other_essentials_cost(N_dishes_in_field)
% cost as [$/m3]

    N = [1,8,80,100,120,160];
    V = [3288,13150,126244,157805,189366,269728];
    cost = interp1(N,V,N_dishes_in_field,'spline');

end





























