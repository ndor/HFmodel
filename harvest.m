function [] = harvest()
% loading tables data - titles

global  INPUT...
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
    error_count...
    note_count

%     if HTS_type>1
%         insulation_table = annulus_insulation_table;
%         HTS_Xb_table = annulus_Xb_table;
%     else
%         insulation_table = pipes_insulation_table;
%         HTS_Xb_table = pipes_Xb_table;
%     end
    

        INPUT.file_name = file_name;
    
        %% Macro Variables:

        INPUT.model_selection = model_type;
        INPUT.site_location = location;
        INPUT.plant_type = plant_type;
        INPUT.system_pressure = system_pressure; % air system pressure [Bar]
        INPUT.HTS_type = HTS_type;
        

        %% Field Variables:
        
        Xfield = cell2mat(field_table(:,2));
        INPUT.field_azimuth_or = Xfield(1); % rectangular field shape azimuthial angle tilt from north [Deg]
        INPUT.field_elevation_or = Xfield(2); % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
        INPUT.N_dishes_per_cluster = Xfield(3);%ceil(X0(3)); % number of dishes at each cluster (minimum 1)
        INPUT.N_clusters_in_field = Xfield(4);%ceil(X0(4)); % number of clusters in the field (minimum 1)
        INPUT.alpha = Xfield(5); % field parralelogram angle (from E-W line, +CW) [deg]
        INPUT.lns = Xfield(6); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
        INPUT.lew = Xfield(7); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
        INPUT.HEX_distance_from_field_center = Xfield(8); % HEX distance from center field - HEX linkage [m]
        
        
        %% Dish & Receiver Variables:
        
        Xdish = cell2mat(dishNreceiver_table(:,2));
        INPUT.Dish_effective_Area = Xdish(1); % effective reflective area of each dish [m2]
        INPUT.Slope_Error = Xdish(2);%sqrt(2.*2.5.^2); % mirror slope error (production quality) [mili.radians]
        INPUT.reflectivity = Xdish(3)./100; % mirror max reflectivity [.%]
        INPUT.receiver_peak_efficiency = Xdish(4)./100; %[.%]
        
        
        %% HTS Ducts Size (N1):
        XHTS = cell2mat(ducts_size_table(:,2));
        INPUT.HTS_hot.duct_R = XHTS(1); % N1 hot pipe radius [m]
        INPUT.HTS_cold.duct_R = XHTS(2); % N1 cold pipe radius [m]
        
        
        %% HTS Insulation Specification:

        for i = 1:6
            Xthick(i) = 0; Xquality(i) = 0;
            if isempty(cell2mat(insulation_table(i,2))) %|| (isnan(cell2mat(insulation_table(i,2))) && ~isempty(insulation_table(i,2)))
                Xthick(i) = 0;
            else
                Xthick(i) = cell2mat(insulation_table(i,2));
            end
            
            Xmaterial(i) = insulation_type(insulation_table{i,3});

            if isempty(cell2mat(insulation_table(i,4))) %|| (isnan(cell2mat(insulation_table(i,4))) && isempty(insulation_table(i,4)))
                Xquality(i) = 0;
            else
                Xquality(i) = cell2mat(insulation_table(i,4));
            end
            
        end

        
        if HTS_type>1
            INPUT.HTS_hot.insulation_occupation = Xthick(1:3)./100;
            INPUT.HTS_hot.insulation_material= Xmaterial(1:3);
            INPUT.HTS_hot.insulation_quality = Xquality(1:3)./100;

            INPUT.HTS_cold.insulation_thickness = Xthick(4:6);
            INPUT.HTS_cold.insulation_material= Xmaterial(4:6);
            INPUT.HTS_cold.insulation_quality = Xquality(4:6)./100;
        else
            INPUT.HTS_hot.insulation_thickness = Xthick(1:3);
            INPUT.HTS_hot.insulation_material= Xmaterial(1:3);
            INPUT.HTS_hot.insulation_quality = Xquality(1:3)./100;

            INPUT.HTS_cold.insulation_thickness = Xthick(4:6);
            INPUT.HTS_cold.insulation_material= Xmaterial(4:6);
            INPUT.HTS_cold.insulation_quality = Xquality(4:6)./100;
        end
        
        
        %% HEX Variables:

        for i = 1:4
            for j = 1:3
                if isempty(cell2mat(HEX_table(i,j+1))) || isnan(cell2mat(HEX_table(i,j+1))) %&& ~isempty(HEX_table(i,j+1)))
                    XHEX(i,j) = 0;
                else
                     XHEX(i,j) = cell2mat(HEX_table(i,j+1));
                end
            end
        end
               
        INPUT.HP_pinch = XHEX(1,1);
        INPUT.HP_approach = XHEX(2,1);
        INPUT.HP_temperature = XHEX(3,1);
        INPUT.HP_pressure = XHEX(4,1);
        INPUT.MP_pinch = XHEX(1,2);
        INPUT.MP_approach = XHEX(2,2);
        INPUT.MP_temperature = XHEX(3,2);
        INPUT.MP_pressure = XHEX(4,2);
        INPUT.LP_pinch = XHEX(1,3);
        INPUT.LP_approach = XHEX(2,3);
        INPUT.LP_temperature = XHEX(3,3);
        INPUT.LP_pressure = XHEX(4,3);        
        
        
        %% Plant Variables:
        
        Xplant = cell2mat(plant_table(:,2));
        INPUT.WS_pressure = Xplant(1);
        INPUT.WS_temperature = Xplant(2);
        INPUT.design_inlet_temperature = Xplant(3);
        INPUT.power_block_peak = Xplant(4)./100;

if model_type>5
    

    %% Optimization Constraints:
    
    for i = 1:8
        if isempty(cell2mat(constraints_table(i,2))) %|| (isnan(cell2mat(constraints_table(i,2))) && ~isempty(constraints_table(i,2)))
            a(i) = 0;
        else
            a(i) = cell2mat(constraints_table(i,2));
        end
        if isempty(cell2mat(constraints_table(i,3))) %|| (isnan(cell2mat(constraints_table(i,3))) && ~isempty(constraints_table(i,3)))
            b(i) = 0;
        else
            b(i) = cell2mat(constraints_table(i,3));
        end
        if isempty(cell2mat(constraints_table(i,4))) %|| (isnan(cell2mat(constraints_table(i,4))) && ~isempty(constraints_table(i,4)))
            c(i) = 0;
        else
            c(i) = cell2mat(constraints_table(i,4));
        end
    end
    
    num_constraints_table = [a',b',c'];
           
            if a(1)>0
                OPT_USER.A_max = c(1);
                OPT_USER.A_min = b(1);
            else
                OPT_USER.A_max = 0;
                OPT_USER.A_min = 0;
            end
            
            if a(2)>0
                OPT_USER.B_max = c(2);
                OPT_USER.B_min = b(2);
            else
                OPT_USER.B_max = 0;
                OPT_USER.B_min = 0;
            end
            
            if a(3)>0
                OPT_USER.C_max = c(3);
                OPT_USER.C_min = b(3);
            else
                OPT_USER.C_max = 0;
                OPT_USER.C_min = 0;
            end
            
            if a(4)>0
                OPT_USER.D_max = c(4);
                OPT_USER.D_min = b(4);
            else
                OPT_USER.D_max = 0;
                OPT_USER.D_min = 0;
            end
            
            if a(5)>0
                OPT_USER.E_max = c(5);
                OPT_USER.E_min = b(5);
            else
                OPT_USER.E_max = 0;
                OPT_USER.E_min = 0;
            end
            
            if a(6)>0
                OPT_USER.F_max = c(6);
                OPT_USER.F_min = b(6);
            else
                OPT_USER.F_max = 0;
                OPT_USER.F_min = 0;
            end
            
            if a(7)>0
                OPT_USER.G_max = c(7);
                OPT_USER.G_min = b(7);
            else
                OPT_USER.G_max = 0;
                OPT_USER.G_min = 0;
            end
            
            if a(8)>0
                OPT_USER.H_max = c(8);
                OPT_USER.H_min = b(8);
            else
                OPT_USER.H_max = 0;
                OPT_USER.H_min = 0;
            end  
            
            INPUT.OPT_USER = OPT_USER;
            
            
            %% Field Bounds Specification:

    for i = 1:7
        if isempty(cell2mat(field_Xb_table(i,2))) %|| (isnan(cell2mat(field_Xb_table(i,2))) && ~isempty(field_Xb_table(i,2)))
            a(i) = 0;
        else
            a(i) = cell2mat(field_Xb_table(i,2));
        end
        if isempty(cell2mat(field_Xb_table(i,3))) %|| (isnan(cell2mat(field_Xb_table(i,3))) && ~isempty(field_Xb_table(i,3)))
            b(i) = 0;
        else
            b(i) = cell2mat(field_Xb_table(i,3));
        end
        if isempty(cell2mat(field_Xb_table(i,4))) %|| (isnan(cell2mat(field_Xb_table(i,4))) && ~isempty(field_Xb_table(i,4)))
            c(i) = 0;
        else
            c(i) = cell2mat(field_Xb_table(i,4));
        end
    end
    num_field_Xb_table = [a',b',c'];
    
            INPUT.UB.field_azimuth_or = c(1); % rectangular field shape azimuthial angle tilt from north [Deg]
            INPUT.UB.field_elevation_or = c(2); % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
            INPUT.UB.N_dishes_per_cluster = c(3); % number of dishes at each cluster (minimum 1)
            INPUT.UB.N_clusters_in_field = c(4); % number of clusters in the field (minimum 1)
            INPUT.UB.alpha = c(5); % field parralelogram angle (from E-W line, +CW) [deg]
            INPUT.UB.lns = c(6); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
            INPUT.UB.lew = c(7); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
            
            INPUT.LB.field_azimuth_or = b(1); % rectangular field shape azimuthial angle tilt from north [Deg]
            INPUT.LB.field_elevation_or = b(2); % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
            INPUT.LB.N_dishes_per_cluster = b(3); % number of dishes at each cluster (minimum 1)
            INPUT.LB.N_clusters_in_field = b(4); % number of clusters in the field (minimum 1)
            INPUT.LB.alpha = b(5); % field parralelogram angle (from E-W line, +CW) [deg]
            INPUT.LB.lns = b(6); % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
            INPUT.LB.lew = b(7); % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
            
            
            %% HTS Bounds Specification:
            
    for i = 1:8
        if isempty(cell2mat(HTS_Xb_table(i,2))) %|| (isnan(cell2mat(HTS_Xb_table(i,2))) && ~isempty((HTS_Xb_table(i,2))))
            a(i) = 0;
        else
            a(i) = cell2mat(HTS_Xb_table(i,2));
        end
        if isempty(cell2mat(HTS_Xb_table(i,3))) %|| (isnan(cell2mat(HTS_Xb_table(i,3))) && ~isempty((HTS_Xb_table(i,3))))
            b(i) = 0;
        else
            b(i) = cell2mat(HTS_Xb_table(i,3));
        end
        if isempty(cell2mat(HTS_Xb_table(i,4))) %|| (isnan(cell2mat(HTS_Xb_table(i,4))) && ~isempty((HTS_Xb_table(i,4))))
            c(i) = 0;
        else
            c(i) = cell2mat(HTS_Xb_table(i,4));
        end
    end
    num_HTS_Xb_table = [a',b',c'];
      
    
            INPUT.UB.HTS_hot.duct_R = c(1); % N1 hot pipe radius [m]
            INPUT.UB.HTS_cold.duct_R = c(2); % N1 cold pipe radius [m]
                       
            INPUT.LB.HTS_hot.duct_R = b(1); % N1 hot pipe radius [m]
            INPUT.LB.HTS_cold.duct_R = b(2); % N1 cold pipe radius [m] 
            
       if HTS_type>1
            INPUT.UB.HTS_hot.insulation_occupation = c(3:5)./100;
            INPUT.UB.HTS_cold.insulation_thickness = c(6:8);          
            INPUT.LB.HTS_hot.insulation_occupation = b(3:5)./100;
            INPUT.LB.HTS_cold.insulation_thickness = b(6:8);
       else
            INPUT.UB.HTS_hot.insulation_thickness = c(3:5);
            INPUT.UB.HTS_cold.insulation_thickness = c(6:8);          
            INPUT.LB.HTS_hot.insulation_thickness = b(3:5);
            INPUT.LB.HTS_cold.insulation_thickness = b(6:8);
        end
            
%             INPUT.UB.HTS_hot.insulation_occupation = c(3:5)./100;
%             INPUT.UB.HTS_hot.insulation_material= c(4);
%             INPUT.UB.HTS_hot.insulation_quality = c(5)./100;
% 
%             INPUT.UB.HTS_cold.insulation_thickness = c(6:8);
%             INPUT.UB.HTS_cold.insulation_material= c(7);
%             INPUT.UB.HTS_cold.insulation_quality = c(8)./100;
%             
% 
%             
%             INPUT.LB.HTS_hot.insulation_occupation = b(3:5)./100;
%             INPUT.LB.HTS_hot.insulation_material= b(4);
%             INPUT.LB.HTS_hot.insulation_quality = b(5)./100;
% 
%             INPUT.LB.HTS_cold.insulation_thickness = b(6:8);
%             INPUT.LB.HTS_cold.insulation_material= b(7);
%             INPUT.LB.HTS_cold.insulation_quality = b(8)./100;
  
end

