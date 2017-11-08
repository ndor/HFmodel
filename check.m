function [] = check()
% loading tables data - titles

global field_table...
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
    
    error_count = 0; note_count = 0; n = 0; k = 0; z = 0; m=0;
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
    z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
    if model_type<1; error_count = error_count+1; k = k+1; e{k} = ['# Please Select Model Type.']; end
    if location<1; error_count = error_count+1; k = k+1; e{k} = ['# Please Select Site Location.']; end
    if plant_type<1; error_count = error_count+1; k = k+1; e{k} = ['# Please Select Plant Type.']; end
    if HTS_type<1; error_count = error_count+1; k = k+1; e{k} = ['# Please Select HTS Type.']; end
    if system_pressure<1; error_count = error_count+1; k = k+1; e{k} = ['# Please Select System Pressure.']; end
    if strcmp(file_name,'') || isempty(file_name); note_count = note_count+1; z = z+1; s{z} = ['# NOTE: Output File Name will consist of model specification only']; end
    if k>2; e{2} = ['Macro Variables:']; end
    if z>2; s{2} = ['Macro Variables:']; end

    
    num_field_table = cell2mat({field_table{:,2}})';
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
%     z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
    if size(num_field_table,1)<7 || sum(isnan(num_field_table))>0
        error_count = error_count+1;
        k = k+1; e{k} = ['# Not all Field Variables are filled.'];
    else
        if abs(num_field_table(1))>90; error_count = error_count+1; k = k+1; e{k} = ['# Azimuth Orientation angle can be -90 to 90 [deg] only.']; end
        if abs(num_field_table(2))>10; error_count = error_count+1; k = k+1; e{k} = ['# Elevation Orientation angle can be -10 to 10 [deg] only.']; end
        if num_field_table(3)>30 || num_field_table(3)<1 || round(num_field_table(3))~=num_field_table(3); error_count = error_count+1; k = k+1; e{k} = ['# There can be 1 to 30 dishes in a cluster, by whole (N) numbers.']; end
        if num_field_table(4)>30 || num_field_table(4)<1 || round(num_field_table(4))~=num_field_table(4) || (num_field_table(4)~=1 && (num_field_table(4)>1 && ceil(num_field_table(4)./2)>num_field_table(4)./2)); k = k+1; e{k} = ['# Ther can be 1 to 30 clusters in a field, by whole (N), and even numbers, besides 1.']; end
        if num_field_table(5)>45 || num_field_table(5)<0; error_count = error_count+1; k = k+1; e{k} = ['# Field Tilt Angle can be 0 to 45 [deg] only.']; end
        if num_field_table(6)<39; error_count = error_count+1; k = k+1; e{k} = ['# N-S Dish Distances can not be lower than 39 [m].']; end
            % minimum dish spacing alert:
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
        if num_field_table(7)<minimum_Lew; error_count = error_count+1; k = k+1; e{k} = ['# Given your N-S Dish Distances, the E-W Dish Distances can not be lower than ' num2str(minimum_Lew) ' [m].']; end
        if num_field_table(8)<0; error_count = error_count+1; k = k+1; e{k} = ['# HEX Distance from Field Center can not be negative [m].']; end
    end
    if k>n; e{n} = ['Field Variables:']; end
    
    
    num_dishNreceiver_table = cell2mat({dishNreceiver_table{:,2}})';
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
%     z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
    if size(num_dishNreceiver_table,1)<4 || sum(isnan(num_dishNreceiver_table))>0
        error_count = error_count+1;
        k = k+1; e{k} = ['# Not all Dish & Receiver Variables are filled.'];
    else
        if num_dishNreceiver_table(1)<0; error_count = error_count+1; k = k+1; e{k} = ['# Collector Size can not be negative [m2].'];end
        if num_dishNreceiver_table(2)<1.5 || num_dishNreceiver_table(2)>5; error_count = error_count+1; k = k+1; e{k} = ['# RMS Slope Error can be 1.5 to 5 [m.rad] only.'];end
        if num_dishNreceiver_table(3)<0 || num_dishNreceiver_table(3)>=100; error_count = error_count+1; k = k+1; e{k} = ['# Base Reflectivity can be 0 to 100 [%] only.'];end
        if num_dishNreceiver_table(4)<1 || num_dishNreceiver_table(4)>=100; error_count = error_count+1; k = k+1; e{k} = ['# Max Receiver Efficiency can be 1 to 100 [%] only.'];end
    end
    if k>n; e{n} = ['Dish & Receiver Variables:']; end
    
    num_ducts_size_table = cell2mat({ducts_size_table{:,2}})';
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
%     z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
    if size(num_ducts_size_table,1)<2 || sum(isnan(num_ducts_size_table))>0
        error_count = error_count+1;
        k = k+1; e{k} = ['# Not all HTS Ducts Size (N1) are filled.'];
    else
        if HTS_type>1
            if num_ducts_size_table(1)<0.03.*Pscalar(system_pressure) || num_ducts_size_table(1)>0.1.*Pscalar(system_pressure); error_count = error_count+1; k = k+1; e{k} = ['# Given system pressure, (N1)s Hot Duct Radius can be  ' num2str(0.03.*Pscalar(system_pressure)) ' to ' num2str(0.1.*Pscalar(system_pressure)) ' [m] only.'];end
            if num_ducts_size_table(2)<num_ducts_size_table(1)+0.01 || num_ducts_size_table(2)>1.*Pscalar(system_pressure); error_count = error_count+1; k = k+1; e{k} = ['# Given system pressure and your (N1)s Hot Duct Radius, Cold Duct Radius can be ' num2str(num_ducts_size_table(1)+0.01) ' to ' num2str(1.*Pscalar(system_pressure)) ' [m] only.'];end
        else
            if num_ducts_size_table(1)<0.03.*Pscalar(system_pressure) || num_ducts_size_table(1)>0.1.*Pscalar(system_pressure); error_count = error_count+1; k = k+1; e{k} = ['# (N1)s Hot Duct Radius can be ' num2str(0.03.*Pscalar(system_pressure)) ' to ' num2str(0.1.*Pscalar(system_pressure)) ' [m] only.'];end
            if num_ducts_size_table(2)<0.03.*Pscalar(system_pressure) || num_ducts_size_table(2)>0.1.*Pscalar(system_pressure); error_count = error_count+1; k = k+1; e{k} = ['# (N1)s Cold Duct Radius can be ' num2str(0.03.*Pscalar(system_pressure)) ' to ' num2str(0.1.*Pscalar(system_pressure)) ' [m] only.'];end
        end
    end
    if k>n; e{n} = ['HTS Ducts Size (N1):']; end


%     num_insulation_thickness_table = cell2mat(insulation_table(:,2));
    num_insulation_material_table = insulation_table(:,3);
%     num_insulation_quality_table = cell2mat(insulation_table(:,4));

    layer = ['Hot','Cold'];
        for i = 1:6
            include(i) = ~strcmp(num_insulation_material_table(i),'NoN / Select Material...');
            if i>3; l = layer(4:6); t = i-3; else l = layer(1:3); t = i; end
            if include(i)>0

                if isempty(cell2mat(insulation_table(i,2)))
                    num_insulation_thickness_table(i) =0;
                else
                    num_insulation_thickness_table(i) = cell2mat(insulation_table(i,2));
                end
                if isempty(num_insulation_thickness_table(i)) || isnan(num_insulation_thickness_table(i)) || num_insulation_thickness_table(i)<0
                    num_insulation_thickness_table(i) = 0;
                    error_count = error_count+1;
                    k = k+1;
                    e{k} = ['# ' l ' Layer ' num2str(t) ' Selected Material has no thickness (Value missing or negative).'];
                end
                
                if isempty(cell2mat(insulation_table(i,4)) )
                    num_insulation_quality_table(i,4)=0;
                else
                    num_insulation_quality_table(i) = cell2mat(insulation_table(i,4));
                end
                if isempty(num_insulation_quality_table(i)) || isnan(num_insulation_quality_table(i)) || num_insulation_quality_table(i)<0
                    num_insulation_quality_table(i) = 0; 
                    error_count = error_count+1;
                    k = k+1;
                    e{k} = ['# ' l ' Layer ' num2str(t) ' Selected Material has no Quality factor (Value missing or negative).'];
                end
            else
                num_insulation_thickness_table(i) = 0;
                num_insulation_quality_table(i) = 0;
            end
        end 



    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
%     z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
    if sum(include(1:3))<1 || sum(include(4:6))<1%size(num_insulation_thickness_table(1:3),1)<1 && size(num_insulation_thickness_table(4:6),1)<1
        if sum(include(1:3))<1
            error_count = error_count+1;
            k = k+1; e{k} = ['# Hot Duct must be insulated at least with one layer (Layer Material missing).'];
        end
        if sum(include(4:6))<1
            error_count = error_count+1;
            k = k+1; e{k} = ['# Cold Duct must be insulated at least with one layer (Layer Material missing).'];
        end
        
    elseif numel(num_insulation_thickness_table)<sum(include) || numel(num_insulation_quality_table)<sum(include)
            if numel(num_insulation_thickness_table)<sum(include)
                if i<4
                    error_count = error_count+1;
                    k = k+1;
                    e{k} = ['# Hot Insulation Layers (chosen material) Thickness are not filled.'];
                else
                    error_count = error_count+1;
                    k = k+1;
                    e{k} = ['# Cold Insulation Layers (chosen material) Thickness are not filled.'];
                end
            end
            if numel(num_insulation_quality_table)<sum(include)
                if i<4
                    error_count = error_count+1;
                    k = k+1;
                    e{k} = ['# Hot Insulation Layers (chosen material) Quality are not filled.'];
                else
                    error_count = error_count+1;
                    k = k+1;
                    e{k} = ['# Cold Insulation Layers (chosen material) Quality are not filled.'];
                end
            end
        else
        hot_thick = 0;
        cold_thick = 0;
        for i = 1:6
            if ~strcmp(num_insulation_material_table(i),'NoN / Select Material...')
                if i<4
                    hot_thick  = num_insulation_thickness_table(i)+hot_thick;
                    if num_insulation_quality_table(i)<1 || num_insulation_quality_table(i)>100; error_count = error_count+1; k = k+1; e{k} = ['# Hot Insulation - Layer ' num2str(i) ' Installation Quality can be from 1 to 100 [%] only.'];end
                else
                    cold_thick  = num_insulation_thickness_table(i)+cold_thick;
                    if num_insulation_quality_table(i)<1 || num_insulation_quality_table(i)>100; error_count = error_count+1; k = k+1; e{k} = ['# Cold Insulation - Layer ' num2str(i-3) ' Installation Quality can be from 1 to 100 [%] only.'];end
                end
            end
        end
            if HTS_type>1 % insulation layer thickness  check
                    if hot_thick>90 || hot_thick<1; error_count = error_count+1; k = k+1; e{k} = ['# Hot insulation Layers Occupation can be 1 to 90 [%] of annular gap only.'];end
                    if cold_thick>2 || cold_thick<0.01; error_count = error_count+1; k = k+1; e{k} = ['# Cold insulation Layers Thickness can be 0.01 to 2 [m] only.'];end
            else
                    if hot_thick>2 || hot_thick<0.01; error_count = error_count+1; k = k+1; e{k} = ['# Hot insulation Layers Thickness can be 0.01 to 2 [m] only.'];end
                    if cold_thick>2 || cold_thick<0.01; error_count = error_count+1; k = k+1; e{k} = ['# Cold insulation Layers Thickness can be 0.01 to 2 [m] only.'];end
            end   
    end
    if k>n; e{n} = ['HTS Insulation Specification:']; end


    for i = 1:4
        for j = 1:3
            if isempty(cell2mat(HEX_table(i,j+1)))
                num_HEX_table(i,j) = NaN;
            else
                 num_HEX_table(i,j) = cell2mat(HEX_table(i,j+1));
            end
        end
    end
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
%     z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
    if (sum(~isnan(num_HEX_table(:,1)))<4 && sum(~isnan(num_HEX_table(:,1)))>0) ||...
          (sum(~isnan(num_HEX_table(:,2)))<4 && sum(~isnan(num_HEX_table(:,2)))>0) ||...
            (sum(~isnan(num_HEX_table(:,3)))<4 && sum(~isnan(num_HEX_table(:,3)))>0)
        error_count = error_count+1;
        k = k+1; e{k} = ['# Not all HEX Variables are filled.'];
    else
        line = ['HP','MP','LP'];
        for i = 1:3
            if flag(i)>0
                error_count = error_count+1; k = k+1; e{k} = ['# ' line((i.*2-1):(i.*2)) ' Steam Line is not fully defined.']; 
                if i>1; k = k+1; e{k} = ['   If line ' line((i.*2-1):(i.*2)) ' is to be discarded, its column should not contain any input.'];end
            else
                if sum(isnan(num_HEX_table(:,1)))==4; error_count = error_count+1; k = k+1; e{k} = ['# HP Steam Line should always be defined.']; k = k+1; e{k} = ['   Even a single line HEX is defined as an HP (only) line.'];end
                if i>1 && num_HEX_table(3,i)-num_HEX_table(1,i)>=num_HEX_table(3,i-1)
                    error_count = error_count+1; k = k+1; e{k} = ['# Pinched ' line((i.*2-1):(i.*2)) ' Steam Line can not be higher than  ' line(((i-1).*2-1):((i-1).*2)) '.'];
                end
                if (num_HEX_table(1,i)<5 || num_HEX_table(1,i)>30) && ~isnan(num_HEX_table(1,i)); error_count = error_count+1; k = k+1; e{k} = ['# ' line((i.*2-1):(i.*2)) ' Steam Line Pinch temperature can be 5 to 30 [*C] only.']; end
                if (num_HEX_table(2,i)<5 || num_HEX_table(2,i)>30) && ~isnan(num_HEX_table(2,i)); error_count = error_count+1; k = k+1; e{k} = ['# ' line((i.*2-1):(i.*2)) ' Steam Line Approach temperature can be 5 to 30 [*C] only.']; end
            end
        end
    end
    if k>n; e{n} = ['HEX Variables:']; end
    
    
	num_plant_table = cell2mat(plant_table(:,2))';
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
%     z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
        if size(num_plant_table,1)<4 && sum(isnan(num_plant_table))>0
            error_count = error_count+1;
            k = k+1; e{k} = ['# Not all Plant Variables are filled.'];
        else
            if sum(isnan(num_plant_table))>0 || numel(num_plant_table)<4; error_count = error_count+1; k = k+1; e{k} = ['# Some of Plant Variables are missing. Plant must be fully defined.']; end
            if numel(num_plant_table)==4
                if num_plant_table(1)>min(num_HEX_table(4,:)); error_count = error_count+1; k = k+1; e{k} = ['# Supply Water Pressure can not be higher than the lowest pressure steam line : ' num2str(min(num_HEX_table(4,:))) ' [BarA].']; end
                if num_plant_table(2)>min(num_HEX_table(3,:)); error_count = error_count+1; k = k+1; e{k} = ['# Supply Water Temperature can not be higher than the lowest temperature steam line : ' num2str(min(num_HEX_table(3,:))) ' [*C].']; end
                if num_plant_table(3)<max(num_HEX_table(3,:)+num_HEX_table(1,:)); error_count = error_count+1; k = k+1; e{k} = ['# HEX Hot Air inlet Temperature can not be lower than the highest pinched temperature steam line : ' num2str(min(num_HEX_table(3,:))) ' [*C].']; end
                if num_plant_table(4)<1 || num_plant_table(4)>100; error_count = error_count+1; k = k+1; e{k} = ['# Turbine (or Power Block) Efficiency can be 1 to 100 [%] only.']; end
            end
        end
        if k>n; e{n} = ['Plant Variables:']; end
    
              
        
        
if model_type>5
    
    
    for i = 1:8
        if isempty(cell2mat(constraints_table(i,2)))
            a(i) = NaN;
        else
            a(i) = cell2mat(constraints_table(i,2));
        end
        if isempty(cell2mat(constraints_table(i,3)))
            b(i) = NaN;
        else
            b(i) = cell2mat(constraints_table(i,3));
        end
        if isempty(cell2mat(constraints_table(i,4)))
            c(i) = NaN;
        else
            c(i) = cell2mat(constraints_table(i,4));
        end
    end
    num_constraints_table = [a',b',c'];

    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
    z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;

        if numel(a(a>0))>0 && (sum(isnan(b(a>0))) || sum(isnan(c(a>0))))
            error_count = error_count+1;
            k = k+1; e{k} = ['# Not all Optimization Constraints are filled.'];
        else
            for i = 1:8
                if a(i)>0 && sum(isnan(num_constraints_table(i,2:3))>0); error_count = error_count+1; k = k+1; e{k} = ['# Marked (checkboxed) Constraint in line ' num2str(i) ', is under defined.';'All marked Constraints must be fully defined - i.e. - Max & Min values.']; end
                if a(i)<1 && sum(~isnan(b(i))>0) || sum(~isnan(c(i)))>0;  note_count = note_count+1; z = z+1; s{z} = ['# NOTE: Filled values of Constraint in line ' num2str(i) ' will not be issued for optimization - Constraint is discarded (uncheckboxed).']; end
                if b(i)>c(i) && a(i)>0; error_count = error_count+1;error_count = error_count+1;  e = ['# Constraint Value of Min is higher than Max in line ' num2str(i) '.']; end
                if i==6
                    if (b(i)<1 || c(i)>900) && a(i)>0; error_count = error_count+1;error_count = error_count+1;  e = ['# The number of Dishes in the Field can be 1 to 900 only.']; end
                else
                    if b(i)<0 && a(i)>0; error_count = error_count+1;error_count = error_count+1;  e = ['# Marked (checkboxed) Min Constraint in line ' num2str(i) ', can not be negative.']; end
                end
            end
        end
        if k>n; e{n} = ['Optimization Constraints:']; end
        if z>m; s{m} = ['Optimization Constraints:']; end
    
        
    for i = 1:7
        if isempty(cell2mat(field_Xb_table(i,2)))
            a(i) = NaN;
        else
            a(i) = cell2mat(field_Xb_table(i,2));
        end
        if isempty(cell2mat(field_Xb_table(i,3)))
            b(i) = NaN;
        else
            b(i) = cell2mat(field_Xb_table(i,3));
        end
        if isempty(cell2mat(field_Xb_table(i,4)))
            c(i) = NaN;
        else
            c(i) = cell2mat(field_Xb_table(i,4));
        end
    end
    num_field_Xb_table = [a',b',c'];
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
    z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
        if sum(a)>0 && (sum(isnan(b(a>0)))>0 || sum(isnan(c(a>0)))>0 || numel(b)<sum(a) || numel(c)<sum(a))
            error_count = error_count+1;
            k = k+1; e{k} = ['# Not all Field Bounds are filled.'];
        else
            for i = 1:7
                if a(i)>0 && (sum(isnan(b(i))>0) || sum(isnan(c(i))>0)); error_count = error_count+1;error_count = error_count+1;  e = ['# Marked (checkboxed) Field Bounds Specification in line ' num2str(i) ', is under defined.']; k = k+1; e{k} = ['All marked Field Bounds Specification must be fully defined - i.e. - Max & Min values.']; end
                if a(i)<1 && (sum(~isnan(b(i))>0) || sum(~isnan(c(i))>0));  note_count = note_count+1; z = z+1; s{z} = ['# NOTE: Filled values of Field Bounds Specification in line ' num2str(i) ' will not be issued for optimization - Field Bounds Specification is discarded (uncheckboxed).']; end
                if b(i)>c(i) && a(i)>0; error_count = error_count+1; k = k+1; e{k} = ['# Field Bounds Specification Value of Min is higher than Max in line ' num2str(i) '.']; end
            end
            if (b(1)<-90 || c(1)>90) && a(1)>0; error_count = error_count+1; k = k+1; e{k} = ['# Azimuth Orientation Angle can be -90 to 90 [deg] only.']; end
            if (b(2)<-10 || c(2)>10) && a(2)>0; error_count = error_count+1; k = k+1; e{k} = ['# Elevation Orientation Angle can be -10 to 10 [deg] only.']; end
            if (b(3)<1 || c(3)>30) && a(3)>0; error_count = error_count+1; k = k+1; e{k} = ['# There can be 1 to 30 [#] Dishes in a Cluster.']; end
            if (b(4)<1 || c(4)>30) && a(4)>0; error_count = error_count+1; k = k+1; e{k} = ['# There can be 1 to 30 [#] Clusters in a Field.']; end
            if (b(5)<0 || c(5)>45) && a(5)>0; error_count = error_count+1; k = k+1; e{k} = ['# Field Tilt Angle can be 0 to 45 [deg] only.']; end
            if (b(6)<39 || c(6)>150) && a(6)>0; error_count = error_count+1; k = k+1; e{k} = ['# N-S Dish distance can be 39 to 250 [m] only.']; end
            if (b(7)<67 || c(7)>250) && a(7)>0; error_count = error_count+1; k = k+1; e{k} = ['# E-W Dish distance can be 69 to 250 [m] only.']; end
        end
        if k>n; e{n} = ['Field Bounds Specification:']; end
        if z>m; s{z} = ['Field Bounds Specification:']; end
        
    for i = 1:8
        if isempty(cell2mat(HTS_Xb_table(i,2)))
            a(i) = NaN;
        else
            a(i) = cell2mat(HTS_Xb_table(i,2));
        end
        if isempty(cell2mat(HTS_Xb_table(i,3)))
            b(i) = NaN;
        else
            b(i) = cell2mat(HTS_Xb_table(i,3));
        end
        if isempty(cell2mat(HTS_Xb_table(i,4)))
            c(i) = NaN;
        else
            c(i) = cell2mat(HTS_Xb_table(i,4));
        end
    end
    num_HTS_Xb_table = [a',b',c'];
    k = k+1; e{k} = ['']; k = k+1; e{k} = ['']; n = k;
    z = z+1; s{z} = ['']; z = z+1; s{z} = ['']; m = z;
        if size(num_HTS_Xb_table,1)<8
            error_count = error_count+1;
            k = k+1; e{k} = ['# Not all Field Bounds are filled.'];
        else
            hot_thick_low = 0; hot_thick_high = 0;
            cold_thick_low = 0; cold_thick_high = 0;
            for i = 1:8
                if a(i)>0 && (sum(isnan(b(i))>0) || sum(isnan(c(i))>0)); error_count = error_count+1; k = k+1; e{k} = ['# Marked (checkboxed) HTS Bounds Specification in line ' num2str(i) ', is under defined.'];k = k+1; e{k} = ['All marked HTS Bounds Specification must be fully defined - i.e. - Max & Min values.']; end
                if a(i)<1 && (sum(~isnan(b(i))>0) || sum(~isnan(c(i))>0)); note_count = note_count+1; z = z+1; s{z} = ['# NOTE: Filled values of HTS Bounds Specification in line ' num2str(i) ' will not be issued for optimization - HTS Bounds Specification is discarded (uncheckboxed).']; end
                if b(i)>c(i) && a(i)>0; error_count = error_count+1; k = k+1; e{k} = ['# HTS Bounds Specification Value of Min is higher than Max in line ' num2str(i) '.']; end

                if i >2
                    if ~strcmp(num_insulation_material_table(i-2),'NoN / Select Material...') && a(i)>0
                        if i<6
                            if a(i)>0
                                hot_thick_low  = b(i)+hot_thick_low;
                                hot_thick_high  = c(i)+hot_thick_high;
                            end
                        else
                            if a(i)>0
                                cold_thick_low  = b(i)+cold_thick_low;
                                cold_thick_high  = c(i)+cold_thick_high;
                            end
                        end  
                    end
                end
            end  
            if HTS_type>1 % insulation layer thickness  check
                    if numel(a(a(3:5)>0))>0 && (hot_thick_high>90 || hot_thick_low<1); error_count = error_count+1; k = k+1; e{k} = ['# Hot insulation Layers Occupation can be 1 to 90 [%], Min Max accordingly, of annular gap only.']; end
                    if numel(a(a(6:8)>0))>0 && (cold_thick_high>1.*Pscalar(system_pressure) || cold_thick_low<0.03.*Pscalar(system_pressure)); error_count = error_count+1; k = k+1; e{k} = ['# Cold insulation Layers Thickness can be ' num2str(0.03.*Pscalar(system_pressure)) ' to ' num2str(1.*Pscalar(system_pressure)) ' [m], Min Max accordingly, only.']; end
            else
                    if numel(a(a(3:5)>0))>0 && (hot_thick_high>1.*Pscalar(system_pressure) || hot_thick_low<0.03.*Pscalar(system_pressure)); error_count = error_count+1; k = k+1; e{k} = ['# Hot insulation Layers Thickness can be ' num2str(0.03.*Pscalar(system_pressure)) ' to ' num2str(1.*Pscalar(system_pressure)) ' [m], Min Max accordingly, only.']; end
                    if numel(a(a(6:8)>0))>0 && (cold_thick_high>1.*Pscalar(system_pressure) || cold_thick_low<0.03.*Pscalar(system_pressure)); error_count = error_count+1; k = k+1; e{k} = ['# Cold insulation Layers Thickness can be ' num2str(0.03.*Pscalar(system_pressure)) ' to ' num2str(1.*Pscalar(system_pressure)) ' [m], Min Max accordingly, only.']; end
            end    
        end    
        if k>n; e{n} = ['HTS Bounds Specification:']; end
        if z>m; s{z} = ['HTS Bounds Specification:']; end
end

            e{1} = ['User Input - ' num2str(error_count) ' Errors Found:'];
            e{2} = ['---------------------------------------------'];
            s{1} = ['User Input - ' num2str(note_count) ' Notes:'];
            s{2} = ['---------------------------------------------'];

            if error_count>0; msgbox(e,'User Input Error','error'); end
            if note_count>0; msgbox(s,'Input Note','warn'); end

error_count
            

function A = Pscalar(P)
% P is system fluid pressure [BarA]
A = (1+(P-4).*(-(1-0.83333)./6)).^1.135;


