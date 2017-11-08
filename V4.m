function [] = V4()
%V4 Summary of this function goes here
%   Detailed explanation goes here


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


        switch location
            
            case 1%'Wuhai'
                INPUT.latitude = 39.65;
                INPUT.altitude = 1109;
                INPUT.site_name = 'Wuhai';
                INPUT.TMY_file_name = 'Wuhai_TMY.xls';
                
            case 2%'Kuqa'
                INPUT.latitude = 41.7167;
                INPUT.altitude = 1100;
                INPUT.site_name = 'Kuqa';
                INPUT.TMY_file_name = 'Kuqa_TMY.xls';      
                
            case 3%'Sde-Boker'
                INPUT.latitude = 30.868889;
                INPUT.altitude = 475;
                INPUT.site_name = 'Sde-Boker';
                INPUT.TMY_file_name = 'Sde-Boker_TMY.xls';
                            
            case 4%'Sdom'
                INPUT.latitude = 31.171109;
                INPUT.altitude = -418;
                INPUT.site_name = 'Sdom';
                INPUT.TMY_file_name = 'Sdom_ TMY.xls';

            case 5%'KJC'
                INPUT.latitude = 35.0316;
                INPUT.altitude = 754.9896;
                INPUT.site_name = 'KJC';
                INPUT.TMY_file_name = 'KJC_TMY.xls';

        end

        if HTS_type > 1 % annul;us HTS
            switch model_type

                case 1

                        % steady state model (V4.1) work points run only:
                        INPUT.TMY_file_name = 'DNI4steady.xls';
                        V41_annulus(); 

                case 2 

                        % transient state model (V4.2) work points run only:
                        INPUT.TMY_file_name = 'DNI4transient.xls';
                        V42_annulus();    

                case 3

                        % steady state model (V4.1) run only:
                        V41_annulus();

                case 4 

                        % transient state model (V4.2) run only:
                        V42_annulus(); 

                case 5

                        % transient & steady state models (V4.1 & V4.2) run only:
                        INPUT.file_name = [INPUT.file_name 'Steady State Heat Transfer'];
                        V41_annulus();
                        INPUT.file_name = [INPUT.file_name 'Transient State Heat Transfer'];
                        V42_annulus(); 

                case 6 

                        % run optimization model and than workpoints by steady state model(V4.3 & V4.1) only:
                        V43_annulus(); % steady state optimization model - V4.3
                        INPUT.TMY_file_name = 'DNI4steady.xls';
                        V41_annulus(); 

                case 7

                        % run optimization (V4.3) and than full run of model steady state (V4.1):
                        V43_annulus();
                        V41_annulus();

                case 8

                        % run optimization (V4.3) and than full run of model steady state (V4.1):
                        V43_annulus();
                        V42_annulus();

                case 9

                        % run optimization (V4.3) and than full run of both models - steady state (V4.1) & transient state (V4.2):
                        V43_annulus();
                        INPUT.file_name = [INPUT.file_name 'Steady State Heat Transfer'];
                        V41_annulus();
                        INPUT.file_name = [INPUT.file_name 'Transient State Heat Transfer'];
                        V42_annulus();

            end
            
        else % HTS pipes

            switch model_type

                case 1

                        % steady state model (V4.1) work points run only:
                        INPUT.TMY_file_name = 'DNI4steady.xls';
                        V41_pipes(); 

                case 2 

                        % transient state model (V4.2) work points run only:
                        INPUT.TMY_file_name = 'DNI4transient.xls';
                        V42_pipes();    

                case 3

                        % steady state model (V4.1) run only:
                        V41_pipes();

                case 4 

                        % transient state model (V4.2) run only:
                        V42_pipes(); 

                case 5

                        % transient & steady state models (V4.1 & V4.2) run only:
                        INPUT.file_name = [INPUT.file_name 'Steady State Heat Transfer'];
                        V41_pipes();
                        INPUT.file_name = [INPUT.file_name 'Transient State Heat Transfer'];
                        V42_pipes(); 

                case 6 

                        % run optimization model and than workpoints by steady state model(V4.3 & V4.1) only:
                        V43_pipes(); % steady state optimization model - V4.3
                        INPUT.TMY_file_name = 'DNI4steady.xls';
                        V41_pipes(); 

                case 7

                        % run optimization (V4.3) and than full run of model steady state (V4.1):
                        V43_pipes();
                        V41_pipes();

                case 8

                        % run optimization (V4.3) and than full run of model steady state (V4.1):
                        V43_pipes();
                        V42_pipes();

                case 9

                        % run optimization (V4.3) and than full run of both models - steady state (V4.1) & transient state (V4.2):
                        V43_pipes();
                        INPUT.file_name = [INPUT.file_name 'Steady State Heat Transfer'];
                        V41_pipes();
                        INPUT.file_name = [INPUT.file_name 'Transient State Heat Transfer'];
                        V42_pipes();

            end
            
        end




