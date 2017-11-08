function insulation_type = insulation_type(type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% global insulation_type
insulation_type = 'NoN / Slect Material...';
    switch type
        case 1
            insulation_type = 'MPS (MicroTherm)';
        case 2
            insulation_type = 'Blanket (MicroTherm)';
        case  3
            insulation_type =  'PyroJell (XT-E)';
        case 4
            insulation_type = 'Majus (MicroTherm)';
        case 5
            insulation_type = 'Rock Wool';
        case 6
            insulation_type = 'Glass Wool';
        case 7
            insulation_type = 'Air Cavity (Low Pressure)';
    end
end

