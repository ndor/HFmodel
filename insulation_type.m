function insulation_type = insulation_type(string)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% global insulation_type

insulation_type = 0;

    switch string
        case 'MPS (MicroTherm)'
            insulation_type = 1;
        case 'Blanket (MicroTherm)'
            insulation_type = 2;
        case 'PyroJell (XT-E)'
            insulation_type = 3;
        case 'Majus (MicroTherm)'
            insulation_type = 4;
        case 'Rock Wool'
            insulation_type = 5;
        case 'Glass Wool'
            insulation_type = 6;
        case 'Air Cavity (Low Pressure)'
            insulation_type = 7;
    end

end

