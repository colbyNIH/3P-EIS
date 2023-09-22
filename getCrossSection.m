function cross_section = getCrossSection(chamber_name)
% this simple script returns the cross sectional area of all known
% measurment chambers

switch chamber_name
    case 'Ussing'
        cross_section = 11408980.317/10000^2; % [cm^2] from ZEISS image of D3C
    case 'Perfusion'
        cross_section = 1.12; % [cm^2]
    case 'EndOhm12'
        cross_section = 1.12; % [cm^2]
    case 'model'
        cross_section = 1; % [cm^2]
    otherwise
        warning('the specified chamber is missing a cross-sectional area!\n')
        cross_section = 1;
end

end