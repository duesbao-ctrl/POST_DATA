function type = pd_normalize_analysis_type(value)
%PD_NORMALIZE_ANALYSIS_TYPE Normalize public aliases to stable result types.

    type = lower(strtrim(pd_to_char(value)));
    switch type
        case {'mass-v','mass_v','massv','mass-vx','mass_vx','massvx'}
            type = 'vx';
        case {'mass-x','mass_x'}
            type = 'massx';
        case {'network','network-2d','network_2d'}
            type = 'network2d';
    end
end
