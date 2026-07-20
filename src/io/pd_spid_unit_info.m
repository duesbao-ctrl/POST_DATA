function info = pd_spid_unit_info(unitSystem)
%PD_SPID_UNIT_INFO Conversion factors for SPID native unit systems.
% Factors convert one native value to the display/analysis units named in
% the output structure.

    if isstruct(unitSystem)
        if isfield(unitSystem, 'unitSystem')
            unitSystem = unitSystem.unitSystem;
        else
            unitSystem = '';
        end
    end
    unitSystem = lower(strtrim(pd_to_char(unitSystem)));

    info = struct();
    info.unitSystem = unitSystem;
    info.isKnown = true;
    info.lengthUmPerUnit = NaN;
    info.velocityKmSPerUnit = NaN;
    info.arealDensityMgCm2PerUnit = NaN;
    info.densityGcm3PerUnit = NaN;
    switch unitSystem
        case 'si'
            info.lengthUmPerUnit = 1e6;
            info.velocityKmSPerUnit = 1e-3;
            info.arealDensityMgCm2PerUnit = 100;
            info.densityGcm3PerUnit = 1e-3;
        case 'centimeter_gram_microsecond'
            info.lengthUmPerUnit = 1e4;
            info.velocityKmSPerUnit = 10;
            info.arealDensityMgCm2PerUnit = 1000;
            info.densityGcm3PerUnit = 1;
        case 'microscale'
            info.lengthUmPerUnit = 10;
            info.velocityKmSPerUnit = 10;
            info.arealDensityMgCm2PerUnit = 1;
            info.densityGcm3PerUnit = 1;
        otherwise
            info.isKnown = false;
    end
end
