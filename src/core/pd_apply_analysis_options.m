function request = pd_apply_analysis_options(request, catalog, tableData)
%PD_APPLY_ANALYSIS_OPTIONS Parse GUI values into a versioned request.

    [values, enabled] = pd_validate_catalog_values(catalog, tableData);
    request.analysisOptions = {};
    for i = 1:numel(catalog)
        if ~enabled(i), continue; end
        value = values{i};
        switch catalog(i).Target
            case 'chunk.dimension'
                request.chunk.dimension = value;
            case 'chunk.variable'
                request.chunk.variable = value;
            case 'chunk.dV'
                request.chunk.dV = value;
            case 'chunk.coordScale'
                request.chunk.coordScale = value;
            case 'chunk.coordRangeX'
                request.chunk.coordRangeX = value;
            case 'chunk.coordRangeY'
                request.chunk.coordRangeY = value;
            case 'chunk.coordRangeZ'
                request.chunk.coordRangeZ = value;
            case 'chunk.gradientVariable'
                request.chunk.gradientVariable = value;
            case 'chunk.gradientSmoothLevel'
                request.chunk.gradientSmoothLevel = value;
            case 'chunk.strainRateVelocityComponent'
                request.chunk.strainRateVelocityComponent = value;
            case 'chunk.strainRateDensityVariable'
                request.chunk.strainRateDensityVariable = value;
            case 'analysis'
                request.analysisOptions(end + 1:end + 2) = {catalog(i).Name, value};
            otherwise
                error('postdata:UnknownOptionTarget', ...
                    'Unknown target for option %s.', catalog(i).Name);
        end
    end
end
