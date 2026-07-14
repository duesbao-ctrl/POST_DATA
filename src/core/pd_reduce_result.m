function result = pd_reduce_result(result, level)
%PD_REDUCE_RESULT Apply summary/standard/full storage policy.
% Summary always retains enough data for the primary GUI visualization.

    level = lower(strtrim(pd_to_char(level)));
    switch level
        case 'full'
            % Preserve every analyzer output.
        case 'standard'
            result = removeIfPresent(result, {'plots'});
        case 'summary'
            result = reduceSummary(result);
        otherwise
            error('postdata:BadResultLevel', 'Unknown result level: %s', level);
    end
    result.storage = struct('level', level, ...
        'note', 'summary retains primary visualization data; full preserves all analyzer fields');
end

function output = reduceSummary(input)
    common = {'filePath','clusterPath','selection','stepIndex','timestep', ...
        'schemaVersion','analysisType','software','preflight','run'};
    switch input.analysisType
        case 'chunk'
            keep = [common, {'variableRequested','variableUsed','sourceType', ...
                'dimensionRequested','dimension','x','y','value','coordScale', ...
                'coordRangeX','coordRangeY'}];
        case 'cluster'
            keep = [common, {'totalRows','selectedRows','diameter','diameterRange', ...
                'sizeModel','diameterEmptyBinMode','diameterPlotStyle','stats','fit', ...
                'meanDefinition','meanByBin','hist'}];
        case 'vx'
            keep = [common, {'velocity','density','cumulativeDensity', ...
                'densityVars','plottedVars','cumulativeDirection','velocityFactor', ...
                'velocityVar','velocityLabel','velocityUnit','negativeDensityMode', ...
                'negativeValueCount','isMonotonicDecreasing'}];
        case 'massx'
            keep = [common, {'coordinate','x','coordinateVar','coordinateFactor', ...
                'coordinateRange','coordinateLabel','coordinateUnit','density', ...
                'cumulativeDensity','densityVars','plottedVars','cumulativeDirection'}];
            keep = [keep, {'negativeDensityMode','negativeValueCount', ...
                'isMonotonicDecreasing'}];
        case 'network2d'
            keep = [common, {'thresholdN','connectivity','boundary','coordScale', ...
                'ncountVar','xCenters','yCenters','dx','dy','spacingSource', ...
                'validMask','poreMask','global','summary','stats','profile'}];
        otherwise
            keep = fieldnames(input).';
    end
    output = keepFields(input, keep);
    if strcmp(input.analysisType, 'network2d') && isfield(input, 'cutCell')
        minimal = struct('geometryMode', input.cutCell.geometryMode, ...
            'enabled', input.cutCell.enabled);
        if isfield(input.cutCell, 'pore') && isfield(input.cutCell.pore, 'fraction')
            minimal.pore = struct('fraction', input.cutCell.pore.fraction);
        end
        output.cutCell = minimal;
    end
end

function output = keepFields(input, names)
    output = struct();
    for i = 1:numel(names)
        if isfield(input, names{i}), output.(names{i}) = input.(names{i}); end
    end
end

function output = removeIfPresent(output, names)
    present = names(isfield(output, names));
    if ~isempty(present), output = rmfield(output, present); end
end
