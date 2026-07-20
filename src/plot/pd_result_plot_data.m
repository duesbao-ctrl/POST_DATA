function output = pd_result_plot_data(result, viewName)
%PD_RESULT_PLOT_DATA Return the exact numeric data used by one result view.
% The contract is shared by rendering, GUI inspection, clipboard copy, and
% CSV export. All table values are numeric and MATLAB R2016b compatible.

    if ~isstruct(result) || ~isfield(result, 'analysisType')
        error('postdata:plotData:BadResult', ...
            'Result must contain analysisType metadata from postdata_run.');
    end
    views = pd_result_plot_views(result);
    if nargin < 2 || isempty(viewName)
        if isempty(views)
            error('postdata:plotData:NoViews', ...
                'The result does not provide any plot views.');
        end
        viewName = views(1).Id;
    end
    viewName = lower(strtrim(pd_to_char(viewName)));
    viewIndex = find(strcmp(viewName, {views.Id}), 1, 'first');
    if isempty(viewIndex)
        error('postdata:plotData:BadView', ...
            'View "%s" is unavailable for this result.', viewName);
    end

    output = emptyContract();
    output.ViewId = viewName;
    output.Label = views(viewIndex).Label;
    output.AnalysisType = lower(strtrim(pd_to_char(result.analysisType)));
    switch output.AnalysisType
        case 'chunk'
            output = chunkData(output, result, viewName);
        case 'cluster'
            output = clusterData(output, result, viewName);
        case 'vx'
            output = massData(output, result, viewName, 'velocity');
        case 'massx'
            output = massXData(output, result, viewName);
        case 'network2d'
            output = networkData(output, result, viewName);
        otherwise
            error('postdata:plotData:UnknownType', ...
                'Unsupported analysis type: %s', output.AnalysisType);
    end
    output.RowCount = size(output.Values, 1);
    output.ColumnCount = size(output.Values, 2);
end

function output = chunkData(output, result, viewName)
    variableName = result.variableUsed;
    if strcmp(viewName, 'histogram')
        values = result.value(isfinite(result.value));
        if isempty(values)
            centers = zeros(0, 1);
            counts = zeros(0, 1);
        else
            nBins = max(5, min(50, round(sqrt(numel(values)))));
            [counts, centers] = hist(values, nBins); %#ok<HIST>
            centers = centers(:);
            counts = counts(:);
        end
        output.Kind = 'bar';
        output.ColumnNames = {variableName, 'count'};
        output.Values = [centers, counts];
        return;
    end
    if any(strcmp(viewName, {'profile-x','profile-y'}))
        if strcmp(viewName, 'profile-x')
            coordinate = result.x;
            coordinateName = 'x';
        else
            coordinate = result.y;
            coordinateName = 'y';
        end
        [centers, means] = coordinateMean(coordinate, result.value);
        output.Kind = 'line';
        output.ColumnNames = {coordinateName, ['mean_', variableName]};
        output.SeriesNames = {['mean_', variableName]};
        output.Values = [centers(:), means(:)];
        return;
    end
    if isfield(result, 'dimension') && strcmpi(result.dimension, '3d')
        output.Kind = 'field3d';
        output.ColumnNames = {'x', 'y', 'z', variableName};
        output.Values = [result.x(:), result.y(:), ...
            result.coordZ(:), result.value(:)];
        return;
    end
    if isempty(result.y)
        output.Kind = 'line';
        output.ColumnNames = {'x', variableName};
        output.SeriesNames = {variableName};
        output.Values = [result.x(:), result.value(:)];
        return;
    end
    output.Kind = 'field2d';
    output.ColumnNames = {'x', 'y', variableName};
    output.Values = [result.x(:), result.y(:), result.value(:)];
    [output.IsRectangularGrid, output.GridX, output.GridY, output.GridZ] = ...
        rectangularGrid(result.x, result.y, result.value);
end

function output = clusterData(output, result, viewName)
    if strcmp(viewName, 'cdf')
        values = sort(result.diameter(:));
        probability = (1:numel(values)).' ./ max(1, numel(values));
        output.Kind = 'line';
        output.ColumnNames = {'equivalentDiameter', 'cumulativeProbability'};
        output.SeriesNames = {'cumulativeProbability'};
        output.Values = [values, probability];
        return;
    end
    if strcmp(viewName, 'mean')
        output.Kind = 'line';
        output.ColumnNames = {'positionBin', 'meanEquivalentDiameter'};
        output.SeriesNames = {'meanEquivalentDiameter'};
        output.Values = [result.meanByBin.centers(:), ...
            result.meanByBin.meanDiameter(:)];
        return;
    end
    if strcmp(viewName, 'probability') && isfield(result, 'hist')
        centers = result.hist.centers(:);
        counts = result.hist.prob(:);
        valueName = 'probability';
    else
        values = result.diameter(:);
        if isempty(values)
            centers = zeros(0, 1);
            counts = zeros(0, 1);
        else
            nBins = max(5, min(50, round(sqrt(numel(values)))));
            [counts, centers] = hist(values, nBins); %#ok<HIST>
            centers = centers(:);
            counts = counts(:);
        end
        valueName = 'count';
    end
    output.Kind = 'bar';
    output.ColumnNames = {'equivalentDiameter', valueName};
    output.Values = [centers, counts];
end

function output = massData(output, result, viewName, coordinateField)
    if strcmp(viewName, 'differential')
        values = result.density;
        valuePrefix = 'density_';
    else
        values = result.cumulativeDensity;
        valuePrefix = 'cumulative_';
    end
    coordinate = result.(coordinateField);
    headers = {coordinateField};
    for i = 1:numel(result.densityVars)
        name = regexprep(result.densityVars{i}, '[^A-Za-z0-9_]', '_');
        headers{end + 1} = [valuePrefix, name]; %#ok<AGROW>
    end
    output.Kind = 'line';
    output.ColumnNames = headers;
    output.SeriesNames = result.densityVars;
    output.Values = [coordinate(:), values];
end

function output = massXData(output, result, viewName)
    if strcmp(viewName, 'particle-count')
        values = result.localParticleCount;
        valuePrefix = 'particle_count_';
    elseif strcmp(viewName, 'differential')
        values = result.density;
        valuePrefix = 'density_';
    else
        values = result.cumulativeDensity;
        valuePrefix = 'cumulative_';
    end
    headers = {'coordinate'};
    for i = 1:numel(result.densityVars)
        name = regexprep(result.densityVars{i}, '[^A-Za-z0-9_]', '_');
        headers{end + 1} = [valuePrefix, name]; %#ok<AGROW>
    end
    output.Kind = 'line';
    output.ColumnNames = headers;
    output.SeriesNames = result.densityVars;
    output.Values = [result.coordinate(:), values];
end

function output = networkData(output, result, viewName)
    if any(strcmp(viewName, {'pore-label','matrix-label'}))
        if strcmp(viewName, 'pore-label')
            phaseName = 'pore';
        else
            phaseName = 'matrix';
        end
        grid = result.(phaseName).labelGrid;
        output = gridData(output, result.xCenters, result.yCenters, ...
            grid, [phaseName, 'Label']);
        output.Kind = 'label-grid';
        return;
    end
    if any(strcmp(viewName, {'pore-diameter','matrix-diameter'}))
        if strcmp(viewName, 'pore-diameter')
            phaseName = 'pore';
        else
            phaseName = 'matrix';
        end
        values = result.(phaseName).components.equivDiameter;
        if isempty(values)
            centers = zeros(0, 1);
            counts = zeros(0, 1);
        else
            nBins = max(3, min(40, round(sqrt(numel(values)))));
            [counts, centers] = hist(values, nBins); %#ok<HIST>
            centers = centers(:);
            counts = counts(:);
        end
        output.Kind = 'bar';
        output.ColumnNames = {'equivalentDiameter', 'count'};
        output.Values = [centers, counts];
        return;
    end
    if any(strcmp(viewName, {'profile-x','profile-y'}))
        axisName = viewName(end);
        profile = result.profile.(axisName);
        output.Kind = 'line';
        output.ColumnNames = {axisName, 'porosity'};
        output.SeriesNames = {'porosity'};
        output.Values = [profile.centers(:), profile.porosity(:)];
        return;
    end
    phase = double(result.poreMask);
    phase(~result.validMask) = NaN;
    if isfield(result, 'cutCell') && isfield(result.cutCell, 'pore') && ...
            isfield(result.cutCell.pore, 'fraction')
        phase = result.cutCell.pore.fraction;
    end
    output = gridData(output, result.xCenters, result.yCenters, ...
        phase, 'poreFraction');
    output.Kind = 'phase-grid';
end

function output = gridData(output, xValues, yValues, grid, valueName)
    [xGrid, yGrid] = meshgrid(xValues, yValues);
    output.ColumnNames = {'x', 'y', valueName};
    output.Values = [xGrid(:), yGrid(:), grid(:)];
    output.IsRectangularGrid = true;
    output.GridX = xValues(:);
    output.GridY = yValues(:);
    output.GridZ = grid;
end

function output = emptyContract()
    output = struct();
    output.AnalysisType = '';
    output.ViewId = '';
    output.Label = '';
    output.Kind = '';
    output.ColumnNames = {};
    output.SeriesNames = {};
    output.Values = zeros(0, 0);
    output.IsRectangularGrid = false;
    output.GridX = [];
    output.GridY = [];
    output.GridZ = [];
    output.RowCount = 0;
    output.ColumnCount = 0;
end

function [centers, means] = coordinateMean(coordinate, values)
    valid = isfinite(coordinate) & isfinite(values);
    [centers, ~, groups] = unique(coordinate(valid));
    means = accumarray(groups, values(valid), [], @mean);
end

function [isGrid, xValues, yValues, valueGrid] = rectangularGrid(x, y, value)
    x = x(:);
    y = y(:);
    value = value(:);
    [xValues, ~, xIndex] = unique(x);
    [yValues, ~, yIndex] = unique(y);
    isGrid = numel(xValues) * numel(yValues) == numel(value);
    valueGrid = [];
    if ~isGrid
        return;
    end
    linearIndex = sub2ind([numel(yValues), numel(xValues)], yIndex, xIndex);
    if numel(unique(linearIndex)) ~= numel(linearIndex)
        isGrid = false;
        return;
    end
    valueGrid = nan(numel(yValues), numel(xValues));
    valueGrid(linearIndex) = value;
end
