function distribution = pd_cumulative_areal_density(step, coordinate, densityVarsInput, direction, callerId, negativeMode, densityFactor)
%PD_CUMULATIVE_AREAL_DENSITY Shared cumulative mass distribution engine.
% Coordinates are sorted ascending. high-to-low sums values at coordinates
% greater than or equal to the current coordinate; low-to-high does the reverse.

    if nargin < 4 || isempty(direction), direction = 'high-to-low'; end
    if nargin < 5 || isempty(callerId), callerId = 'postdata:cumulative'; end
    if nargin < 6 || isempty(negativeMode), negativeMode = 'clip'; end
    if nargin < 7 || isempty(densityFactor), densityFactor = 1; end
    if ~(isnumeric(densityFactor) && isscalar(densityFactor) && ...
            isfinite(densityFactor) && densityFactor > 0)
        error([callerId, ':BadDensityFactor'], ...
            'DensityFactor must be a positive finite scalar.');
    end
    direction = lower(strtrim(pd_to_char(direction)));
    if ~any(strcmp(direction, {'high-to-low','low-to-high'}))
        error([callerId, ':BadDirection'], ...
            'CumulativeDirection must be high-to-low or low-to-high.');
    end
    negativeMode = lower(strtrim(pd_to_char(negativeMode)));
    if ~any(strcmp(negativeMode, {'clip','error'}))
        error([callerId, ':BadNegativeDensityMode'], ...
            'NegativeDensityMode must be clip or error.');
    end
    coordinate = coordinate(:);
    if numel(coordinate) ~= size(step.data, 1)
        error([callerId, ':CoordinateSizeMismatch'], ...
            'Coordinate length must equal the selected timestep row count.');
    end
    finiteCoordinate = isfinite(coordinate);
    if ~any(finiteCoordinate)
        error([callerId, ':NoFiniteCoordinate'], ...
            'No finite coordinate values are available.');
    end
    if any(~finiteCoordinate)
        warning([callerId, ':NonFiniteCoordinate'], ...
            'Ignoring %d row(s) with non-finite coordinates.', nnz(~finiteCoordinate));
    end

    densityVarsRequested = resolveDensityVars(step, densityVarsInput, callerId);
    col = step.colIndex;
    density = zeros(size(step.data, 1), numel(densityVarsRequested));
    keep = false(1, numel(densityVarsRequested));
    for k = 1:numel(densityVarsRequested)
        name = densityVarsRequested{k};
        if ~isfield(col, name)
            warning([callerId, ':MissingDensityVar'], ...
                'Density variable "%s" is absent and will be skipped.', name);
            continue;
        end
        density(:, k) = step.data(:, col.(name));
        keep(k) = true;
    end
    densityVars = densityVarsRequested(keep);
    density = density(:, keep) .* densityFactor;
    if isempty(densityVars)
        error([callerId, ':NoAvailableDensityVar'], ...
            'No requested ArealDensity column is available.');
    end

    coordinate = coordinate(finiteCoordinate);
    density = density(finiteCoordinate, :);
    nonFiniteDensity = ~isfinite(density);
    if any(nonFiniteDensity(:))
        warning([callerId, ':NonFiniteDensity'], ...
            'Treating %d non-finite density value(s) as zero.', nnz(nonFiniteDensity));
        density(nonFiniteDensity) = 0;
    end

    negativeValueCount = nnz(density < 0);
    if negativeValueCount > 0
        if strcmp(negativeMode, 'error')
            error([callerId, ':NegativeDensity'], ...
                'Found %d negative areal-density value(s).', negativeValueCount);
        end
        warning([callerId, ':NegativeDensityClipped'], ...
            'Clipping %d negative areal-density value(s) to zero.', ...
            negativeValueCount);
        density(density < 0) = 0;
    end

    [idxTotal, idx1, idx2] = findMassDensityIndices(densityVars);
    if idxTotal > 0 && idx1 > 0 && idx2 > 0
        density(:, idxTotal) = density(:, idx1) + density(:, idx2);
    end

    [coordinate, order] = sort(coordinate, 'ascend');
    density = density(order, :);
    [coordinate, ~, groups] = unique(coordinate, 'sorted');
    if numel(groups) ~= numel(coordinate)
        aggregated = zeros(numel(coordinate), size(density, 2));
        for k = 1:size(density, 2)
            aggregated(:, k) = accumarray(groups, density(:, k), ...
                [numel(coordinate), 1], @sum, 0);
        end
        density = aggregated;
    end
    if strcmp(direction, 'high-to-low')
        cumulative = flipud(cumsum(flipud(density), 1));
    else
        cumulative = cumsum(density, 1);
    end
    monotonicTolerance = 64 .* eps(max(1, max(abs(cumulative(:)))));
    isMonotonicDecreasing = all(diff(cumulative, 1, 1) <= monotonicTolerance, 1);
    if strcmp(direction, 'high-to-low') && ~all(isMonotonicDecreasing)
        error([callerId, ':InternalNonMonotonicCumulative'], ...
            'High-to-low cumulative density must be monotonic decreasing.');
    end

    isNonZero = false(1, size(cumulative, 2));
    for k = 1:size(cumulative, 2)
        finiteValues = cumulative(isfinite(cumulative(:, k)), k);
        isNonZero(k) = ~isempty(finiteValues) && any(abs(finiteValues) > 0);
    end
    plotMask = isNonZero;
    if idxTotal > 0 && idx1 > 0 && idx2 > 0 && ...
            xor(isNonZero(idx1), isNonZero(idx2)) && isNonZero(idxTotal)
        plotMask(:) = false;
        plotMask(idxTotal) = true;
    end

    distribution = struct();
    distribution.coordinate = coordinate;
    distribution.density = density;
    distribution.cumulativeDensity = cumulative;
    distribution.densityVars = densityVars;
    distribution.plotMask = plotMask;
    distribution.plottedVars = densityVars(plotMask);
    distribution.direction = direction;
    distribution.negativeDensityMode = negativeMode;
    distribution.negativeValueCount = negativeValueCount;
    distribution.isMonotonicDecreasing = isMonotonicDecreasing;
    distribution.densityFactor = densityFactor;
end

function densityVars = resolveDensityVars(step, densityVarsInput, callerId)
    if isstring(densityVarsInput)
        densityVarsInput = cellstr(densityVarsInput);
    elseif ischar(densityVarsInput)
        densityVarsInput = {densityVarsInput};
    end
    if ~isempty(densityVarsInput)
        densityVars = cellfun(@matlab.lang.makeValidName, densityVarsInput, ...
            'UniformOutput', false);
        return;
    end
    densityVars = {};
    for i = 1:numel(step.validVarNames)
        name = step.validVarNames{i};
        if ~isempty(regexp(name, 'ArealDensity$', 'once'))
            densityVars{end + 1} = name; %#ok<AGROW>
        end
    end
    if isempty(densityVars)
        error([callerId, ':NoDensityCols'], 'No ArealDensity columns were found.');
    end
end

function [idxTotal, idx1, idx2] = findMassDensityIndices(names)
    idxTotal = findSuffix(names, 'massArealDensity');
    idx1 = findSuffix(names, 'mass1ArealDensity');
    idx2 = findSuffix(names, 'mass2ArealDensity');
    if isempty(idxTotal), idxTotal = 0; end
    if isempty(idx1), idx1 = 0; end
    if isempty(idx2), idx2 = 0; end
end

function index = findSuffix(names, suffix)
    index = 0;
    expression = [regexptranslate('escape', suffix), '$'];
    for i = 1:numel(names)
        if ~isempty(regexp(names{i}, expression, 'once'))
            index = i;
            return;
        end
    end
end
