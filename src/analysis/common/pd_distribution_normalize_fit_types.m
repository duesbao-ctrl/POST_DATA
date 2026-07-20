function fitTypes = pd_distribution_normalize_fit_types(value, name, errorId)
%PD_DISTRIBUTION_NORMALIZE_FIT_TYPES Normalize requested fit model names.

    if nargin < 2
        name = 'FitTypes';
    end
    if nargin < 3
        errorId = 'pd_distribution_normalize_fit_types:BadFitTypes';
    end
    if isTextScalar(value)
        textValue = lower(strtrim(toChar(value)));
        if isempty(textValue) || strcmp(textValue, 'none') || strcmp(textValue, 'off')
            fitTypes = {};
            return;
        end
        if strcmp(textValue, 'all')
            fitTypes = allTypes();
            return;
        end
        raw = strsplit(textValue, {',', ';', '|', ' '});
    elseif iscell(value)
        raw = value;
    else
        error(errorId, '%s must be a string or a cell array of strings.', name);
    end

    fitTypes = {};
    for i = 1:numel(raw)
        if isempty(raw{i})
            continue;
        end
        item = lower(strtrim(toChar(raw{i})));
        if isempty(item)
            continue;
        end
        if strcmp(item, 'all')
            fitTypes = allTypes();
            return;
        elseif strcmp(item, 'none') || strcmp(item, 'off')
            fitTypes = {};
            return;
        end
        fitTypes{end + 1} = normalizeTypeName(item, name, errorId); %#ok<AGROW>
    end
    fitTypes = uniqueStable(fitTypes);
end

function typeName = normalizeTypeName(item, name, errorId)
    switch item
        case {'powerlaw', 'power-law', 'power', 'pareto', 'powerexponent', 'power-exponent'}
            typeName = 'powerlaw';
        case {'gamma', 'gam'}
            typeName = 'gamma';
        case {'lognormal', 'log-normal', 'lognorm', 'ln'}
            typeName = 'lognormal';
        otherwise
            error(errorId, ...
                'Unsupported %s entry "%s". Use powerlaw/gamma/lognormal/all/none.', ...
                name, item);
    end
end

function values = allTypes()
    values = {'powerlaw', 'gamma', 'lognormal'};
end

function out = uniqueStable(in)
    out = {};
    for i = 1:numel(in)
        if ~any(strcmp(in{i}, out))
            out{end + 1} = in{i}; %#ok<AGROW>
        end
    end
end

function value = isTextScalar(input)
    value = ischar(input);
    if ~value && exist('isstring', 'builtin') == 5
        value = isstring(input) && isscalar(input);
    end
end

function out = toChar(value)
    if ischar(value)
        out = value;
    else
        out = char(value);
    end
end
