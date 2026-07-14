function [x, y, ncount] = pd_network_extract_columns(data, col, ncountVar)
%EXTRACTREQUIREDCOLUMNS Extract and validate 2D coordinates and Ncount.
% MATLAB R2016b compatible.

    x = [];
    y = [];
    if isfield(col, 'Coord1')
        x = data(:, col.Coord1);
    elseif isfield(col, 'c_x')
        x = data(:, col.c_x);
    end

    if isfield(col, 'Coord2')
        y = data(:, col.Coord2);
    elseif isfield(col, 'c_y')
        y = data(:, col.c_y);
    end

    if isempty(x) || isempty(y)
        error('analyze_chunk_network2d:Not2DChunk', ...
            '2D chunk data requires x/y coordinates (Coord1/Coord2 or c_x/c_y).');
    end
    if ~isfield(col, ncountVar)
        error('analyze_chunk_network2d:MissingNcount', ...
            'Missing Ncount variable "%s".', ncountVar);
    end

    ncount = data(:, col.(ncountVar));
    if any(~isfinite(x)) || any(~isfinite(y)) || any(~isfinite(ncount))
        error('analyze_chunk_network2d:NonFiniteInput', ...
            'x, y, and Ncount must be finite for every present chunk cell.');
    end
end
