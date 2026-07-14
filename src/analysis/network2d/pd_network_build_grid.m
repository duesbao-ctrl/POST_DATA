function [ncountGrid, validMask, xCenters, yCenters, presentCount] = pd_network_build_grid(x, y, ncount)
%BUILDGRID Convert scattered chunk rows to a validated rectangular grid.
% Missing coordinate pairs remain invalid cells; duplicate pairs are errors.

    xCenters = unique(x(:)).';
    yCenters = unique(y(:)).';
    nx = numel(xCenters);
    ny = numel(yCenters);

    [~, ix] = ismember(x(:), xCenters(:));
    [~, iy] = ismember(y(:), yCenters(:));
    lin = sub2ind([ny, nx], iy, ix);

    if numel(unique(lin)) ~= numel(lin)
        error('analyze_chunk_network2d:DuplicateGridCell', ...
            'Duplicate coordinate pairs detected in the 2D chunk data.');
    end

    ncountGrid = nan(ny, nx);
    validMask = false(ny, nx);
    ncountGrid(lin) = ncount(:);
    validMask(lin) = true;
    presentCount = numel(lin);
end
