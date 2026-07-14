function morphology = pd_network_build_morphology_stats( ...
        poreMask, matrixMask, poreArea, matrixArea, dx, dy, boundary, ...
        enableSkeletonGraph, toolboxOverride)
%PD_NETWORK_BUILD_MORPHOLOGY_STATS Build skeleton and thickness metrics.
%   TOOLBOXOVERRIDE is optional and supports deterministic degradation tests.

    if nargin < 9
        toolboxOverride = [];
    end
    if isempty(toolboxOverride)
        toolbox = detectImageToolboxAvailability();
    else
        toolbox = toolboxOverride;
    end

    morphology = struct();
    morphology.toolbox = toolbox;
    morphology.network = struct( ...
        'pore', buildPhaseSkeletonStats(poreMask, poreArea, dx, dy, ...
            boundary, enableSkeletonGraph, toolbox), ...
        'matrix', buildPhaseSkeletonStats(matrixMask, matrixArea, dx, dy, ...
            boundary, enableSkeletonGraph, toolbox));
    morphology.thickness = struct( ...
        'matrix', buildThicknessStats(matrixMask, dx, dy, boundary, toolbox));
end
function net = buildPhaseSkeletonStats(phaseMask, phaseArea, dx, dy, boundary, enabled, toolbox)
    net = struct();
    net.enabled = enabled;
    net.available = toolbox.skeletonAvailable;
    net.skeletonLength = NaN;
    net.branchPoints = NaN;
    net.endPoints = NaN;
    net.branchDensity = NaN;
    net.note = '';

    if ~enabled
        net.note = 'Skeleton graph is disabled (EnableSkeletonGraph=false).';
        return;
    end
    if ~toolbox.skeletonAvailable
        net.note = 'Skeleton analysis requires Image Processing Toolbox (bwskel or bwmorph).';
        return;
    end
    if ~any(phaseMask(:))
        net.skeletonLength = 0;
        net.branchPoints = 0;
        net.endPoints = 0;
        net.branchDensity = 0;
        net.note = 'Phase is empty.';
        return;
    end

    [workMask, centerMask] = buildPeriodicMorphologyMask(phaseMask, boundary);
    skel = skeletonizeBinary(workMask, toolbox);
    branchMask = bwmorph(skel, 'branchpoints');
    endMask = bwmorph(skel, 'endpoints');

    net.skeletonLength = computeSkeletonLength(skel, centerMask, dx, dy);
    net.branchPoints = nnz(branchMask & centerMask);
    net.endPoints = nnz(endMask & centerMask);
    net.branchDensity = safeDivide(net.branchPoints, phaseArea);
    if ~all(centerMask(:))
        net.note = 'Periodic morphology metrics were computed on a tiled domain and cropped back to the center tile.';
    end
end

function thick = buildThicknessStats(matrixMask, dx, dy, boundary, toolbox)
    thick = struct('min', NaN, 'mean', NaN, 'p1', NaN, 'p5', NaN, ...
        'sampleCount', 0, 'note', '');
    if ~any(matrixMask(:))
        thick.note = 'Matrix phase is empty.';
        return;
    end
    if ~(toolbox.distanceAvailable && toolbox.skeletonAvailable)
        thick.note = 'Thickness analysis requires Image Processing Toolbox (bwdist plus bwskel/bwmorph).';
        return;
    end

    [workMask, centerMask] = buildPeriodicMorphologyMask(matrixMask, boundary);
    dist = bwdist(~workMask);
    scaleLen = effectivePixelLength(dx, dy);
    skel = skeletonizeBinary(workMask, toolbox);
    sampleMask = skel & centerMask;
    values = 2 .* dist(sampleMask) .* scaleLen;
    if isempty(values)
        values = 2 .* dist(workMask & centerMask) .* scaleLen;
        thick.note = 'Skeleton sampling was empty; thickness fell back to all matrix pixels in the center tile.';
    elseif anyPeriodicBoundary(boundary)
        thick.note = 'Thickness was sampled on the center tile of a periodic tiling.';
    end
    values = values(isfinite(values) & (values > 0));
    if isempty(values)
        return;
    end

    q = pd_stats_quantiles(values, [0.01, 0.05]);
    thick.min = min(values);
    thick.mean = mean(values);
    thick.p1 = q(1);
    thick.p5 = q(2);
    thick.sampleCount = numel(values);
end

function toolbox = detectImageToolboxAvailability()
    toolbox = struct();
    toolbox.bwskel = (exist('bwskel', 'file') == 2);
    toolbox.bwmorph = (exist('bwmorph', 'file') == 2);
    toolbox.distanceAvailable = (exist('bwdist', 'file') == 2);
    toolbox.skeletonAvailable = toolbox.bwskel || toolbox.bwmorph;
    toolbox.imageProcessingAvailable = toolbox.distanceAvailable && toolbox.skeletonAvailable;
end

function [workMask, centerMask] = buildPeriodicMorphologyMask(mask, boundary)
    [ny, nx] = size(mask);
    repY = 1;
    repX = 1;
    rowStart = 1;
    colStart = 1;
    if pd_network_is_periodic_axis(boundary, 'y')
        repY = 3;
        rowStart = ny + 1;
    end
    if pd_network_is_periodic_axis(boundary, 'x')
        repX = 3;
        colStart = nx + 1;
    end
    workMask = repmat(mask, repY, repX);
    centerMask = false(size(workMask));
    centerMask(rowStart:(rowStart + ny - 1), colStart:(colStart + nx - 1)) = true;
end

function skel = skeletonizeBinary(mask, toolbox)
    if toolbox.bwskel
        skel = bwskel(mask);
    else
        skel = bwmorph(mask, 'skel', Inf);
    end
    skel = logical(skel);
end

function len = computeSkeletonLength(skel, centerMask, dx, dy)
    [ny, nx] = size(skel);
    len = 0;
    [rows, cols] = find(skel & centerMask);
    for i = 1:numel(rows)
        r = rows(i);
        c = cols(i);
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue;
                end
                nr = r + dr;
                nc = c + dc;
                if nr < 1 || nr > ny || nc < 1 || nc > nx
                    continue;
                end
                if ~skel(nr, nc)
                    continue;
                end
                if centerMask(nr, nc)
                    if sub2ind([ny, nx], nr, nc) <= sub2ind([ny, nx], r, c)
                        continue;
                    end
                end
                len = len + neighborLength(dr, dc, dx, dy);
            end
        end
    end
end

function len = neighborLength(dr, dc, dx, dy)
    if dr == 0
        len = dx;
    elseif dc == 0
        len = dy;
    else
        len = hypot(dx, dy);
    end
end

function scaleLen = effectivePixelLength(dx, dy)
    scaleLen = sqrt(dx * dy);
end

function tf = anyPeriodicBoundary(boundary)
    tf = pd_network_is_periodic_axis(boundary, 'x') || ...
        pd_network_is_periodic_axis(boundary, 'y');
end

function y = safeDivide(a, b)
    if isempty(b) || ~isfinite(b) || b == 0
        y = NaN;
    else
        y = a ./ b;
    end
end
