function [comp, rank, status] = pd_network_compute_component_stats( ...
        labelGrid, phaseMask, validMask, xCenters, yCenters, dx, dy, ...
        boundary, wrapXComp, wrapYComp, phaseGeom)
%PD_NETWORK_COMPUTE_COMPONENT_STATS Measure labeled phase components.
%   The label topology is supplied by PD_NETWORK_LABEL_COMPONENTS. This
%   function owns component geometry, ranking, and boundary connectivity.

    if nargin < 11
        phaseGeom = [];
    end
    numComponents = max(labelGrid(:));
    comp = emptyComponentStats();
    comp.ownerGrid = zeros(size(labelGrid));

    if numComponents > 0
        comp = measureComponents(comp, labelGrid, numComponents, phaseMask, ...
            validMask, xCenters, yCenters, dx, dy, boundary, wrapXComp, ...
            wrapYComp, phaseGeom);
    end

    rank = buildComponentSizeRank(comp);
    status = struct();
    status.percolatesX = computePercolates(comp, 'x', boundary);
    status.percolatesY = computePercolates(comp, 'y', boundary);
    status.wrapsX = computeWraps(comp, 'x', boundary);
    status.wrapsY = computeWraps(comp, 'y', boundary);
end

function comp = measureComponents(comp, labelGrid, numComponents, phaseMask, ...
        validMask, xCenters, yCenters, dx, dy, boundary, wrapXComp, ...
        wrapYComp, phaseGeom)
    [ny, nx] = size(labelGrid);
    useCutCell = isstruct(phaseGeom) && isfield(phaseGeom, 'areaGrid') && phaseGeom.enabled;
    if useCutCell
        ownerGrid = buildFractionOwnerGrid(labelGrid, phaseGeom.areaGrid > 0, ...
            validMask, boundary);
    else
        ownerGrid = labelGrid;
    end
    comp.ownerGrid = ownerGrid;
    comp.label = (1:numComponents).';
    comp.cellCount = zeros(numComponents, 1);
    comp.area = zeros(numComponents, 1);
    comp.equivDiameter = zeros(numComponents, 1);
    comp.centroidX = zeros(numComponents, 1);
    comp.centroidY = zeros(numComponents, 1);
    comp.bboxWidth = zeros(numComponents, 1);
    comp.bboxHeight = zeros(numComponents, 1);
    comp.perimeterOpen = zeros(numComponents, 1);
    comp.interfacePerimeter = zeros(numComponents, 1);
    comp.shapeFactor = nan(numComponents, 1);
    comp.touchesLeft = false(numComponents, 1);
    comp.touchesRight = false(numComponents, 1);
    comp.touchesBottom = false(numComponents, 1);
    comp.touchesTop = false(numComponents, 1);
    comp.wrapX = false(numComponents, 1);
    comp.wrapY = false(numComponents, 1);
    comp.percolatesX = false(numComponents, 1);
    comp.percolatesY = false(numComponents, 1);

    for k = 1:numComponents
        mask = (labelGrid == k);
        ownerMask = (ownerGrid == k);
        if ~any(ownerMask(:))
            ownerMask = mask;
        end
        [rows, cols] = find(ownerMask);
        [binaryRows, binaryCols] = find(mask);
        comp.cellCount(k) = numel(rows);
        if useCutCell
            comp.area(k) = sum(phaseGeom.areaGrid(ownerMask));
        else
            comp.area(k) = comp.cellCount(k) * dx * dy;
        end
        comp.equivDiameter(k) = 2 * sqrt(comp.area(k) / pi);
        if useCutCell
            [comp.centroidX(k), comp.centroidY(k)] = weightedComponentCentroid( ...
                ownerMask, phaseGeom.areaGrid, phaseGeom.centroidXGrid, ...
                phaseGeom.centroidYGrid, xCenters, yCenters);
        else
            comp.centroidX(k) = mean(xCenters(cols));
            comp.centroidY(k) = mean(yCenters(rows));
        end
        comp.bboxWidth(k) = (max(cols) - min(cols) + 1) * dx;
        comp.bboxHeight(k) = (max(rows) - min(rows) + 1) * dy;

        comp.touchesLeft(k) = any(binaryCols == 1);
        comp.touchesRight(k) = any(binaryCols == nx);
        comp.touchesBottom(k) = any(binaryRows == 1);
        comp.touchesTop(k) = any(binaryRows == ny);
        comp.wrapX(k) = wrapXComp(k);
        comp.wrapY(k) = wrapYComp(k);
        comp.percolatesX(k) = comp.touchesLeft(k) && comp.touchesRight(k);
        comp.percolatesY(k) = comp.touchesBottom(k) && comp.touchesTop(k);

        [binaryPerimeterOpen, binaryInterfacePerimeter] = ...
            computePerimeterForMask(mask, phaseMask, validMask, dx, dy, boundary);
        comp.perimeterOpen(k) = binaryPerimeterOpen;
        comp.interfacePerimeter(k) = binaryInterfacePerimeter;
        if useCutCell
            cutInterface = sum(phaseGeom.interfaceLengthGrid(ownerMask));
            comp.interfacePerimeter(k) = cutInterface;
            comp.perimeterOpen(k) = max(0, binaryPerimeterOpen - ...
                binaryInterfacePerimeter) + cutInterface;
        end
        if comp.perimeterOpen(k) > 0
            comp.shapeFactor(k) = 4 * pi * comp.area(k) / (comp.perimeterOpen(k) ^ 2);
        end
    end
end

function ownerGrid = buildFractionOwnerGrid(labelGrid, fractionMask, validMask, boundary)
    [ny, nx] = size(labelGrid);
    ownerGrid = zeros(ny, nx);
    seedMask = labelGrid > 0;
    ownerGrid(seedMask) = labelGrid(seedMask);
    workMask = fractionMask & validMask;
    if ~any(seedMask(:)) || ~any(workMask(:))
        return;
    end

    maxQueue = max(1, nnz(workMask) + nnz(seedMask));
    queue = zeros(maxQueue, 1);
    seedLin = find(seedMask);
    queue(1:numel(seedLin)) = seedLin;
    head = 1;
    tail = numel(seedLin);

    while head <= tail
        cur = queue(head);
        head = head + 1;
        [r, c] = ind2sub([ny, nx], cur);
        curLabel = ownerGrid(cur);
        for d = 1:4
            [nr, nc, hasNeighbor] = pd_network_neighbor(r, c, d, ny, nx, boundary);
            if ~hasNeighbor || ~workMask(nr, nc)
                continue;
            end
            nextLin = sub2ind([ny, nx], nr, nc);
            if ownerGrid(nextLin) ~= 0
                continue;
            end
            tail = tail + 1;
            if tail > numel(queue)
                queue = [queue; zeros(maxQueue, 1)]; %#ok<AGROW>
            end
            queue(tail) = nextLin;
            ownerGrid(nextLin) = curLabel;
        end
    end
end

function [cx, cy] = weightedComponentCentroid(ownerMask, areaGrid, ...
        centroidXGrid, centroidYGrid, xCenters, yCenters)
    area = areaGrid(ownerMask);
    cxVals = centroidXGrid(ownerMask);
    cyVals = centroidYGrid(ownerMask);
    valid = isfinite(area) & area > 0 & isfinite(cxVals) & isfinite(cyVals);
    if any(valid)
        w = area(valid);
        cx = sum(w .* cxVals(valid)) / sum(w);
        cy = sum(w .* cyVals(valid)) / sum(w);
        return;
    end

    [rows, cols] = find(ownerMask);
    if isempty(rows)
        cx = NaN;
        cy = NaN;
    else
        cx = mean(xCenters(cols));
        cy = mean(yCenters(rows));
    end
end

function [perimeterOpen, interfacePerimeter] = computePerimeterForMask( ...
        componentMask, phaseMask, validMask, dx, dy, boundary)
    [ny, nx] = size(componentMask);
    perimeterOpen = 0;
    interfacePerimeter = 0;
    [rows, cols] = find(componentMask);

    for i = 1:numel(rows)
        r = rows(i);
        c = cols(i);
        for d = 1:4
            [nr, nc, hasNeighbor] = pd_network_neighbor(r, c, d, ny, nx, boundary);
            sideLength = sideLengthByDir(d, dx, dy);
            if ~hasNeighbor || ~validMask(nr, nc)
                perimeterOpen = perimeterOpen + sideLength;
                continue;
            end
            if phaseMask(nr, nc)
                continue;
            end
            perimeterOpen = perimeterOpen + sideLength;
            interfacePerimeter = interfacePerimeter + sideLength;
        end
    end
end

function rank = buildComponentSizeRank(comp)
    labels = comp.label(:);
    areas = comp.area(:);
    cellCounts = comp.cellCount(:);
    rank = struct('label', zeros(0, 1), 'area', zeros(0, 1), ...
        'cellCount', zeros(0, 1), 'rank', zeros(0, 1), ...
        'rankByLabel', zeros(0, 1));
    if isempty(labels)
        return;
    end
    [~, order] = sortrows([-areas, labels]);
    rank.label = labels(order);
    rank.area = areas(order);
    rank.cellCount = cellCounts(order);
    rank.rank = (1:numel(order)).';
    rank.rankByLabel = zeros(max(labels), 1);
    rank.rankByLabel(rank.label) = rank.rank;
end

function value = computePercolates(comp, axisName, boundary)
    if strcmp(axisName, 'x')
        if pd_network_is_periodic_axis(boundary, 'x')
            value = NaN;
        else
            value = any(comp.percolatesX);
        end
    else
        if pd_network_is_periodic_axis(boundary, 'y')
            value = NaN;
        else
            value = any(comp.percolatesY);
        end
    end
end

function value = computeWraps(comp, axisName, boundary)
    if strcmp(axisName, 'x')
        value = pd_network_is_periodic_axis(boundary, 'x') && any(comp.wrapX);
    else
        value = pd_network_is_periodic_axis(boundary, 'y') && any(comp.wrapY);
    end
end

function comp = emptyComponentStats()
    comp = struct( ...
        'label', zeros(0, 1), ...
        'cellCount', zeros(0, 1), ...
        'area', zeros(0, 1), ...
        'equivDiameter', zeros(0, 1), ...
        'centroidX', zeros(0, 1), ...
        'centroidY', zeros(0, 1), ...
        'bboxWidth', zeros(0, 1), ...
        'bboxHeight', zeros(0, 1), ...
        'perimeterOpen', zeros(0, 1), ...
        'interfacePerimeter', zeros(0, 1), ...
        'shapeFactor', zeros(0, 1), ...
        'touchesLeft', false(0, 1), ...
        'touchesRight', false(0, 1), ...
        'touchesBottom', false(0, 1), ...
        'touchesTop', false(0, 1), ...
        'wrapX', false(0, 1), ...
        'wrapY', false(0, 1), ...
        'percolatesX', false(0, 1), ...
        'percolatesY', false(0, 1));
    comp.ownerGrid = zeros(0, 0);
end

function len = sideLengthByDir(dirId, dx, dy)
    if dirId <= 2
        len = dx;
    else
        len = dy;
    end
end
