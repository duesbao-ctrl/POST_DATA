function profiles = pd_network_build_directional_profiles(poreMask, validMask, xCenters, yCenters, dx, dy, opt, porePhase, cutCell)
    if nargin < 9
        cutCell = [];
    end
    profiles = struct();
    if strcmp(opt.ProfileAxis, 'x') || strcmp(opt.ProfileAxis, 'both')
        profiles.x = buildDirectionalProfile(poreMask, validMask, porePhase, ...
            xCenters, yCenters, dx, dy, 'x', round(opt.ProfileNumBins), opt, cutCell);
    end
    if strcmp(opt.ProfileAxis, 'y') || strcmp(opt.ProfileAxis, 'both')
        profiles.y = buildDirectionalProfile(poreMask, validMask, porePhase, ...
            xCenters, yCenters, dx, dy, 'y', round(opt.ProfileNumBins), opt, cutCell);
    end
end

function profile = buildDirectionalProfile(poreMask, validMask, porePhase, ...
        xCenters, yCenters, dx, dy, axisName, nBins, opt, cutCell)
    if nargin < 11
        cutCell = [];
    end
    cellArea = dx * dy;
    [axisCenters, coordMin, coordMax, perpAxis, axisStep] = resolveProfileAxis(axisName, xCenters, yCenters, dx, dy);
    coordRange = getProfileRange(opt, axisName);
    edges = buildAxisEdges(coordRange, coordMin, coordMax, axisStep, nBins);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    binIdCenter = discretize(axisCenters(:), edges);
    useCutCell = isstruct(cutCell) && isfield(cutCell, 'enabled') && cutCell.enabled;
    if useCutCell
        interfaceLength = computeInterfaceLengthProfileCutCell(cutCell, axisName, edges);
    else
        interfaceLength = computeInterfaceLengthProfile(poreMask, validMask, xCenters, yCenters, dx, dy, axisName, edges);
    end

    validArea = zeros(1, nBins);
    poreArea = zeros(1, nBins);
    matrixArea = zeros(1, nBins);
    porosity = nan(1, nBins);
    specificInterfaceBulk = nan(1, nBins);
    specificInterfacePore = nan(1, nBins);
    specificInterfaceMatrix = nan(1, nBins);

    for i = 1:nBins
        idxAxis = find(binIdCenter == i);
        if isempty(idxAxis)
            continue;
        end

        [poreSub, validSub] = extractSliceSubgrid(poreMask, validMask, axisName, idxAxis);
        validArea(i) = nnz(validSub) * cellArea;
        if useCutCell
            [poreAreaSub, matrixAreaSub] = extractSliceCutCellAreas(cutCell, axisName, idxAxis);
            poreArea(i) = sum(poreAreaSub(:));
            matrixArea(i) = sum(matrixAreaSub(:));
        else
            poreArea(i) = nnz(poreSub) * cellArea;
            matrixArea(i) = nnz(validSub & ~poreSub) * cellArea;
        end
        porosity(i) = safeDivide(poreArea(i), validArea(i));
        specificInterfaceBulk(i) = safeDivide(interfaceLength(i), validArea(i));
        specificInterfacePore(i) = safeDivide(interfaceLength(i), poreArea(i));
        specificInterfaceMatrix(i) = safeDivide(interfaceLength(i), matrixArea(i));

    end

    connProfile = buildComponentConnectivityProfile(porePhase.components, porePhase.connectivity, ...
        axisName, perpAxis, edges);

    profile = struct();
    profile.axis = axisName;
    profile.perpendicularAxis = perpAxis;
    profile.range = [edges(1), edges(end)];
    profile.edges = edges;
    profile.centers = centers;
    profile.validArea = validArea;
    profile.poreArea = poreArea;
    profile.matrixArea = matrixArea;
    profile.porosity = porosity;
    profile.interfaceLength = interfaceLength;
    profile.specificInterfaceBulk = specificInterfaceBulk;
    profile.specificInterfacePore = specificInterfacePore;
    profile.specificInterfaceMatrix = specificInterfaceMatrix;
    profile.connectivityFlag = connProfile.flag;
    profile.connectivityFraction = connProfile.fraction;
    profile.connectivityArea = connProfile.area;
    profile.connectivityComponentCount = connProfile.count;
    profile.componentPoreArea = connProfile.componentPoreArea;
end

function [axisCenters, coordMin, coordMax, perpAxis, axisStep] = resolveProfileAxis(axisName, xCenters, yCenters, dx, dy)
    if strcmp(axisName, 'x')
        axisCenters = xCenters;
        coordMin = xCenters(1) - 0.5 * dx;
        coordMax = xCenters(end) + 0.5 * dx;
        perpAxis = 'y';
        axisStep = dx;
    else
        axisCenters = yCenters;
        coordMin = yCenters(1) - 0.5 * dy;
        coordMax = yCenters(end) + 0.5 * dy;
        perpAxis = 'x';
        axisStep = dy;
    end
end

function interfaceLength = computeInterfaceLengthProfile(poreMask, validMask, xCenters, yCenters, dx, dy, axisName, edges)
    [ny, nx] = size(poreMask);
    nBins = numel(edges) - 1;
    interfaceLength = zeros(1, nBins);

    for r = 1:ny-1
        for c = 1:nx
            if ~(validMask(r, c) && validMask(r+1, c))
                continue;
            end
            if poreMask(r, c) == poreMask(r+1, c)
                continue;
            end
            if strcmp(axisName, 'x')
                coord = xCenters(c);
            else
                coord = 0.5 * (yCenters(r) + yCenters(r+1));
            end
            binId = pd_network_locate_bin(coord, edges);
            if binId > 0
                interfaceLength(binId) = interfaceLength(binId) + dx;
            end
        end
    end

    for r = 1:ny
        for c = 1:nx-1
            if ~(validMask(r, c) && validMask(r, c+1))
                continue;
            end
            if poreMask(r, c) == poreMask(r, c+1)
                continue;
            end
            if strcmp(axisName, 'x')
                coord = 0.5 * (xCenters(c) + xCenters(c+1));
            else
                coord = yCenters(r);
            end
            binId = pd_network_locate_bin(coord, edges);
            if binId > 0
                interfaceLength(binId) = interfaceLength(binId) + dy;
            end
        end
    end
end

function [poreSub, validSub] = extractSliceSubgrid(poreMask, validMask, axisName, idxAxis)
    if strcmp(axisName, 'x')
        poreSub = poreMask(:, idxAxis);
        validSub = validMask(:, idxAxis);
    else
        poreSub = poreMask(idxAxis, :);
        validSub = validMask(idxAxis, :);
    end
end

function range = getProfileRange(opt, axisName)
    if strcmp(axisName, 'x')
        range = opt.ProfileRangeX;
    else
        range = opt.ProfileRangeY;
    end
end

function edges = buildAxisEdges(coordRange, coordMin, coordMax, axisStep, nBins)
    if isempty(coordRange)
        minEdge = coordMin;
        maxEdge = coordMax;
    else
        minEdge = coordRange(1);
        maxEdge = coordRange(2);
    end
    if minEdge == maxEdge
        pad = 0.5 * axisStep;
        minEdge = minEdge - pad;
        maxEdge = maxEdge + pad;
    end
    edges = linspace(minEdge, maxEdge, nBins + 1);
end

function profile = buildComponentConnectivityProfile(comp, connectivity, axisName, perpAxis, edges)
    nBins = numel(edges) - 1;
    profile = struct('flag', false(1, nBins), ...
        'fraction', zeros(1, nBins), ...
        'area', zeros(1, nBins), ...
        'count', zeros(1, nBins), ...
        'componentPoreArea', zeros(1, nBins));
    if isempty(comp.label)
        return;
    end

    if strcmp(axisName, 'x')
        centroid = comp.centroidX;
    else
        centroid = comp.centroidY;
    end
    binId = nan(size(centroid));
    for k = 1:numel(centroid)
        binId(k) = pd_network_locate_bin(centroid(k), edges);
    end
    connectedIds = getConnectivityIds(connectivity, perpAxis);
    connectedMask = ismember(comp.label, connectedIds);

    for i = 1:nBins
        idxBin = (binId == i);
        if ~any(idxBin)
            continue;
        end
        profile.componentPoreArea(i) = sum(comp.area(idxBin));
        idxConnected = idxBin & connectedMask;
        profile.flag(i) = any(idxConnected);
        profile.count(i) = sum(idxConnected);
        profile.area(i) = sum(comp.area(idxConnected));
        profile.fraction(i) = safeDivide(profile.area(i), profile.componentPoreArea(i));
        if isnan(profile.fraction(i))
            profile.fraction(i) = 0;
        end
    end
end

function ids = getConnectivityIds(connectivity, axisName)
    if strcmp(axisName, 'x')
        ids = connectivity.x.componentIds;
    else
        ids = connectivity.y.componentIds;
    end
end

function interfaceLength = computeInterfaceLengthProfileCutCell(cutCell, axisName, edges)
    nBins = numel(edges) - 1;
    interfaceLength = zeros(1, nBins);
    if isempty(cutCell.interfaceSegmentLength)
        return;
    end
    if strcmp(axisName, 'x')
        coord = cutCell.interfaceSegmentX;
    else
        coord = cutCell.interfaceSegmentY;
    end
    for i = 1:numel(cutCell.interfaceSegmentLength)
        binId = pd_network_locate_bin(coord(i), edges);
        if binId > 0
            interfaceLength(binId) = interfaceLength(binId) + cutCell.interfaceSegmentLength(i);
        end
    end
end

function [poreAreaSub, matrixAreaSub] = extractSliceCutCellAreas(cutCell, axisName, idxAxis)
    if strcmp(axisName, 'x')
        poreAreaSub = cutCell.pore.areaGrid(:, idxAxis);
        matrixAreaSub = cutCell.matrix.areaGrid(:, idxAxis);
    else
        poreAreaSub = cutCell.pore.areaGrid(idxAxis, :);
        matrixAreaSub = cutCell.matrix.areaGrid(idxAxis, :);
    end
end

function y = safeDivide(a, b)
    if isempty(b) || ~isfinite(b) || b == 0
        y = NaN;
    else
        y = a ./ b;
    end
end
