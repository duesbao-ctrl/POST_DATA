function cutCell = pd_network_build_geometry(ncountGrid, validMask, poreMask, matrixMask, ...
        xCenters, yCenters, dx, dy, opt)
    cellArea = dx * dy;
    [ny, nx] = size(validMask);

    cutCell = struct();
    cutCell.geometryMode = opt.GeometryMode;
    cutCell.enabled = strcmp(opt.GeometryMode, 'cutcell');
    if cutCell.enabled
        cutCell.method = opt.CutCellMethod;
    else
        cutCell.method = 'none';
    end
    cutCell.fallback = opt.CutCellFallback;
    cutCell.plotRefinement = max(1, round(opt.CutCellPlotRefinement));
    cutCell.note = '';
    cutCell.fallbackCellCount = 0;
    cutCell.interfaceLength = 0;
    cutCell.interfaceSegmentX = zeros(0, 1);
    cutCell.interfaceSegmentY = zeros(0, 1);
    cutCell.interfaceSegmentLength = zeros(0, 1);

    poreFraction = zeros(ny, nx);
    poreFraction(poreMask) = 1;
    poreFraction(~validMask) = NaN;
    matrixFraction = zeros(ny, nx);
    matrixFraction(matrixMask) = 1;
    matrixFraction(~validMask) = NaN;

    poreCentroidX = nan(ny, nx);
    poreCentroidY = nan(ny, nx);
    matrixCentroidX = nan(ny, nx);
    matrixCentroidY = nan(ny, nx);
    for r = 1:ny
        for c = 1:nx
            if validMask(r, c)
                poreCentroidX(r, c) = xCenters(c);
                poreCentroidY(r, c) = yCenters(r);
                matrixCentroidX(r, c) = xCenters(c);
                matrixCentroidY(r, c) = yCenters(r);
            end
        end
    end

    interfaceLengthGrid = zeros(ny, nx);

    if cutCell.enabled
        phiCenter = opt.ThresholdN - ncountGrid;
        cornerPhi = estimateCornerPhi(phiCenter, validMask, opt.Boundary);
        poreFraction(:) = NaN;
        matrixFraction(:) = NaN;
        poreCentroidX(:) = NaN;
        poreCentroidY(:) = NaN;
        matrixCentroidX(:) = NaN;
        matrixCentroidY(:) = NaN;

        segX = zeros(nnz(validMask) * 4, 1);
        segY = zeros(nnz(validMask) * 4, 1);
        segLen = zeros(nnz(validMask) * 4, 1);
        segCount = 0;

        for r = 1:ny
            for c = 1:nx
                if ~validMask(r, c)
                    continue;
                end
                x0 = xCenters(c) - 0.5 * dx;
                x1 = xCenters(c) + 0.5 * dx;
                y0 = yCenters(r) - 0.5 * dy;
                y1 = yCenters(r) + 0.5 * dy;
                vertexXY = [ ...
                    x0, y0; ...
                    x1, y0; ...
                    x1, y1; ...
                    x0, y1];
                vertexPhi = [ ...
                    cornerPhi(r, c); ...
                    cornerPhi(r, c + 1); ...
                    cornerPhi(r + 1, c + 1); ...
                    cornerPhi(r + 1, c)];
                centerXY = [xCenters(c), yCenters(r)];
                centerPhi = phiCenter(r, c);

                if any(~isfinite(vertexPhi)) || ~isfinite(centerPhi)
                    if strcmp(cutCell.fallback, 'error')
                        error('analyze_chunk_network2d:CutCellFallbackRequired', ...
                            'Cut-cell reconstruction needs finite center/corner values for every valid cell.');
                    end
                    cutCell.fallbackCellCount = cutCell.fallbackCellCount + 1;
                    poreFraction(r, c) = double(poreMask(r, c));
                    matrixFraction(r, c) = double(matrixMask(r, c));
                    poreCentroidX(r, c) = xCenters(c);
                    poreCentroidY(r, c) = yCenters(r);
                    matrixCentroidX(r, c) = xCenters(c);
                    matrixCentroidY(r, c) = yCenters(r);
                    continue;
                end

                [pArea, pCx, pCy, mArea, mCx, mCy, iLen, iX, iY] = ...
                    cutCellFromCenterCornerTriangulation(vertexXY, vertexPhi, centerXY, centerPhi);
                poreFraction(r, c) = clamp01(pArea / cellArea);
                matrixFraction(r, c) = clamp01(mArea / cellArea);
                if pArea > 0
                    poreCentroidX(r, c) = pCx;
                    poreCentroidY(r, c) = pCy;
                end
                if mArea > 0
                    matrixCentroidX(r, c) = mCx;
                    matrixCentroidY(r, c) = mCy;
                end
                interfaceLengthGrid(r, c) = iLen;
                if iLen > 0
                    segCount = segCount + 1;
                    if segCount > numel(segLen)
                        segX = [segX; zeros(numel(segX), 1)]; %#ok<AGROW>
                        segY = [segY; zeros(numel(segY), 1)]; %#ok<AGROW>
                        segLen = [segLen; zeros(numel(segLen), 1)]; %#ok<AGROW>
                    end
                    segX(segCount) = iX;
                    segY(segCount) = iY;
                    segLen(segCount) = iLen;
                end
            end
        end

        cutCell.interfaceSegmentX = segX(1:segCount);
        cutCell.interfaceSegmentY = segY(1:segCount);
        cutCell.interfaceSegmentLength = segLen(1:segCount);
        cutCell.interfaceLength = sum(cutCell.interfaceSegmentLength);
        cutCell.note = ['Cut-cell geometry uses piecewise-linear interpolation of ', ...
            'ThresholdN-Ncount over center-to-corner triangles. It changes geometric ', ...
            'area/interface/diameter/profile estimates, not binary topology/connectivity.'];
    else
        cutCell.interfaceLength = computeBinaryInterfaceLength(poreMask, validMask, dx, dy, opt.Boundary);
        cutCell.note = 'Cut-cell geometry is disabled; binary cell areas and grid-edge interface lengths are used.';
    end

    [poreFraction, matrixFraction] = normalizePhaseFractions(poreFraction, matrixFraction, validMask);

    cutCell.phaseGrid = poreFraction;
    cutCell.phaseGrid(~validMask) = NaN;
    cutCell.pore = buildCutPhaseGeometry('pore', poreFraction, poreCentroidX, poreCentroidY, ...
        interfaceLengthGrid, cellArea, cutCell.enabled, cutCell.note);
    cutCell.matrix = buildCutPhaseGeometry('matrix', matrixFraction, matrixCentroidX, matrixCentroidY, ...
        interfaceLengthGrid, cellArea, cutCell.enabled, cutCell.note);
end

function [poreFraction, matrixFraction] = normalizePhaseFractions(poreFraction, matrixFraction, validMask)
    idx = find(validMask);
    for i = 1:numel(idx)
        lin = idx(i);
        pf = poreFraction(lin);
        mf = matrixFraction(lin);
        if ~(isfinite(pf) && isfinite(mf))
            continue;
        end
        pf = clamp01(pf);
        mf = clamp01(mf);
        s = pf + mf;
        if s > 0
            pf = pf / s;
            mf = mf / s;
        end
        poreFraction(lin) = pf;
        matrixFraction(lin) = mf;
    end
end
function phaseGeom = buildCutPhaseGeometry(name, fraction, centroidX, centroidY, interfaceLengthGrid, ...
        cellArea, enabled, note)
    areaGrid = fraction .* cellArea;
    areaGrid(~isfinite(areaGrid)) = 0;
    phaseGeom = struct();
    phaseGeom.name = name;
    phaseGeom.enabled = enabled;
    phaseGeom.note = note;
    phaseGeom.fraction = fraction;
    phaseGeom.areaGrid = areaGrid;
    phaseGeom.centroidXGrid = centroidX;
    phaseGeom.centroidYGrid = centroidY;
    phaseGeom.interfaceLengthGrid = interfaceLengthGrid;
end

function cornerPhi = estimateCornerPhi(phiCenter, validMask, boundary)
    [ny, nx] = size(phiCenter);
    cornerPhi = nan(ny + 1, nx + 1);
    for rr = 1:(ny + 1)
        for cc = 1:(nx + 1)
            cells = adjacentCellsForCorner(rr, cc, ny, nx, boundary);
            vals = nan(size(cells, 1), 1);
            n = 0;
            for k = 1:size(cells, 1)
                r = cells(k, 1);
                c = cells(k, 2);
                if r >= 1 && r <= ny && c >= 1 && c <= nx && validMask(r, c)
                    n = n + 1;
                    vals(n) = phiCenter(r, c);
                end
            end
            if n > 0
                cornerPhi(rr, cc) = mean(vals(1:n));
            end
        end
    end
end
function cells = adjacentCellsForCorner(rr, cc, ny, nx, boundary)
    rows = [rr - 1, rr];
    cols = [cc - 1, cc];
    cells = zeros(4, 2);
    n = 0;
    for i = 1:2
        r = rows(i);
        if r < 1
            if pd_network_is_periodic_axis(boundary, 'y')
                r = ny;
            else
                continue;
            end
        elseif r > ny
            if pd_network_is_periodic_axis(boundary, 'y')
                r = 1;
            else
                continue;
            end
        end
        for j = 1:2
            c = cols(j);
            if c < 1
                if pd_network_is_periodic_axis(boundary, 'x')
                    c = nx;
                else
                    continue;
                end
            elseif c > nx
                if pd_network_is_periodic_axis(boundary, 'x')
                    c = 1;
                else
                    continue;
                end
            end
            n = n + 1;
            cells(n, :) = [r, c];
        end
    end
    cells = cells(1:n, :);
end

function [pArea, pCx, pCy, mArea, mCx, mCy, iLen, iX, iY] = ...
        cutCellFromCenterCornerTriangulation(vertexXY, vertexPhi, centerXY, centerPhi)
    pArea = 0;
    pMx = 0;
    pMy = 0;
    mArea = 0;
    mMx = 0;
    mMy = 0;
    iLen = 0;
    iMx = 0;
    iMy = 0;

    triCorner = [1 2; 2 3; 3 4; 4 1];
    for i = 1:4
        pts = [centerXY; vertexXY(triCorner(i, 1), :); vertexXY(triCorner(i, 2), :)];
        vals = [centerPhi; vertexPhi(triCorner(i, 1)); vertexPhi(triCorner(i, 2))];
        [areaPos, cxPos, cyPos] = clippedTriangleMoment(pts, vals, true);
        [areaNeg, cxNeg, cyNeg] = clippedTriangleMoment(pts, vals, false);
        [segLen, segX, segY] = triangleInterfaceSegment(pts, vals);

        pArea = pArea + areaPos;
        pMx = pMx + areaPos * cxPos;
        pMy = pMy + areaPos * cyPos;
        mArea = mArea + areaNeg;
        mMx = mMx + areaNeg * cxNeg;
        mMy = mMy + areaNeg * cyNeg;
        iLen = iLen + segLen;
        iMx = iMx + segLen * segX;
        iMy = iMy + segLen * segY;
    end

    pCx = safeMomentCenter(pMx, pArea, centerXY(1));
    pCy = safeMomentCenter(pMy, pArea, centerXY(2));
    mCx = safeMomentCenter(mMx, mArea, centerXY(1));
    mCy = safeMomentCenter(mMy, mArea, centerXY(2));
    iX = safeMomentCenter(iMx, iLen, centerXY(1));
    iY = safeMomentCenter(iMy, iLen, centerXY(2));
end

function [area, cx, cy] = clippedTriangleMoment(pts, vals, keepPositive)
    polyPts = pts;
    polyVals = vals(:);
    [polyPts, ~] = clipScalarPolygon(polyPts, polyVals, keepPositive);
    [area, cx, cy] = polygonAreaCentroid(polyPts);
end

function [outPts, outVals] = clipScalarPolygon(inPts, inVals, keepPositive)
    outPts = zeros(0, 2);
    outVals = zeros(0, 1);
    n = size(inPts, 1);
    if n == 0
        return;
    end

    for i = 1:n
        j = i + 1;
        if j > n
            j = 1;
        end
        p1 = inPts(i, :);
        p2 = inPts(j, :);
        v1 = inVals(i);
        v2 = inVals(j);
        in1 = scalarInside(v1, keepPositive);
        in2 = scalarInside(v2, keepPositive);

        if in1
            outPts(end + 1, :) = p1; %#ok<AGROW>
            outVals(end + 1, 1) = v1; %#ok<AGROW>
        end
        if xor(in1, in2)
            t = v1 / (v1 - v2);
            t = min(max(t, 0), 1);
            p = p1 + t .* (p2 - p1);
            outPts(end + 1, :) = p; %#ok<AGROW>
            outVals(end + 1, 1) = 0; %#ok<AGROW>
        end
    end
end

function inside = scalarInside(v, keepPositive)
    tol = 1e-12;
    if keepPositive
        inside = (v > tol);
    else
        inside = (v <= tol);
    end
end

function [area, cx, cy] = polygonAreaCentroid(pts)
    area = 0;
    cx = NaN;
    cy = NaN;
    n = size(pts, 1);
    if n < 3
        return;
    end

    x = pts(:, 1);
    y = pts(:, 2);
    x2 = x([2:end, 1]);
    y2 = y([2:end, 1]);
    crossVal = x .* y2 - x2 .* y;
    signedArea = 0.5 * sum(crossVal);
    if abs(signedArea) <= eps(max(max(abs(x)), max(abs(y))) + 1)
        return;
    end
    cxSigned = sum((x + x2) .* crossVal) / (6 * signedArea);
    cySigned = sum((y + y2) .* crossVal) / (6 * signedArea);
    area = abs(signedArea);
    cx = cxSigned;
    cy = cySigned;
end

function [segLen, segX, segY] = triangleInterfaceSegment(pts, vals)
    crossPts = zeros(3, 2);
    nCross = 0;
    edgePair = [1 2; 2 3; 3 1];
    for e = 1:3
        i = edgePair(e, 1);
        j = edgePair(e, 2);
        v1 = vals(i);
        v2 = vals(j);
        if (v1 > 0 && v2 <= 0) || (v1 <= 0 && v2 > 0)
            t = v1 / (v1 - v2);
            t = min(max(t, 0), 1);
            nCross = nCross + 1;
            crossPts(nCross, :) = pts(i, :) + t .* (pts(j, :) - pts(i, :));
        end
    end

    if nCross < 2
        segLen = 0;
        segX = NaN;
        segY = NaN;
        return;
    end

    p1 = crossPts(1, :);
    p2 = crossPts(2, :);
    segLen = hypot(p2(1) - p1(1), p2(2) - p1(2));
    segX = 0.5 * (p1(1) + p2(1));
    segY = 0.5 * (p1(2) + p2(2));
end

function value = safeMomentCenter(moment, weight, fallback)
    if weight > 0 && isfinite(moment)
        value = moment / weight;
    else
        value = fallback;
    end
end

function value = clamp01(value)
    value = min(max(value, 0), 1);
end

function len = computeBinaryInterfaceLength(poreMask, validMask, dx, dy, boundary)
    [~, len] = computePerimeterForMask(poreMask, poreMask, validMask, dx, dy, boundary);
end


function [perimeterOpen, interfacePerimeter] = computePerimeterForMask(componentMask, phaseMask, validMask, dx, dy, boundary)
    [ny, nx] = size(componentMask);
    perimeterOpen = 0;
    interfacePerimeter = 0;
    [rows, cols] = find(componentMask);
    for i = 1:numel(rows)
        r = rows(i);
        c = cols(i);
        for direction = 1:4
            [nr, nc, hasNeighbor] = pd_network_neighbor(r, c, direction, ny, nx, boundary);
            if direction == 1 || direction == 2, sideLength = dy; else, sideLength = dx; end
            if ~hasNeighbor || ~validMask(nr, nc)
                perimeterOpen = perimeterOpen + sideLength;
            elseif ~phaseMask(nr, nc)
                perimeterOpen = perimeterOpen + sideLength;
                interfacePerimeter = interfacePerimeter + sideLength;
            end
        end
    end
end
