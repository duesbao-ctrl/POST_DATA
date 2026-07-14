function [dx, dy, source] = pd_network_resolve_spacing(options, chunkFile, xCenters, yCenters)
%PD_NETWORK_RESOLVE_SPACING Resolve grid spacing from options, filename, or coordinates.

    meta = parseFilenameMetadata(chunkFile);
    [dx, srcDx] = resolveOneSpacing(options.Dx, 'dx', xCenters, meta, options.CoordScale);
    [dy, srcDy] = resolveOneSpacing(options.Dy, 'dy', yCenters, meta, options.CoordScale);
    source = struct('dx', srcDx, 'dy', srcDy);
end

function [spacing, source] = resolveOneSpacing(userValue, key, coords, meta, coordScale)
    if ~isempty(userValue)
        spacing = userValue .* coordScale;
        source = 'option';
        return;
    end
    if isfield(meta, key)
        spacing = meta.(key) * coordScale;
        if isfinite(spacing) && spacing > 0
            source = 'filename';
            return;
        end
    end
    spacing = inferSpacingFromCoords(coords, key);
    source = 'coordinates';
end

function spacing = inferSpacingFromCoords(coords, key)
    coords = coords(:);
    if numel(coords) < 2
        error('analyze_chunk_network2d:NeedSpacing', ...
            'Cannot infer %s from coordinates with fewer than 2 unique centers.', key);
    end
    differences = diff(sort(coords));
    if any(~isfinite(differences)) || any(differences <= 0)
        error('analyze_chunk_network2d:BadCoordinateSpacing', ...
            'Cannot infer %s from non-monotone coordinates.', key);
    end
    spacing = median(differences);
    tolerance = max(1e-9, 1e-2 * max(abs(spacing), 1));
    if any(abs(differences - spacing) > tolerance)
        error('analyze_chunk_network2d:IrregularGrid', ...
            'Cannot infer %s from irregular coordinate spacing.', key);
    end
end

function meta = parseFilenameMetadata(chunkFile)
    [~, nameOnly, ~] = fileparts(chunkFile);
    parts = strsplit(nameOnly, '_');
    meta = struct();
    for i = 2:(numel(parts) - 1)
        key = matlab.lang.makeValidName(parts{i});
        value = str2double(parts{i + 1});
        if ~isnan(value), meta.(key) = value; end
    end
end
