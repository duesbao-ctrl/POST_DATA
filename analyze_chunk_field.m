function out = analyze_chunk_field(chunkFile, varargin)
%ANALYZE_CHUNK_FIELD Analyze one variable on 1D/2D chunk data at selected timestep.
%   out = ANALYZE_CHUNK_FIELD(chunkFile, ...)
%
% Core behavior:
% 1) If variable exists in raw chunk data, plot it directly.
% 2) If not, allow derived variables:
%    T, vx, vy, vz, velocity, pressure, density,
%    Sxx, Syy, Szz, Sxy, Sxz, Syz, vonMisesS
%    (all based on compute_temp_stress_chunk results).
% 3) Otherwise throw error.

    p = inputParser;
    p.addRequired('chunkFile', @isTextScalar);
    p.addParameter('ChunkDim', '2d', @isTextScalar);   % '1d' or '2d'
    p.addParameter('Variable', 'c_rho', @isTextScalar);

    p.addParameter('SelectBy', 'Index', @isTextScalar);
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);

    p.addParameter('dV', [], @isnumeric);        % required for stress/vonMises
    p.addParameter('DoPlot', true, @islogical);
    p.addParameter('PlotOptions', {}, @iscell);  % passed to plot_cloud2d/plot_line1d
    p.parse(chunkFile, varargin{:});
    opt = p.Results;
    chunkFile = toChar(chunkFile);
    opt.ChunkDim = toChar(opt.ChunkDim);
    opt.Variable = toChar(opt.Variable);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex, ...
                    'ProgressMode', opt.ProgressMode};

    step = read_chunk_step_fast(chunkFile, selectorArgs{:});
    D = step.data;
    col = step.colIndex;

    varReqRaw = strtrim(opt.Variable);
    varReq = matlab.lang.makeValidName(varReqRaw);

    sourceType = 'raw';
    computedPack = [];

    if isfield(col, varReq)
        z = D(:, col.(varReq));
        yLabel = varReqRaw;
    else
        sourceType = 'derived';
        [z, yLabel, computedPack] = computeDerived(varReqRaw, chunkFile, selectorArgs, opt.dV);
    end

    [x, y, coordOk, coordMsg] = extractCoordinates(D, col, opt.ChunkDim);

    fig = [];
    ax = [];
    plotOut = [];
    if opt.DoPlot
        if ~coordOk
            warning('analyze_chunk_field:SkipPlotMissingCoord', '%s', coordMsg);
        else
            ttl = sprintf('%s @ timestep %g', yLabel, step.timestep);
            dimTag = lower(strtrim(opt.ChunkDim));

            if strcmp(dimTag, '1d')
                validPlot = isfinite(x) & isfinite(z);
                if any(validPlot)
                    [fig, ax, plotOut] = plot_line1d(x, z, ...
                        'Title', ttl, 'XLabel', 'x', 'YLabel', yLabel, opt.PlotOptions{:});
                else
                    warning('analyze_chunk_field:SkipPlotNoFiniteData', ...
                        'No finite x/z data found. Skip 1D plotting.');
                end
            elseif strcmp(dimTag, '2d')
                validPlot = isfinite(x) & isfinite(y) & isfinite(z);
                if any(validPlot)
                    [fig, ax, plotOut] = plot_cloud2d(x, y, z, ...
                        'Title', ttl, 'XLabel', 'x', 'YLabel', 'y', 'ColorbarLabel', yLabel, opt.PlotOptions{:});
                else
                    warning('analyze_chunk_field:SkipPlotNoFiniteData', ...
                        'No finite x/y/z data found. Skip 2D plotting.');
                end
            else
                error('analyze_chunk_field:BadChunkDim', 'ChunkDim must be ''1d'' or ''2d''.');
            end
        end
    end

    out = struct();
    out.filePath = chunkFile;
    out.selection = struct('SelectBy', opt.SelectBy, 'Index', opt.Index, ...
        'TimeStep', opt.TimeStep, 'Time', opt.Time);
    out.stepIndex = step.stepIndex;
    out.timestep = step.timestep;
    out.variableRequested = varReqRaw;
    out.variableUsed = yLabel;
    out.sourceType = sourceType;
    out.x = x;
    out.y = y;
    out.value = z;
    out.derivedSource = computedPack;
    out.figure = fig;
    out.axes = ax;
    out.plotOut = plotOut;
end

function [z, yLabel, pack] = computeDerived(varReqRaw, chunkFile, selectorArgs, dV)
    v = lower(strtrim(varReqRaw));
    pack = [];

    isT = strcmp(v, 't');
    isVelocity = any(strcmp(v, {'vx','vy','vz','velocity','speed'}));
    isPressure = any(strcmp(v, {'pressure','p'}));
    isDensity = any(strcmp(v, {'density','rho'}));
    isStress = any(strcmp(v, {'sxx','syy','szz','sxy','sxz','syz','vonmisess'}));

    if ~(isT || isVelocity || isPressure || isDensity || isStress)
        error('analyze_chunk_field:UnknownVariable', ...
            ['Variable "%s" not found in raw data and not in supported derived set: ', ...
             'T,vx,vy,vz,velocity,pressure,density,Sxx,Syy,Szz,Sxy,Sxz,Syz,vonMisesS.'], ...
             varReqRaw);
    end

    needsDV = isStress || isPressure || isDensity;
    if needsDV
        if isempty(dV) || ~isscalar(dV) || ~isfinite(dV) || dV <= 0
            error('analyze_chunk_field:NeedDV', ...
                'Derived variable "%s" requires a positive scalar dV.', varReqRaw);
        end
        c = compute_temp_stress_chunk(chunkFile, selectorArgs{:}, 'dV', dV);
    else
        % Velocity and T do not strictly need dV. Pass a valid placeholder.
        c = compute_temp_stress_chunk(chunkFile, selectorArgs{:}, 'dV', 1.0);
    end

    switch v
        case 't'
            z = c.T;
            yLabel = 'T';

        case 'vx'
            z = c.vx;
            yLabel = 'Vx';

        case 'vy'
            z = c.vy;
            yLabel = 'Vy';

        case 'vz'
            z = c.vz;
            yLabel = 'Vz';

        case {'velocity','speed'}
            z = c.speed;
            yLabel = 'Velocity';

        case {'pressure','p'}
            assertPressureComputed(c, 'Pressure');
            z = c.pressure;
            yLabel = 'Pressure (GPa)';

        case {'density','rho'}
            z = c.density;
            yLabel = 'Density (g/cm^3)';

        case 'sxx'
            assertComputed(c, 'xx', 'Sxx');
            z = c.sigma.xx;
            yLabel = 'Sxx';

        case 'syy'
            assertComputed(c, 'yy', 'Syy');
            z = c.sigma.yy;
            yLabel = 'Syy';

        case 'szz'
            assertComputed(c, 'zz', 'Szz');
            z = c.sigma.zz;
            yLabel = 'Szz';

        case 'sxy'
            assertComputed(c, 'xy', 'Sxy');
            z = c.sigma.xy;
            yLabel = 'Sxy';

        case 'sxz'
            assertComputed(c, 'xz', 'Sxz');
            z = c.sigma.xz;
            yLabel = 'Sxz';

        case 'syz'
            assertComputed(c, 'yz', 'Syz');
            z = c.sigma.yz;
            yLabel = 'Syz';

        case 'vonmisess'
            need = {'xx','yy','zz','xy','xz','yz'};
            for i = 1:numel(need)
                assertComputed(c, need{i}, ['vonMisesS (', need{i}, ')']);
            end
            sx = c.sigma.xx; sy = c.sigma.yy; sz = c.sigma.zz;
            txy = c.sigma.xy; txz = c.sigma.xz; tyz = c.sigma.yz;
            z = sqrt(0.5 .* ((sx-sy).^2 + (sy-sz).^2 + (sz-sx).^2 + 6.*(txy.^2 + txz.^2 + tyz.^2)));
            yLabel = 'vonMisesS';

        otherwise
            error('analyze_chunk_field:InternalVarSwitch', 'Unhandled derived variable: %s', varReqRaw);
    end

    pack = c;
end

function assertComputed(c, key, alias)
    if ~isfield(c, 'computed') || ~isfield(c.computed, key) || ~c.computed.(key)
        error('analyze_chunk_field:MissingStressComponent', ...
            'Cannot compute %s because required stress component "%s" is unavailable.', alias, key);
    end
end

function assertPressureComputed(c, alias)
    if ~isfield(c, 'computed') || ~isfield(c.computed, 'pressure') || ~c.computed.pressure
        error('analyze_chunk_field:MissingPressure', ...
            'Cannot compute %s because required diagonal stress components are unavailable.', alias);
    end
end

function [x, y, ok, msg] = extractCoordinates(D, col, chunkDim)
    ok = true;
    msg = '';
    if isfield(col, 'Coord1')
        x = D(:, col.Coord1);
    elseif isfield(col, 'c_x')
        x = D(:, col.c_x);
    elseif isfield(col, 'Chunk')
        x = D(:, col.Chunk);
    else
        x = [];
        y = [];
        ok = false;
        msg = 'Missing x column (Coord1/c_x/Chunk). Skip plotting.';
        return;
    end

    dimTag = lower(strtrim(chunkDim));
    if strcmp(dimTag, '2d')
        if isfield(col, 'Coord2')
            y = D(:, col.Coord2);
        elseif isfield(col, 'c_y')
            y = D(:, col.c_y);
        else
            y = [];
            ok = false;
            msg = '2D mode requires y column (Coord2/c_y). Skip plotting.';
        end
    else
        y = [];
    end
end

function tf = isTextScalar(v)
    tf = ischar(v) || (isstring(v) && isscalar(v));
end

function s = toChar(v)
    if isstring(v)
        s = char(v);
    else
        s = v;
    end
end
