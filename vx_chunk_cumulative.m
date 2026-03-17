function out = vx_chunk_cumulative(vxChunkPath, varargin)
%VX_CHUNK_CUMULATIVE Cumulative areal-density vs velocity for vx_chunk data.
%   Supports SelectBy = Index | TimeStep | Time.
%
% Example:
% out = vx_chunk_cumulative('d:\...\vx_chunk.txt', ...
%     'SelectBy','Time','Time',33.0,'SlurmPath','d:\...\slurm_9.log');

    p = inputParser;
    p.addRequired('vxChunkPath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);      % Index | TimeStep | Time
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);

    p.addParameter('VelocityFactor', 0.001, @isnumeric);  % v = Chunk * VelocityFactor
    p.addParameter('DensityVars', {}, @(x) iscell(x) || isTextScalar(x));
    p.addParameter('MakePlot', true, @islogical);
    p.parse(vxChunkPath, varargin{:});
    opt = p.Results;
    vxChunkPath = toChar(vxChunkPath);
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

    S = read_chunk_step_fast(vxChunkPath, selectorArgs{:});
    stepIdx = S.stepIndex;
    selectionInfo = buildSelectionInfo(S, opt);
    col = S.colIndex;

    if ~isfield(col, 'Chunk')
        error('vx_chunk_cumulative:MissingChunk', 'Column "Chunk" is required.');
    end

    densityVarsReq = resolveDensityVars(S, opt.DensityVars);

    v = S.data(:, col.Chunk) .* opt.VelocityFactor;
    densAll = zeros(size(S.data,1), numel(densityVarsReq));
    keep = false(1, numel(densityVarsReq));
    for k = 1:numel(densityVarsReq)
        vn = densityVarsReq{k};
        if ~isfield(col, vn)
            warning('vx_chunk_cumulative:MissingDensityVar', ...
                'Density variable "%s" not found in this file. Skip this curve.', vn);
            continue;
        end
        densAll(:, k) = S.data(:, col.(vn));
        keep(k) = true;
    end

    densityVars = densityVarsReq(keep);
    dens = densAll(:, keep);
    nGroup = numel(densityVars);
    if nGroup == 0
        warning('vx_chunk_cumulative:NoAvailableDensityVar', ...
            'No available ArealDensity columns to process. Nothing will be plotted.');
        [vSortedOnly, ~] = sort(v, 'ascend');
        out = struct();
        out.filePath = vxChunkPath;
        out.selection = selectionInfo;
        out.stepIndex = stepIdx;
        out.timestep = S.timestep;
        out.velocity = vSortedOnly;
        out.density = [];
        out.cumulativeDensity = [];
        out.densityVars = {};
        out.plottedVars = {};
        out.figure = [];
        return;
    end

    % Enforce: massArealDensity = mass1ArealDensity + mass2ArealDensity (if all exist).
    [idxTotal, idx1, idx2] = findMassDensityIndices(densityVars);
    if idxTotal > 0 && idx1 > 0 && idx2 > 0
        dens(:, idxTotal) = dens(:, idx1) + dens(:, idx2);
    end

    [vSorted, idx] = sort(v, 'ascend');
    densSorted = dens(idx, :);

    % Cumulative from max velocity to min velocity:
    % cum(i,:) = sum of density for all j with v(j) >= v(i)
    cumDesc = flipud(cumsum(flipud(densSorted), 1));

    fig = [];
    if opt.MakePlot
        fig = figure('Name', 'Cumulative Areal Density vs Velocity');
        hold on;

        % First pass: check which curves are non-zero.
        isNonZero = false(1, nGroup);
        for k = 1:nGroup
            yk = cumDesc(:, k);
            ymax = max(yk(isfinite(yk)));
            isNonZero(k) = ~(isempty(ymax) || ymax <= 0);
        end

        % Rule requested by user:
        % if only one of mass1/mass2 is non-zero, plot only massArealDensity.
        plotMask = isNonZero;
        if idxTotal > 0 && idx1 > 0 && idx2 > 0
            oneSubNonZero = xor(isNonZero(idx1), isNonZero(idx2));
            if oneSubNonZero && isNonZero(idxTotal)
                plotMask(:) = false;
                plotMask(idxTotal) = true;
            end
        end

        hList = [];
        nameList = {};
        for k = 1:nGroup
            if ~plotMask(k)
                continue;
            end
            yk = cumDesc(:, k);
            hk = plot(vSorted, yk, 'LineWidth', 1.6);
            hList = [hList; hk]; %#ok<AGROW>
            nameList{end+1} = densityVars{k}; %#ok<AGROW>
        end
        xlabel('Velocity (km/s)');
        ylabel('Cumulative Areal Density (mg/cm^2)');
        title(sprintf('vx\\_chunk cumulative distribution @ timestep %g', S.timestep));
        if ~isempty(hList)
            legend(hList, nameList, 'Location', 'best', 'Interpreter', 'none');
        else
            warning('vx_chunk_cumulative:NoNonZeroCurve', ...
                'All cumulative curves are zero. Skip drawing data curves.');
            text(0.5, 0.5, 'No non-zero cumulative curves to display', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
        grid on;
    end

    out = struct();
    out.filePath = vxChunkPath;
    out.selection = selectionInfo;
    out.stepIndex = stepIdx;
    out.timestep = S.timestep;
    out.velocity = vSorted;
    out.density = densSorted;
    out.cumulativeDensity = cumDesc;
    out.densityVars = densityVars;
    if exist('nameList', 'var')
        out.plottedVars = nameList;
    else
        out.plottedVars = {};
    end
    out.figure = fig;
end

function densityVars = resolveDensityVars(vxData, densityVarsInput)
    if isstring(densityVarsInput)
        densityVarsInput = cellstr(densityVarsInput);
    elseif ischar(densityVarsInput)
        densityVarsInput = {densityVarsInput};
    end

    if ~isempty(densityVarsInput)
        densityVars = cellfun(@matlab.lang.makeValidName, densityVarsInput, 'UniformOutput', false);
        return;
    end

    % Auto-detect all areal-density columns, keep original order.
    vn = vxData.validVarNames;
    densityVars = {};
    for i = 1:numel(vn)
        name = vn{i};
        if ~isempty(strfind(name, 'ArealDensity'))
            densityVars{end+1} = name; %#ok<AGROW>
        end
    end

    if isempty(densityVars)
        error('vx_chunk_cumulative:NoDensityCols', 'No ArealDensity columns found.');
    end
end

function info = buildSelectionInfo(step, opt)
    mode = lower(strtrim(opt.SelectBy));
    info = struct();
    info.mode = mode;
    info.stepIndex = step.stepIndex;
    info.mappedTimeStep = step.timestep;

    switch mode
        case 'index'
            info.requestedIndex = opt.Index;

        case 'timestep'
            info.requestedTimeStep = opt.TimeStep;

        case 'time'
            info.requestedTime = opt.Time;

        otherwise
            % read_chunk_step_fast has already validated selector.
    end
end

function [idxTotal, idx1, idx2] = findMassDensityIndices(densityVars)
    idxTotal = 0;
    idx1 = 0;
    idx2 = 0;
    for i = 1:numel(densityVars)
        n = densityVars{i};
        if ~isempty(strfind(n, 'mass1ArealDensity'))
            idx1 = i;
        elseif ~isempty(strfind(n, 'mass2ArealDensity'))
            idx2 = i;
        elseif ~isempty(strfind(n, 'massArealDensity'))
            idxTotal = i;
        end
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
