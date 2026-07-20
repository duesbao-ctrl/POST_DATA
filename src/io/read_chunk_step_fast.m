function out = read_chunk_step_fast(filePath, varargin)
%READ_CHUNK_STEP_FAST Fast read one timestep block from chunk text file.
%   out = READ_CHUNK_STEP_FAST(filePath, 'SelectBy','TimeStep','TimeStep',5000)
%   out = READ_CHUNK_STEP_FAST(filePath, 'SelectBy','Index','Index',12)
%   out = READ_CHUNK_STEP_FAST(filePath, 'SelectBy','Time','Time',33.0,'SlurmPath','slurm_9.log')
%   out = READ_CHUNK_STEP_FAST(filePath, 'ProgressMode','console')
%
% Key idea:
%   Build/reuse a lightweight index of block header positions, then seek to
%   target block and read only one block instead of the whole file.

    hasExplicitSelectBy = hasExplicitNameValue(varargin, 'SelectBy');

    p = inputParser;
    p.addRequired('filePath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);      % Index | TimeStep | Time
    p.addParameter('Index', [], @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);
    p.addParameter('CancelCallback', @() false, @(x) isa(x, 'function_handle'));
    p.addParameter('ProgressCallback', @(fraction, message) [], ...
        @(x) isa(x, 'function_handle'));
    p.parse(filePath, varargin{:});
    opt = p.Results;

    filePath = toChar(filePath);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);
    validateReaderOptions(opt);
    if ~exist(filePath, 'file')
        error('read_chunk_step_fast:FileNotFound', 'File not found: %s', filePath);
    end

    reportProgress(opt.CancelCallback, opt.ProgressCallback, 0, ...
        'Preparing chunk input');
    idx = buildOrLoadIndex(filePath, opt.ProgressMode, ...
        opt.CancelCallback, opt.ProgressCallback);
    nBlocks = numel(idx.timesteps);
    if nBlocks == 0
        error('read_chunk_step_fast:NoBlocks', 'No timestep blocks found in file: %s', filePath);
    end

    mode = lower(strtrim(opt.SelectBy));
    % Backward compatibility: infer mode when SelectBy keeps default but
    % caller passes only Time or TimeStep.
    if ~hasExplicitSelectBy && strcmp(mode, 'index')
        if ~isempty(opt.Time)
            mode = 'time';
        elseif ~isempty(opt.TimeStep) && isempty(opt.Index)
            mode = 'timestep';
        end
    end
    switch mode
        case 'index'
            if isempty(opt.Index)
                error('read_chunk_step_fast:MissingIndex', ...
                    'SelectBy=Index requires Index.');
            end
        k = round(opt.Index);
        if k < 1 || k > nBlocks
            error('read_chunk_step_fast:BadIndex', 'Index out of range: %d (1..%d)', k, nBlocks);
        end
        sel = k;
        reqStep = idx.timesteps(sel);
        case 'time'
            if isempty(opt.Time)
                error('read_chunk_step_fast:MissingTime', ...
                    'SelectBy=Time requires Time.');
            end
        validEmbeddedTime = isfinite(idx.physicalTimes);
        if any(validEmbeddedTime)
            candidates = find(validEmbeddedTime);
            [~, nearest] = min(abs(idx.physicalTimes(candidates) - opt.Time));
            sel = candidates(nearest);
            reqStep = idx.timesteps(sel);
        else
            slurmPath = opt.SlurmPath;
            if isempty(slurmPath)
                try
                    slurmPath = get_slurm_txt_fullpath(fileparts(filePath));
                catch discoveryError
                    if strcmp(discoveryError.identifier, ...
                            'get_slurm_txt_fullpath:NoMatch')
                        error('read_chunk_step_fast:MissingSlurmPath', ...
                            ['SlurmPath is required when selecting by Time ' ...
                             'from a legacy chunk file without "# Time" ' ...
                             'metadata.']);
                    end
                    rethrow(discoveryError);
                end
            end
            S = read_slurm_stepcpu(slurmPath);
            m = round(opt.SlurmModuleIndex);
            if m < 1 || m > numel(S.modules)
                error('read_chunk_step_fast:BadSlurmModule', ...
                    'SlurmModuleIndex out of range: %d', m);
            end
            scol = S.modules(m).colIndex;
            if ~isfield(scol, 'Time') || ~isfield(scol, 'Step')
                error('read_chunk_step_fast:SlurmMissingCols', ...
                    'Slurm block missing Time/Step columns.');
            end
            slurmTimeAll = S.modules(m).data(:, scol.Time);
            slurmStepAll = S.modules(m).data(:, scol.Step);
            validSlurm = isfinite(slurmTimeAll) & isfinite(slurmStepAll);
            if ~any(validSlurm)
                error('read_chunk_step_fast:SlurmNoFiniteRows', ...
                    'No finite Time/Step rows found in slurm block %d.', m);
            end
            slurmTime = slurmTimeAll(validSlurm);
            slurmStep = slurmStepAll(validSlurm);
            [~, it] = min(abs(slurmTime - opt.Time));
            req = slurmStep(it);
            [sel, reqStep] = locateByStep(idx, req);
        end
        case 'timestep'
            if isempty(opt.TimeStep)
                error('read_chunk_step_fast:MissingTimeStep', ...
                    'SelectBy=TimeStep requires TimeStep.');
            end
        req = opt.TimeStep;
        [sel, reqStep] = locateByStep(idx, req);
        otherwise
        error('read_chunk_step_fast:NeedSelector', ...
            'SelectBy must be Index/TimeStep/Time.');
    end

    fid = fopen(filePath, 'r');
    if fid < 0
        error('read_chunk_step_fast:FileOpenFailed', 'Cannot open file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    preamble = pd_read_chunk_preamble(fid, 'read_chunk_step_fast');
    varNames = preamble.varNames;
    validVarNames = preamble.validVarNames;
    numVars = numel(varNames);

    reportProgress(opt.CancelCallback, opt.ProgressCallback, 0.65, ...
        sprintf('Reading timestep block %d of %d', sel, nBlocks));
    seekStatus = fseek(fid, idx.blockPos(sel), 'bof');
    if seekStatus ~= 0
        error('read_chunk_step_fast:SeekFailed', ...
            'Cannot seek to block %d in file: %s', sel, filePath);
    end
    bh = fgetl(fid);
    h = pd_parse_chunk_summary_line(bh, 'read_chunk_step_fast');
    timestep = h(1);
    numChunks = round(h(2));
    totalCount = h(3);

    rowProgress = @(rowsRead, totalRows) reportProgress( ...
        opt.CancelCallback, opt.ProgressCallback, ...
        0.65 + 0.34 .* rowsRead ./ max(totalRows, 1), ...
        'Reading selected timestep data');
    data = pd_read_chunk_data_rows(fid, numChunks, numVars, ...
        'read_chunk_step_fast', timestep, rowProgress);

    colIndex = struct();
    for i = 1:numel(validVarNames)
        colIndex.(validVarNames{i}) = i;
    end

    data = applyEmptyCellNaN(data, validVarNames);

    out = struct();
    out.filePath = filePath;
    out.headerLines = preamble.headerLines;
    out.varNames = varNames;
    out.validVarNames = validVarNames;
    out.colIndex = colIndex;
    out.metadata = preamble.metadata;
    out.inputFormat = preamble.metadata.inputFormat;
    out.unitSystem = preamble.metadata.unitSystem;
    out.taskName = preamble.metadata.name;
    out.chunkKind = preamble.metadata.kind;
    out.requestedTimeStep = reqStep;
    out.stepIndex = sel;
    out.timestep = timestep;
    out.physicalTime = idx.physicalTimes(sel);
    out.numChunks = numChunks;
    out.totalCount = totalCount;
    out.data = data;
    out.index = idx;
    reportProgress(opt.CancelCallback, opt.ProgressCallback, 1, ...
        'Chunk input ready');
end

function data = applyEmptyCellNaN(data, varNames)
% For 1d/2d bin-style chunk files:
% if Ncount==0 in a cell, set non-coordinate/value columns to NaN.

    idxN = find(strcmp(varNames, 'Ncount'), 1, 'first');
    idxC1 = find(strcmp(varNames, 'Coord1'), 1, 'first');
    if isempty(idxN) || isempty(idxC1)
        return;
    end

    keep = false(1, numel(varNames));
    for k = 1:numel(varNames)
        vn = varNames{k};
        if strcmp(vn, 'Chunk') || strcmp(vn, 'Ncount') || strncmp(vn, 'Coord', 5)
            keep(k) = true;
        end
    end

    m = (data(:, idxN) == 0);
    if any(m) && any(~keep)
        data(m, ~keep) = NaN;
    end
end

function [sel, reqStep] = locateByStep(idx, req)
    nBlocks = numel(idx.timesteps);
    % Always use nearest search with uniform-step shortcut.
    if idx.isUniform
        k0 = round((req - idx.timesteps(1)) / idx.dt) + 1;
        k0 = max(1, min(nBlocks, k0));
        candidates = unique(max(1, min(nBlocks, [k0-1, k0, k0+1])));
        [~, ii] = min(abs(idx.timesteps(candidates) - req));
        sel = candidates(ii);
    else
        [~, sel] = min(abs(idx.timesteps - req));
    end
    reqStep = req;
end

function idx = buildOrLoadIndex(filePath, progressMode, cancelCallback, progressCallback)
    d = dir(filePath);
    cachePath = [filePath, '.stepidx.mat'];
    curSig = buildQuickSignature(filePath);

    if exist(cachePath, 'file')
        try
            S = load(cachePath, 'idx');
            if isfield(S, 'idx') && isfield(S.idx, 'fileSize') && ...
                    isfield(S.idx, 'fileDatenum')
                hasSig = isfield(S.idx, 'quickSig');
                hasCurrentFormat = isfield(S.idx, 'formatVersion') && ...
                    isequal(S.idx.formatVersion, chunkIndexFormatVersion()) && ...
                    isfield(S.idx, 'physicalTimes') && ...
                    isfield(S.idx, 'metadata');
                sameSig = hasSig && isSameQuickSig(S.idx.quickSig, curSig);
                if isequal(S.idx.fileSize, d.bytes) && ...
                        isequal(S.idx.fileDatenum, d.datenum) && sameSig && ...
                        hasCurrentFormat
                    idx = S.idx;
                    reportProgress(cancelCallback, progressCallback, 0.6, ...
                        'Validated chunk index cache');
                    return;
                end
            end
        catch cacheError
            warning('read_chunk_step_fast:InvalidIndexCache', ...
                'Ignoring invalid index cache %s (%s).', ...
                cachePath, cacheError.message);
        end
    end

    fid = fopen(filePath, 'r');
    if fid < 0
        error('read_chunk_step_fast:FileOpenFailed', 'Cannot open file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));
    tracker = make_file_progress(progressMode, 'Building index', filePath);
    progressCleanup = onCleanup(@() tracker.close()); %#ok<NASGU>
    totalBytes = max(d.bytes, 1);

    preamble = pd_read_chunk_preamble(fid, 'read_chunk_step_fast');
    tracker.update(ftell(fid) / totalBytes);

    blockPos = zeros(1000, 1);
    dataPos  = zeros(1000, 1);
    timesteps = zeros(1000, 1);
    physicalTimes = nan(1000, 1);
    numChunks = zeros(1000, 1);
    n = 0;

    while true
        frame = pd_read_next_chunk_frame_header(fid, 'read_chunk_step_fast');
        if frame.eof
            break;
        end

        n = n + 1;
        if n > numel(blockPos)
            blockPos = [blockPos; zeros(numel(blockPos),1)]; %#ok<AGROW>
            dataPos  = [dataPos;  zeros(numel(dataPos),1)]; %#ok<AGROW>
            timesteps = [timesteps; zeros(numel(timesteps),1)]; %#ok<AGROW>
            physicalTimes = [physicalTimes; nan(numel(physicalTimes),1)]; %#ok<AGROW>
            numChunks = [numChunks; zeros(numel(numChunks),1)]; %#ok<AGROW>
        end

        blockPos(n) = frame.blockPosition;
        dataPos(n) = frame.dataPosition;
        timesteps(n) = frame.timestep;
        physicalTimes(n) = frame.physicalTime;
        numChunks(n) = frame.numChunks;

        pd_read_chunk_data_rows(fid, frame.numChunks, ...
            numel(preamble.varNames), 'read_chunk_step_fast', ...
            frame.timestep);
        tracker.update(ftell(fid) / totalBytes);
        if mod(n, 16) == 0
            reportProgress(cancelCallback, progressCallback, ...
                0.6 * min(1, ftell(fid) / totalBytes), ...
                sprintf('Indexing chunk timestep %d', n));
        end
    end

    blockPos = blockPos(1:n);
    dataPos = dataPos(1:n);
    timesteps = timesteps(1:n);
    physicalTimes = physicalTimes(1:n);
    numChunks = numChunks(1:n);

    dt = NaN;
    isUniform = false;
    if numel(timesteps) >= 2
        dts = diff(timesteps);
        dt = dts(1);
        if all(abs(dts - dt) <= max(1e-12, 1e-9*max(abs(dt),1)))
            isUniform = true;
        end
    end

    idx = struct();
    idx.formatVersion = chunkIndexFormatVersion();
    idx.filePath = filePath;
    idx.fileSize = d.bytes;
    idx.fileDatenum = d.datenum;
    idx.quickSig = curSig;
    idx.blockPos = blockPos;
    idx.dataPos = dataPos;
    idx.timesteps = timesteps;
    idx.physicalTimes = physicalTimes;
    idx.numChunks = numChunks;
    idx.headerLines = preamble.headerLines;
    idx.varNames = preamble.varNames;
    idx.validVarNames = preamble.validVarNames;
    idx.metadata = preamble.metadata;
    idx.dt = dt;
    idx.isUniform = isUniform;

    try
        writeIndexCache(cachePath, idx);
    catch cacheError
        warning('read_chunk_step_fast:IndexCacheWriteFailed', ...
            'Index cache could not be updated: %s', cacheError.message);
    end
    tracker.finish();
end

function value = chunkIndexFormatVersion()
    value = 3;
end

function writeIndexCache(cachePath, idx)
    folder = fileparts(cachePath);
    temporaryPath = [tempname(folder), '.mat'];
    cleanupObj = onCleanup(@() deleteIfExists(temporaryPath)); %#ok<NASGU>
    save(temporaryPath, 'idx');
    [ok, message] = movefile(temporaryPath, cachePath, 'f');
    if ~ok
        error('read_chunk_step_fast:IndexCacheMoveFailed', '%s', message);
    end
end

function deleteIfExists(filePath)
    if exist(filePath, 'file')
        delete(filePath);
    end
end

function reportProgress(cancelCallback, progressCallback, fraction, message)
    drawnow;
    cancelled = cancelCallback();
    if ~(islogical(cancelled) && isscalar(cancelled))
        error('read_chunk_step_fast:BadCancelCallbackResult', ...
            'CancelCallback must return a logical scalar.');
    end
    if cancelled
        error('postdata:UserCancelled', 'Analysis cancelled by user.');
    end
    progressCallback(max(0, min(1, fraction)), message);
end

function validateReaderOptions(options)
    mode = strtrim(options.ProgressMode);
    if ~any(strcmpi(mode, {'auto','console','waitbar','off'}))
        error('read_chunk_step_fast:BadProgressMode', ...
            'ProgressMode must be auto, console, waitbar, or off.');
    end
    validateOptionalScalar(options.Index, 'Index');
    validateOptionalScalar(options.TimeStep, 'TimeStep');
    validateOptionalScalar(options.Time, 'Time');
    if ~(isnumeric(options.SlurmModuleIndex) && ...
            isscalar(options.SlurmModuleIndex) && ...
            isfinite(options.SlurmModuleIndex) && ...
            options.SlurmModuleIndex >= 1 && ...
            abs(options.SlurmModuleIndex - round(options.SlurmModuleIndex)) <= 1e-12)
        error('read_chunk_step_fast:BadSlurmModuleIndex', ...
            'SlurmModuleIndex must be a positive integer.');
    end
    if ~isempty(options.Index) && ...
            (options.Index < 1 || abs(options.Index - round(options.Index)) > 1e-12)
        error('read_chunk_step_fast:BadIndex', ...
            'Index must be a positive integer.');
    end
end

function validateOptionalScalar(value, name)
    if ~isempty(value) && ~(isnumeric(value) && isscalar(value) && isfinite(value))
        error('read_chunk_step_fast:BadSelectorValue', ...
            '%s must be empty or a finite numeric scalar.', name);
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

function tf = hasExplicitNameValue(args, name)
    tf = false;
    if isempty(args)
        return;
    end

    for i = 1:2:numel(args)
        key = args{i};
        if ~(ischar(key) || (isstring(key) && isscalar(key)))
            continue;
        end
        if strcmpi(char(key), name)
            tf = true;
            return;
        end
    end
end

function sig = buildQuickSignature(filePath)
% Build a lightweight signature from file size + head/tail bytes.

    info = dir(filePath);
    fileSize = info.bytes;

    sig = struct('fileSize', fileSize, ...
                 'headLen', 0, 'tailLen', 0, ...
                 'headSum', 0, 'headWeighted', 0, ...
                 'tailSum', 0, 'tailWeighted', 0);

    fid = fopen(filePath, 'rb');
    if fid < 0
        return;
    end
    cleanupObj = onCleanup(@() fclose(fid));

    sampleN = 1024;
    headN = min(sampleN, fileSize);
    if headN > 0
        head = fread(fid, headN, 'uint8=>double');
    else
        head = [];
    end

    tailN = min(sampleN, fileSize);
    if tailN > 0
        fseek(fid, fileSize - tailN, 'bof');
        tail = fread(fid, tailN, 'uint8=>double');
    else
        tail = [];
    end

    sig.headLen = numel(head);
    sig.tailLen = numel(tail);
    sig.headSum = sum(head);
    sig.tailSum = sum(tail);
    if ~isempty(head)
        sig.headWeighted = sum((1:numel(head))' .* head);
    end
    if ~isempty(tail)
        sig.tailWeighted = sum((1:numel(tail))' .* tail);
    end
end

function tf = isSameQuickSig(a, b)
    req = {'fileSize','headLen','tailLen','headSum','headWeighted','tailSum','tailWeighted'};
    for i = 1:numel(req)
        k = req{i};
        if ~isfield(a, k) || ~isfield(b, k)
            tf = false;
            return;
        end
        if ~isequal(a.(k), b.(k))
            tf = false;
            return;
        end
    end
    tf = true;
end
