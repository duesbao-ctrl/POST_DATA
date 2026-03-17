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

    p = inputParser;
    p.addRequired('filePath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);      % Index | TimeStep | Time
    p.addParameter('Index', [], @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);
    p.parse(filePath, varargin{:});
    opt = p.Results;

    filePath = toChar(filePath);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);
    if ~exist(filePath, 'file')
        error('read_chunk_step_fast:FileNotFound', 'File not found: %s', filePath);
    end

    idx = buildOrLoadIndex(filePath, opt.ProgressMode);
    nBlocks = numel(idx.timesteps);
    if nBlocks == 0
        error('read_chunk_step_fast:NoBlocks', 'No timestep blocks found in file: %s', filePath);
    end

    mode = lower(strtrim(opt.SelectBy));
    % Backward compatibility: infer mode when SelectBy keeps default but
    % caller passes only Time or TimeStep.
    if strcmp(mode, 'index')
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
        if isempty(opt.SlurmPath)
            error('read_chunk_step_fast:MissingSlurmPath', ...
                'SlurmPath is required when selecting by Time.');
        end
        S = read_slurm_stepcpu(opt.SlurmPath);
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

    header1 = fgetl(fid);
    header2 = fgetl(fid);
    header3 = fgetl(fid);

    varNames = parseVarNames(header3);
    validVarNames = matlab.lang.makeValidName(varNames);
    numVars = numel(varNames);

    fseek(fid, idx.blockPos(sel), 'bof');
    bh = fgetl(fid);
    h = sscanf(bh, '%f').';
    if numel(h) < 3
        error('read_chunk_step_fast:BadBlockHeader', ...
            'Invalid block header near selected block: "%s"', bh);
    end
    timestep = h(1);
    numChunks = round(h(2));
    totalCount = h(3);

    fmt = repmat('%f', 1, numVars);
    C = textscan(fid, fmt, numChunks, 'CollectOutput', true, ...
        'Delimiter', ' \t', 'MultipleDelimsAsOne', true);
    data = C{1};

    if size(data, 1) ~= numChunks
        % Fallback robust path if textscan is short-read on malformed lines.
        data = zeros(numChunks, numVars);
        fseek(fid, idx.dataPos(sel), 'bof');
        for i = 1:numChunks
            ln = fgetl(fid);
            if ~ischar(ln)
                error('read_chunk_step_fast:UnexpectedEOF', ...
                    'Unexpected EOF while reading block at timestep %g', timestep);
            end
            v = sscanf(ln, '%f').';
            if numel(v) < numVars
                error('read_chunk_step_fast:BadDataRow', ...
                    'Data row has %d columns, expected %d', numel(v), numVars);
            end
            data(i, :) = v(1:numVars);
        end
    end

    colIndex = struct();
    for i = 1:numel(validVarNames)
        colIndex.(validVarNames{i}) = i;
    end

    data = applyEmptyCellNaN(data, validVarNames);

    out = struct();
    out.filePath = filePath;
    out.headerLines = {header1, header2, header3};
    out.varNames = varNames;
    out.validVarNames = validVarNames;
    out.colIndex = colIndex;
    out.requestedTimeStep = reqStep;
    out.stepIndex = sel;
    out.timestep = timestep;
    out.numChunks = numChunks;
    out.totalCount = totalCount;
    out.data = data;
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

function idx = buildOrLoadIndex(filePath, progressMode)
    d = dir(filePath);
    cachePath = [filePath, '.stepidx.mat'];
    curSig = buildQuickSignature(filePath);

    if exist(cachePath, 'file')
        S = load(cachePath, 'idx');
        if isfield(S, 'idx') && isfield(S.idx, 'fileSize') && isfield(S.idx, 'fileDatenum')
            hasSig = isfield(S.idx, 'quickSig');
            sameSig = hasSig && isSameQuickSig(S.idx.quickSig, curSig);
            if isequal(S.idx.fileSize, d.bytes) && isequal(S.idx.fileDatenum, d.datenum) && sameSig
                idx = S.idx;
                return;
            end
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

    fgetl(fid); fgetl(fid); fgetl(fid);
    tracker.update(ftell(fid) / totalBytes);

    blockPos = zeros(1000, 1);
    dataPos  = zeros(1000, 1);
    timesteps = zeros(1000, 1);
    numChunks = zeros(1000, 1);
    n = 0;

    while true
        pos = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if isempty(strtrim(line))
            continue;
        end

        h = sscanf(line, '%f').';
        if numel(h) < 3
            continue;
        end

        n = n + 1;
        if n > numel(blockPos)
            blockPos = [blockPos; zeros(numel(blockPos),1)]; %#ok<AGROW>
            dataPos  = [dataPos;  zeros(numel(dataPos),1)]; %#ok<AGROW>
            timesteps = [timesteps; zeros(numel(timesteps),1)]; %#ok<AGROW>
            numChunks = [numChunks; zeros(numel(numChunks),1)]; %#ok<AGROW>
        end

        blockPos(n) = pos;
        dataPos(n) = ftell(fid);
        timesteps(n) = h(1);
        numChunks(n) = round(h(2));

        for i = 1:numChunks(n)
            if ~ischar(fgetl(fid))
                break;
            end
        end
        tracker.update(ftell(fid) / totalBytes);
    end

    blockPos = blockPos(1:n);
    dataPos = dataPos(1:n);
    timesteps = timesteps(1:n);
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
    idx.filePath = filePath;
    idx.fileSize = d.bytes;
    idx.fileDatenum = d.datenum;
    idx.quickSig = curSig;
    idx.blockPos = blockPos;
    idx.dataPos = dataPos;
    idx.timesteps = timesteps;
    idx.numChunks = numChunks;
    idx.dt = dt;
    idx.isUniform = isUniform;

    try
        save(cachePath, 'idx');
    catch
        % ignore cache write failure
    end
    tracker.finish();
end

function varNames = parseVarNames(headerLine)
    txt = strtrim(headerLine);
    if ~isempty(txt) && txt(1) == '#'
        txt = strtrim(txt(2:end));
    end
    varNames = strsplit(txt);
    if isempty(varNames)
        error('read_chunk_step_fast:NoVarNames', 'No variable names found in third header line.');
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
