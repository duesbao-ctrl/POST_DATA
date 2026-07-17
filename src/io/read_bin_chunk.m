function out = read_bin_chunk(filePath, varargin)
%READ_BIN_CHUNK Read chunk-averaged data blocks from LAMMPS-style text output.
%   out = READ_BIN_CHUNK(filePath)
%   out = READ_BIN_CHUNK(filePath, 'ProgressMode', 'console')
%
% File format:
%   extensible comment preamble ending in "# Chunk ..."
%   then repeated blocks, optionally prefixed by "# Time ...":
%     timestep number_of_chunks total_count
%     <number_of_chunks> lines of numeric data
%
% Output fields:
%   out.filePath
%   out.headerLines
%   out.varNames
%   out.validVarNames
%   out.colIndex
%   out.steps(i).timestep
%   out.steps(i).numChunks
%   out.steps(i).totalCount
%   out.steps(i).data

    p = inputParser;
    p.addRequired('filePath', @isTextScalar);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);
    p.parse(filePath, varargin{:});

    filePath = toChar(p.Results.filePath);
    progressMode = toChar(p.Results.ProgressMode);
    info = dir(filePath);
    if isempty(info)
        fileSize = 1;
    else
        fileSize = max(info(1).bytes, 1);
    end

    fid = fopen(filePath, 'r');
    if fid < 0
        error('read_bin_chunk:FileOpenFailed', 'Cannot open file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); 
    tracker = make_file_progress(progressMode, 'Reading file', filePath);
    progressCleanup = onCleanup(@() tracker.close()); %#ok<NASGU>

    preamble = pd_read_chunk_preamble(fid, 'read_bin_chunk');
    tracker.update(ftell(fid) / fileSize);

    varNames = preamble.varNames;
    validVarNames = preamble.validVarNames;
    numVars = numel(varNames);

    stepsCell = cell(128, 1);
    stepCount = 0;

    while true
        frame = pd_read_next_chunk_frame_header(fid, 'read_bin_chunk');
        if frame.eof
            break;
        end

        timestep = frame.timestep;
        numChunks = frame.numChunks;
        totalCount = frame.totalCount;
        data = zeros(numChunks, numVars);

        row = 1;
        while row <= numChunks
            dataLine = fgetl(fid);
            if ~ischar(dataLine)
                error('read_bin_chunk:UnexpectedEOF', ...
                    'Unexpected EOF while reading timestep %g in %s', timestep, filePath);
            end

            if isempty(strtrim(dataLine))
                continue;
            end

            values = sscanf(dataLine, '%f').';
            if numel(values) < numVars
                error('read_bin_chunk:BadDataRow', ...
                    'Data row has %d columns but expected %d. Line: "%s"', ...
                    numel(values), numVars, dataLine);
            end

            data(row, :) = values(1:numVars);
            row = row + 1;
        end

        data = applyEmptyCellNaN(data, validVarNames);

        stepCount = stepCount + 1;
        if stepCount > numel(stepsCell)
            stepsCell = [stepsCell; cell(numel(stepsCell), 1)]; %#ok<AGROW>
        end
        stepsCell{stepCount, 1} = struct( ...
            'timestep', timestep, ...
            'physicalTime', frame.physicalTime, ...
            'numChunks', numChunks, ...
            'totalCount', totalCount, ...
            'data', data);
        tracker.update(ftell(fid) / fileSize);
    end

    if stepCount == 0
        steps = struct('timestep', {}, 'physicalTime', {}, ...
            'numChunks', {}, 'totalCount', {}, 'data', {});
    else
        steps = [stepsCell{1:stepCount}];
    end

    colIndex = struct();
    for k = 1:numVars
        colIndex.(validVarNames{k}) = k;
    end

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
    out.steps = steps;
    tracker.finish();
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
