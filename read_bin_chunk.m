function out = read_bin_chunk(filePath)
%READ_BIN_CHUNK Read chunk-averaged data blocks from LAMMPS-style text output.
%   out = READ_BIN_CHUNK(filePath)
%
% File format:
%   line 1-2: comments
%   line 3  : variable names (usually prefixed by '#')
%   then repeated blocks:
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

    if nargin < 1 || isempty(filePath)
        error('read_bin_chunk:MissingInput', 'filePath is required.');
    end
    if isstring(filePath)
        filePath = char(filePath);
    end
    if ~ischar(filePath)
        error('read_bin_chunk:InvalidInputType', 'filePath must be char or string.');
    end

    fid = fopen(filePath, 'r');
    if fid < 0
        error('read_bin_chunk:FileOpenFailed', 'Cannot open file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); 

    header1 = fgetl(fid);
    header2 = fgetl(fid);
    header3 = fgetl(fid);

    if ~ischar(header1) || ~ischar(header2) || ~ischar(header3)
        error('read_bin_chunk:HeaderTooShort', 'File header is incomplete: %s', filePath);
    end

    varNames = parseVarNames(header3);
    numVars = numel(varNames);
    validVarNames = matlab.lang.makeValidName(varNames);

    stepsCell = {};
    stepCount = 0;

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        if isempty(strtrim(line))
            continue;
        end

        blockHeader = sscanf(line, '%f').';
        if numel(blockHeader) < 3
            error('read_bin_chunk:BadBlockHeader', ...
                'Invalid block header in file %s: "%s"', filePath, line);
        end

        timestep = blockHeader(1);
        numChunks = blockHeader(2);
        totalCount = blockHeader(3);

        if numChunks < 0 || abs(numChunks - round(numChunks)) > 0
            error('read_bin_chunk:BadChunkCount', ...
                'Number-of-chunks must be a non-negative integer. Got: %g', numChunks);
        end

        numChunks = round(numChunks);
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
        stepsCell{stepCount, 1} = struct( ... %#ok<AGROW>
            'timestep', timestep, ...
            'numChunks', numChunks, ...
            'totalCount', totalCount, ...
            'data', data);
    end

    if stepCount == 0
        steps = struct('timestep', {}, 'numChunks', {}, 'totalCount', {}, 'data', {});
    else
        steps = [stepsCell{:}];
    end

    colIndex = struct();
    for k = 1:numVars
        colIndex.(validVarNames{k}) = k;
    end

    out = struct();
    out.filePath = filePath;
    out.headerLines = {header1, header2, header3};
    out.varNames = varNames;
    out.validVarNames = validVarNames;
    out.colIndex = colIndex;
    out.steps = steps;
end

function varNames = parseVarNames(headerLine)
%PARSEVARNAMES Extract variable names from the third header line.

    txt = strtrim(headerLine);
    if ~isempty(txt) && txt(1) == '#'
        txt = strtrim(txt(2:end));
    end

    parts = strsplit(txt);
    if isempty(parts)
        error('read_bin_chunk:NoVarNames', 'No variable names found in third header line.');
    end

    varNames = parts;
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
