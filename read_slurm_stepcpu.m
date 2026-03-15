function out = read_slurm_stepcpu(filePath)
%READ_SLURM_STEPCPU Read thermo table blocks starting with "Step CPU".
%   out = read_slurm_stepcpu(filePath)
%
% Output:
%   out.filePath
%   out.modules(i).headerLine
%   out.modules(i).varNames
%   out.modules(i).validVarNames
%   out.modules(i).colIndex
%   out.modules(i).data
%   out.modules(i).numRows
%   out.modules(i).acceptedRows
%   out.modules(i).skippedRows

    if nargin < 1 || isempty(filePath)
        error('read_slurm_stepcpu:MissingInput', 'filePath is required.');
    end
    if isstring(filePath)
        filePath = char(filePath);
    end
    if ~ischar(filePath)
        error('read_slurm_stepcpu:InvalidInputType', 'filePath must be char or string.');
    end

    fid = fopen(filePath, 'r');
    if fid < 0
        error('read_slurm_stepcpu:FileOpenFailed', 'Cannot open file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    modulesCell = {};
    moduleCount = 0;

    hasPending = false;
    pendingLine = '';

    while true
        if hasPending
            line = pendingLine;
            hasPending = false;
        else
            line = fgetl(fid);
        end

        if ~ischar(line)
            break;
        end

        if ~isStepCpuHeader(line)
            continue;
        end

        headerNorm = normalizeHeaderLine(line);
        varNames = strsplit(headerNorm);
        nCols = numel(varNames);
        validVarNames = matlab.lang.makeValidName(varNames);

        colIndex = struct();
        for k = 1:nCols
            colIndex.(validVarNames{k}) = k;
        end

        data = zeros(256, nCols);
        rowCount = 0;
        skippedRows = 0;

        while true
            nextLine = fgetl(fid);
            if ~ischar(nextLine)
                break;
            end

            if isStepCpuHeader(nextLine)
                pendingLine = nextLine;
                hasPending = true;
                break;
            end

            [ok, values] = parseNumericRow(nextLine, nCols);
            if ok
                rowCount = rowCount + 1;
                if rowCount > size(data, 1)
                    data = [data; zeros(size(data, 1), nCols)]; %#ok<AGROW>
                end
                data(rowCount, :) = values;
            else
                skippedRows = skippedRows + 1;
            end
        end

        if rowCount == 0
            data = zeros(0, nCols);
        else
            data = data(1:rowCount, :);
        end

        moduleCount = moduleCount + 1;
        modulesCell{moduleCount, 1} = struct( ... %#ok<AGROW>
            'headerLine', line, ...
            'varNames', {varNames}, ...
            'validVarNames', {validVarNames}, ...
            'colIndex', colIndex, ...
            'data', data, ...
            'numRows', rowCount, ...
            'acceptedRows', rowCount, ...
            'skippedRows', skippedRows);
    end

    if moduleCount == 0
        error('read_slurm_stepcpu:NoStepCpuBlock', ...
            ['No "Step CPU" data block found in file: %s\n' ...
             'Check file path/encoding or log content.'], filePath);
    end
    modules = [modulesCell{:}];

    out = struct();
    out.filePath = filePath;
    out.modules = modules;
end

function tf = isStepCpuHeader(line)
%ISSTEPCPUHEADER Check whether a line is a "Step CPU" header line.

    txt = normalizeHeaderLine(line);
    if isempty(txt)
        tf = false;
        return;
    end
    parts = strsplit(txt);
    tf = (numel(parts) >= 2) && strcmpi(parts{1}, 'Step') && strcmpi(parts{2}, 'CPU');
end

function txt = normalizeHeaderLine(line)
% Remove leading markers and normalize spaces for header parsing.

    txt = strtrim(line);
    while ~isempty(txt) && (txt(1) == '#' || txt(1) == '%' || txt(1) == '*')
        txt = strtrim(txt(2:end));
    end
end

function [ok, values] = parseNumericRow(line, nCols)
%PARSENUMERICROW Parse one row; valid only if all columns are numeric.

    values = [];
    txt = strtrim(line);
    if isempty(txt)
        ok = false;
        return;
    end

    parts = strsplit(txt);
    if numel(parts) ~= nCols
        ok = false;
        return;
    end

    nums = str2double(parts);
    if any(isnan(nums)) || any(~isfinite(nums))
        ok = false;
        return;
    end

    values = nums;
    ok = true;
end
