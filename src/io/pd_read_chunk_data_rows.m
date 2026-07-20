function data = pd_read_chunk_data_rows(fid, numRows, numVars, callerId, timestep, rowCallback)
%PD_READ_CHUNK_DATA_ROWS Read exact-width, line-oriented chunk data rows.
% When called without an output, rows are validated without being stored.

    if nargin < 4 || isempty(callerId)
        callerId = 'postdata:chunkRows';
    end
    if nargin < 5, timestep = NaN; end
    if nargin < 6, rowCallback = []; end
    callerId = pd_to_char(callerId);
    storeData = nargout > 0;
    if storeData
        data = zeros(numRows, numVars);
    else
        data = [];
    end

    row = 1;
    while row <= numRows
        line = fgetl(fid);
        if ~ischar(line)
            error([callerId, ':UnexpectedEOF'], ...
                'Unexpected EOF while reading timestep %g.', timestep);
        end
        if isempty(strtrim(line))
            continue;
        end
        [values, ~, ~, nextIndex] = sscanf(line, '%f');
        values = values.';
        trailingText = strtrim(line(nextIndex:end));
        if numel(values) ~= numVars || ~isempty(trailingText)
            error([callerId, ':BadDataRow'], ...
                ['Data row does not contain exactly %d numeric columns ' ...
                 '(parsed %d). Timestep %g, line: "%s"'], ...
                numVars, numel(values), timestep, line);
        end
        if storeData
            data(row, :) = values;
        end
        if ~isempty(rowCallback) && ...
                (mod(row - 1, 256) == 0 || row == numRows)
            rowCallback(row, numRows);
        end
        row = row + 1;
    end
end
