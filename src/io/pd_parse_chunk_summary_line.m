function values = pd_parse_chunk_summary_line(line, callerId)
%PD_PARSE_CHUNK_SUMMARY_LINE Parse an exact three-value frame summary.

    if nargin < 2 || isempty(callerId)
        callerId = 'postdata:chunkFrame';
    end
    callerId = pd_to_char(callerId);
    if ~ischar(line)
        error([callerId, ':BadBlockHeader'], ...
            'Chunk frame summary must be a text line.');
    end
    [values, ~, ~, nextIndex] = sscanf(line, '%f');
    values = values.';
    trailingText = strtrim(line(nextIndex:end));
    if numel(values) ~= 3 || ~isempty(trailingText) || ...
            ~all(isfinite(values)) || values(2) < 0 || ...
            abs(values(2) - round(values(2))) > 1e-12
        error([callerId, ':BadBlockHeader'], ...
            ['Invalid block header. Expected exactly a finite timestep, ' ...
             'non-negative integer row count, and finite total count. ' ...
             'Line: "%s"'], line);
    end
end
