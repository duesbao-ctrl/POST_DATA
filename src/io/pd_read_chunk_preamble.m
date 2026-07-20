function preamble = pd_read_chunk_preamble(fid, callerId)
%PD_READ_CHUNK_PREAMBLE Read an extensible chunk-file preamble.
% Supports legacy three-line LAMMPS-style headers and current SPID headers
% containing Units, Name/Kind, and other comment metadata.

    if nargin < 2 || isempty(callerId)
        callerId = 'postdata:chunkPreamble';
    end
    callerId = pd_to_char(callerId);

    headerLines = {};
    metadata = emptyMetadata();
    columnHeaderLine = '';
    columnHeaderPosition = NaN;

    while true
        linePosition = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line)
            error([callerId, ':MissingColumnHeader'], ...
                'No "# Chunk ..." column header was found before end of file.');
        end
        headerLines{end + 1} = line; %#ok<AGROW>
        text = strtrim(line);
        if isempty(text)
            continue;
        end
        if text(1) ~= '#'
            error([callerId, ':UnexpectedPreambleLine'], ...
                ['Expected a comment preamble ending with a "# Chunk ..." ' ...
                 'column header. Found: "%s"'], line);
        end

        content = strtrim(text(2:end));
        tokens = strsplit(content);
        if ~isempty(tokens) && strcmp(tokens{1}, 'Chunk')
            columnHeaderLine = line;
            columnHeaderPosition = linePosition;
            break;
        end
        metadata = parseMetadataLine(metadata, content);
    end

    [varNames, validVarNames] = pd_parse_chunk_variables( ...
        columnHeaderLine, callerId);
    metadata.columnHeaderLine = columnHeaderLine;
    metadata.preambleLineCount = numel(headerLines);
    if metadata.isSpid
        metadata.inputFormat = 'spid-chunk';
    else
        metadata.inputFormat = 'legacy-chunk';
    end

    preamble = struct();
    preamble.headerLines = headerLines;
    preamble.columnHeaderLine = columnHeaderLine;
    preamble.columnHeaderPosition = columnHeaderPosition;
    preamble.dataStartPosition = ftell(fid);
    preamble.varNames = varNames;
    preamble.validVarNames = validVarNames;
    preamble.metadata = metadata;
end

function metadata = parseMetadataLine(metadata, content)
    if strcmpi(content, 'Chunk-averaged data for SPID')
        metadata.marker = content;
        metadata.isSpid = true;
        return;
    end

    value = regexp(content, '^Units\s+(\S+)\s*$', 'tokens', 'once');
    if ~isempty(value)
        metadata.unitSystem = lower(strtrim(value{1}));
        metadata.isSpid = true;
        return;
    end

    value = regexp(content, ...
        '^Name\s+(.+?)\s+Kind\s+(\S+)\s*$', 'tokens', 'once');
    if ~isempty(value)
        metadata.name = strtrim(value{1});
        metadata.kind = lower(strtrim(value{2}));
        metadata.isSpid = true;
    end
end

function metadata = emptyMetadata()
    metadata = struct();
    metadata.inputFormat = 'legacy-chunk';
    metadata.marker = '';
    metadata.unitSystem = '';
    metadata.name = '';
    metadata.kind = '';
    metadata.isSpid = false;
    metadata.columnHeaderLine = '';
    metadata.preambleLineCount = 0;
end
