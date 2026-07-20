function frame = pd_read_next_chunk_frame_header(fid, callerId)
%PD_READ_NEXT_CHUNK_FRAME_HEADER Find the next numeric chunk frame header.
% Blank lines and comment lines are accepted between frames. Current SPID
% files provide "# Time <physical-time>" immediately before each header.

    if nargin < 2 || isempty(callerId)
        callerId = 'postdata:chunkFrame';
    end
    callerId = pd_to_char(callerId);

    frame = emptyFrame();
    physicalTime = NaN;
    prefixLines = {};

    while true
        linePosition = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line)
            frame.eof = true;
            frame.prefixLines = prefixLines;
            return;
        end

        text = strtrim(line);
        if isempty(text)
            prefixLines{end + 1} = line; %#ok<AGROW>
            continue;
        end
        if text(1) == '#'
            prefixLines{end + 1} = line; %#ok<AGROW>
            content = strtrim(text(2:end));
            timeToken = regexp(content, ...
                '^Time\s+([^\s]+)\s*$', 'tokens', 'once');
            if ~isempty(timeToken)
                value = str2double(timeToken{1});
                if ~(isscalar(value) && isfinite(value))
                    error([callerId, ':BadPhysicalTime'], ...
                        'Invalid SPID physical time comment: "%s"', line);
                end
                physicalTime = value;
            elseif ~isempty(regexp(content, '^Time(?:\s|$)', 'once'))
                error([callerId, ':BadPhysicalTime'], ...
                    'Invalid SPID physical time comment: "%s"', line);
            end
            continue;
        end

        values = pd_parse_chunk_summary_line(line, callerId);
        frame.eof = false;
        frame.blockPosition = linePosition;
        frame.dataPosition = ftell(fid);
        frame.summaryLine = line;
        frame.timestep = values(1);
        frame.numChunks = round(values(2));
        frame.totalCount = values(3);
        frame.physicalTime = physicalTime;
        frame.prefixLines = prefixLines;
        return;
    end
end

function frame = emptyFrame()
    frame = struct();
    frame.eof = false;
    frame.blockPosition = NaN;
    frame.dataPosition = NaN;
    frame.summaryLine = '';
    frame.timestep = NaN;
    frame.numChunks = 0;
    frame.totalCount = NaN;
    frame.physicalTime = NaN;
    frame.prefixLines = {};
end
