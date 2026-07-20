function [varNames, validVarNames] = pd_parse_chunk_variables(headerLine, callerId)
%PD_PARSE_CHUNK_VARIABLES Parse and validate a chunk variable header.
% Rejects names that collapse to the same MATLAB structure field because
% silently overwriting one column would corrupt every downstream analysis.

    if nargin < 2 || isempty(callerId)
        callerId = 'postdata:chunkHeader';
    end
    callerId = pd_to_char(callerId);
    if isstring(headerLine) && isscalar(headerLine)
        headerLine = char(headerLine);
    end
    if ~ischar(headerLine)
        error([callerId, ':MissingHeader'], ...
            'The chunk variable header is missing or is not text.');
    end

    text = strtrim(headerLine);
    if ~isempty(text) && text(1) == '#'
        text = strtrim(text(2:end));
    end
    varNames = strsplit(text);
    if isempty(varNames) || (isscalar(varNames) && isempty(varNames{1}))
        error([callerId, ':NoVariableNames'], ...
            'No variable names were found in the chunk column header.');
    end

    validVarNames = matlab.lang.makeValidName(varNames);
    duplicateIndex = firstDuplicate(validVarNames);
    if duplicateIndex > 0
        previous = find(strcmp(validVarNames{duplicateIndex}, ...
            validVarNames(1:duplicateIndex - 1)), 1, 'first');
        error([callerId, ':AmbiguousVariableName'], ...
            ['Header variables "%s" and "%s" both normalize to "%s". ' ...
             'Rename one input column before analysis.'], ...
            varNames{previous}, varNames{duplicateIndex}, ...
            validVarNames{duplicateIndex});
    end
end

function index = firstDuplicate(values)
    index = 0;
    for i = 2:numel(values)
        if any(strcmp(values{i}, values(1:i - 1)))
            index = i;
            return;
        end
    end
end
