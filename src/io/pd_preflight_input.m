function report = pd_preflight_input(request, filePath)
%PD_PREFLIGHT_INPUT Validate input structure before expensive calculations.

    before = dir(filePath);
    if isempty(before) || before.isdir
        error('postdata:PreflightMissingFile', 'Input file does not exist: %s', filePath);
    end
    fid = fopen(filePath, 'r');
    if fid < 0
        error('postdata:PreflightOpenFailed', 'Cannot open input file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    header = cell(1, 3);
    for i = 1:3
        header{i} = fgetl(fid);
        if ~ischar(header{i})
            error('postdata:PreflightShortHeader', ...
                'Input must contain at least three header lines.');
        end
    end
    blockLine = '';
    while ischar(blockLine)
        blockLine = fgetl(fid);
        if ~ischar(blockLine), break; end
        if ~isempty(strtrim(blockLine)), break; end
    end
    blockValues = sscanf(blockLine, '%f').';
    if numel(blockValues) < 3 || ~all(isfinite(blockValues(1:3))) || ...
            blockValues(2) < 0 || ...
            abs(blockValues(2) - round(blockValues(2))) > 1e-12
        error('postdata:PreflightBadBlockHeader', ...
            ['First block header must contain a finite timestep, a non-negative ' ...
             'integer row count, and a finite total count.']);
    end

    [variables, validNames] = pd_parse_chunk_variables( ...
        header{3}, 'postdata:preflight');
    firstDataColumns = 0;
    if blockValues(2) > 0
        firstDataLine = fgetl(fid);
        while ischar(firstDataLine) && isempty(strtrim(firstDataLine))
            firstDataLine = fgetl(fid);
        end
        if ~ischar(firstDataLine)
            error('postdata:PreflightMissingDataRow', ...
                'The first block declares rows but contains no data.');
        end
        firstData = sscanf(firstDataLine, '%f').';
        firstDataColumns = numel(firstData);
        if firstDataColumns ~= numel(variables)
            error('postdata:PreflightDataColumnMismatch', ...
                'First data row has %d columns; the header declares %d.', ...
                firstDataColumns, numel(variables));
        end
    end
    required = requiredVariables(request);
    for i = 1:numel(required)
        alternatives = required{i};
        if ~any(ismember(alternatives, validNames))
            error('postdata:PreflightMissingVariable', ...
                'Input is missing required variable. Expected one of: %s.', ...
                strjoin(alternatives, ', '));
        end
    end

    after = dir(filePath);
    warnings = {};
    if before.bytes ~= after.bytes || before.datenum ~= after.datenum
        warnings{end + 1} = 'Input file changed during preflight and may still be written.'; %#ok<AGROW>
    end
    report = struct();
    report.filePath = filePath;
    report.fileSize = after.bytes;
    report.fileDatenum = after.datenum;
    report.analysisType = request.analysisType;
    report.headerLines = header;
    report.variables = variables;
    report.validVariableNames = validNames;
    report.firstTimestep = blockValues(1);
    report.firstBlockRows = round(blockValues(2));
    report.firstDataColumns = firstDataColumns;
    report.warnings = warnings;
    report.passed = true;
end

function required = requiredVariables(request)
    required = {};
    switch request.analysisType
        case 'chunk'
            if strcmpi(request.chunk.dimension, '2d')
                required = {{'Coord1','c_x','Chunk'}, {'Coord2','c_y'}};
            else
                required = {{'Coord1','c_x','Chunk'}};
            end
        case 'cluster'
            required = {{pd_get_analysis_option(request.analysisOptions, 'NcountVar', 'Ncount')}};
        case 'vx'
            velocityVar = pd_get_analysis_option(request.analysisOptions, ...
                'VelocityVar', 'Chunk');
            required = {{matlab.lang.makeValidName(pd_to_char(velocityVar))}};
        case 'massx'
            coordinateVar = pd_get_analysis_option(request.analysisOptions, ...
                'CoordinateVar', 'auto');
            coordinateVar = pd_to_char(coordinateVar);
            if strcmpi(strtrim(coordinateVar), 'auto')
                required = {{'Coord1','c_x','Chunk'}};
            else
                required = {{matlab.lang.makeValidName(coordinateVar)}};
            end
        case 'network2d'
            ncount = pd_get_analysis_option(request.analysisOptions, 'NcountVar', 'Ncount');
            required = {{'Coord1','c_x'}, {'Coord2','c_y'}, {ncount}};
    end
end
