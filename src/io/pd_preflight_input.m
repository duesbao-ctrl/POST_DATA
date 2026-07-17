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
    preamble = pd_read_chunk_preamble(fid, 'postdata:preflight');
    frame = pd_read_next_chunk_frame_header(fid, 'postdata:preflight');
    if frame.eof
        error('postdata:PreflightMissingBlock', ...
            'Input contains a chunk header but no timestep blocks.');
    end
    variables = preamble.varNames;
    validNames = preamble.validVarNames;
    firstDataColumns = 0;
    if frame.numChunks > 0
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
    report.headerLines = preamble.headerLines;
    report.metadata = preamble.metadata;
    report.inputFormat = preamble.metadata.inputFormat;
    report.unitSystem = preamble.metadata.unitSystem;
    report.taskName = preamble.metadata.name;
    report.chunkKind = preamble.metadata.kind;
    report.variables = variables;
    report.validVariableNames = validNames;
    report.firstTimestep = frame.timestep;
    report.firstPhysicalTime = frame.physicalTime;
    report.firstBlockRows = frame.numChunks;
    report.firstDataColumns = firstDataColumns;
    report.warnings = warnings;
    report.passed = true;
end

function required = requiredVariables(request)
    required = {};
    switch request.analysisType
        case 'chunk'
            if strcmpi(request.chunk.dimension, '3d')
                required = {{'Coord1','c_x','Chunk'}, ...
                    {'Coord2','c_y'}, {'Coord3','c_z'}};
            elseif strcmpi(request.chunk.dimension, '2d')
                required = {{'Coord1','c_x','Chunk'}, {'Coord2','c_y'}};
            else
                required = {{'Coord1','c_x','Chunk'}};
            end
        case 'cluster'
            required = {{pd_get_analysis_option(request.analysisOptions, 'NcountVar', 'Ncount')}};
        case 'vx'
            velocityVar = pd_get_analysis_option(request.analysisOptions, ...
                'VelocityVar', 'auto');
            velocityVar = pd_to_char(velocityVar);
            if strcmpi(strtrim(velocityVar), 'auto')
                required = {{'Coord1','Chunk'}};
            else
                required = {{matlab.lang.makeValidName(velocityVar)}};
            end
        case 'massx'
            coordinateVar = pd_get_analysis_option(request.analysisOptions, ...
                'CoordinateVar', 'auto');
            coordinateVar = pd_to_char(coordinateVar);
            if strcmpi(strtrim(coordinateVar), 'auto')
                required = {{'Coord1','c_x','Chunk'}};
            else
                required = {{matlab.lang.makeValidName(coordinateVar)}};
            end
            ncount = pd_get_analysis_option(request.analysisOptions, ...
                'NcountVar', 'Ncount');
            required{end + 1} = {matlab.lang.makeValidName(pd_to_char(ncount))};
            chunkDim = pd_get_analysis_option(request.analysisOptions, ...
                'ChunkDim', 'auto');
            if strcmpi(strtrim(pd_to_char(chunkDim)), '2d')
                coordinateVarY = pd_get_analysis_option( ...
                    request.analysisOptions, 'CoordinateVarY', 'auto');
                coordinateVarY = pd_to_char(coordinateVarY);
                if strcmpi(strtrim(coordinateVarY), 'auto')
                    required{end + 1} = {'Coord2','c_y'};
                else
                    required{end + 1} = { ...
                        matlab.lang.makeValidName(coordinateVarY)};
                end
            end
        case 'network2d'
            ncount = pd_get_analysis_option(request.analysisOptions, 'NcountVar', 'Ncount');
            required = {{'Coord1','c_x'}, {'Coord2','c_y'}, {ncount}};
    end
end
