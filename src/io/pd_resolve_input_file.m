function filePath = pd_resolve_input_file(request)
%PD_RESOLVE_INPUT_FILE Resolve an explicit or automatically discovered input.

    filePath = pd_to_char(request.filePath);
    if isempty(filePath)
        patterns = defaultPatterns(request);
        files = [];
        for i = 1:numel(patterns)
            matches = dir(fullfile(request.baseDir, patterns{i}));
            files = [files; matches(~[matches.isdir])]; %#ok<AGROW>
        end
        if numel(files) > 1
            [~, uniqueIndex] = unique({files.name}, 'stable');
            files = files(uniqueIndex);
        end
        files = files(~[files.isdir]);
        if isempty(files)
            error('postdata:NoInputFile', ...
                'No file matches %s in %s.', strjoin(patterns, ' or '), request.baseDir);
        end
        if numel(files) > 1
            error('postdata:MultipleInputFiles', ...
                'Multiple files match %s; select one explicitly.', strjoin(patterns, ' or '));
        end
        filePath = fullfile(request.baseDir, files(1).name);
    elseif ~isAbsolutePath(filePath)
        filePath = fullfile(request.baseDir, filePath);
    end

    if ~exist(filePath, 'file')
        error('postdata:InputFileNotFound', 'Input file not found: %s', filePath);
    end
end

function patterns = defaultPatterns(request)
    switch request.analysisType
        case 'chunk'
            if strcmpi(request.chunk.dimension, '1d')
                patterns = {'bin1d*.txt'};
            elseif strcmpi(request.chunk.dimension, '2d')
                patterns = {'bin2d*.txt'};
            else
                patterns = {'bin1d*.txt','bin2d*.txt'};
            end
        case 'network2d'
            patterns = {'bin2d*.txt'};
        case 'cluster'
            patterns = {'cluster_chunk*.txt'};
        case 'vx'
            patterns = {'mass_v*.txt','massv*.txt','vx_chunk*.txt'};
        case 'massx'
            patterns = {'mass_x*.txt','massx*.txt','bin1d*.txt'};
    end
end

function tf = isAbsolutePath(pathValue)
    tf = false;
    if numel(pathValue) >= 2 && pathValue(2) == ':'
        tf = true;
    elseif numel(pathValue) >= 2 && strcmp(pathValue(1:2), '\\')
        tf = true;
    end
end
