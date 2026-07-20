function pd_write_manifest(filePath, result, options)
%PD_WRITE_MANIFEST Save human-readable reproducibility metadata.

    fid = fopen(filePath, 'w');
    if fid < 0, error('postdata:ManifestOpenFailed', 'Cannot create manifest: %s', filePath); end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'POST_DATA2 reproducibility manifest\n');
    fprintf(fid, 'analysisType=%s\n', result.analysisType);
    fprintf(fid, 'softwareVersion=%s\n', result.software.version);
    fprintf(fid, 'minimumMatlabRelease=%s\n', result.software.minimumMatlabRelease);
    fprintf(fid, 'matlabVersion=%s\n', version);
    if isfield(result, 'run')
        fprintf(fid, 'startedAt=%s\n', result.run.startedAt);
        fprintf(fid, 'durationSeconds=%.16g\n', result.run.durationSeconds);
        fprintf(fid, 'resultLevel=%s\n', result.storage.level);
    end
    if isfield(result, 'preflight')
        fprintf(fid, 'inputFile=%s\n', result.preflight.filePath);
        fprintf(fid, 'inputBytes=%d\n', result.preflight.fileSize);
        fprintf(fid, 'inputDatenum=%.16g\n', result.preflight.fileDatenum);
        fprintf(fid, 'inputSignature=bytes:%d;datenum:%.16g\n', ...
            result.preflight.fileSize, result.preflight.fileDatenum);
        fprintf(fid, 'variables=%s\n', strjoin(result.preflight.variables, ','));
    end
    fprintf(fid, '\n[Request]\n');
    if isfield(result, 'run')
        writeScalarStruct(fid, result.run.request, 'request');
        analysisOptions = result.run.request.analysisOptions;
        for i = 1:2:numel(analysisOptions)
            fprintf(fid, 'request.analysisOptions.%s=%s\n', ...
                pd_to_char(analysisOptions{i}), pd_format_option_value(analysisOptions{i + 1}));
        end
    end
    fprintf(fid, '\n[Output]\n');
    names = fieldnames(options);
    for i = 1:numel(names)
        fprintf(fid, 'output.%s=%s\n', names{i}, pd_format_option_value(options.(names{i})));
    end
end

function writeScalarStruct(fid, value, prefix)
    if ~isstruct(value) || numel(value) ~= 1, return; end
    names = fieldnames(value);
    for i = 1:numel(names)
        item = value.(names{i});
        path = [prefix, '.', names{i}];
        if isstruct(item)
            writeScalarStruct(fid, item, path);
        elseif (isnumeric(item) || islogical(item)) && (isscalar(item) || isempty(item))
            fprintf(fid, '%s=%s\n', path, pd_format_option_value(item));
        elseif ischar(item)
            fprintf(fid, '%s=%s\n', path, item);
        end
    end
end
