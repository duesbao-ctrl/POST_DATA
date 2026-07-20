function result = pd_add_result_metadata(result, analysisType)
%ADDRESULTMETADATA Add stable software/schema metadata to an analysis result.
% Existing result fields are preserved for backward compatibility.

    if ~isstruct(result)
        error('postdata:addResultMetadata:BadResult', ...
            'Analysis result must be a structure.');
    end

    info = pd_version();
    result.schemaVersion = info.resultSchemaVersion;
    result.analysisType = lower(strtrim(pd_to_char(analysisType)));
    result.software = struct('name', info.name, ...
                             'version', info.version, ...
                             'minimumMatlabRelease', info.minimumMatlabRelease);
end
