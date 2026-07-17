function request = pd_create_request(analysisType)
%CREATEREQUEST Create a versioned POST_DATA analysis request structure.
% MATLAB R2016b compatible.

    if nargin < 1 || ~pd_is_text_scalar(analysisType)
        error('postdata:createRequest:BadType', 'analysisType must be text.');
    end

    analysisType = pd_normalize_analysis_type(analysisType);
    validTypes = {'chunk', 'cluster', 'vx', 'massx', 'network2d'};
    if ~any(strcmp(analysisType, validTypes))
        error('postdata:createRequest:BadType', ...
            'analysisType must be chunk/cluster/vx/massx/network2d.');
    end

    request = struct();
    request.schemaVersion = 1;
    request.analysisType = analysisType;
    request.baseDir = pwd;
    request.filePath = '';
    request.selection = struct('mode', 'Index', 'value', 1, ...
        'slurmPath', '', 'slurmModuleIndex', 1);
    request.progressMode = 'auto';
    request.makePlots = true;
    request.plotOptions = struct();
    request.resultLevel = 'standard';
    request.execution = struct('cancelCallback', @() false, ...
        'progressCallback', @(fraction, message) [], ...
        'logFile', '');
    request.chunk = struct('dimension', 'auto', 'variable', 'c_rho', ...
        'dV', [], 'coordScale', 1, 'coordRangeX', [], 'coordRangeY', [], ...
        'coordRangeZ', [], ...
        'gradientVariable', '', 'gradientSmoothLevel', 0, ...
        'strainRateVelocityComponent', 'vz', ...
        'strainRateDensityVariable', 'c_rho', ...
        'plotOptions', {{}});
    request.analysisOptions = {};
end
