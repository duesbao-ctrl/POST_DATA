function result = postdata_run(request)
%POSTDATA_RUN Execute a versioned POST_DATA analysis request.
% This is the single public analysis API. MATLAB R2016b compatible.

    postdata_startup();
    request = pd_validate_request(request);
    request.plotOptions = pd_normalize_plot_options(request.plotOptions);
    logFile = request.execution.logFile;
    if isempty(logFile), logFile = pd_default_log_file(); end
    startedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
    timerValue = tic;
    pd_log(logFile, 'info', 'run', ['Starting ', request.analysisType]);
    try
        pd_execution_checkpoint(request, 0, 'Resolving input file');
        filePath = pd_resolve_input_file(request);
        pd_log(logFile, 'info', 'input', ['Resolved ', filePath]);
        preflight = pd_preflight_input(request, filePath);
        for warningIndex = 1:numel(preflight.warnings)
            pd_log(logFile, 'warning', 'preflight', preflight.warnings{warningIndex});
        end
        pd_execution_checkpoint(request, 0.05, 'Input preflight passed');
        selArgs = pd_build_selection_args(request);
        analysisProgress = @(fraction, message) forwardAnalysisProgress( ...
            request, fraction, message);

        switch request.analysisType
            case 'chunk'
                result = analyze_chunk_field(filePath, ...
                    'ChunkDim', request.chunk.dimension, ...
                    'Variable', request.chunk.variable, ...
                    selArgs{:}, ...
                    'dV', request.chunk.dV, ...
                    'CoordScale', request.chunk.coordScale, ...
                    'CoordRangeX', request.chunk.coordRangeX, ...
                    'CoordRangeY', request.chunk.coordRangeY, ...
                    'CoordRangeZ', request.chunk.coordRangeZ, ...
                    'GradientVariable', request.chunk.gradientVariable, ...
                    'GradientSmoothLevel', request.chunk.gradientSmoothLevel, ...
                    'StrainRateVelocityComponent', request.chunk.strainRateVelocityComponent, ...
                    'StrainRateDensityVariable', request.chunk.strainRateDensityVariable, ...
                    'CancelCallback', request.execution.cancelCallback, ...
                    'ProgressCallback', analysisProgress, ...
                    'PlotOptions', request.chunk.plotOptions, ...
                    'DoPlot', request.makePlots);
            case 'cluster'
                options = pd_upsert_option(request.analysisOptions, 'MakePlots', request.makePlots);
                options = pd_upsert_option(options, 'CancelCallback', request.execution.cancelCallback);
                options = pd_upsert_option(options, 'ProgressCallback', analysisProgress);
                result = cluster_postprocess(filePath, selArgs{:}, options{:});
            case 'vx'
                options = pd_upsert_option(request.analysisOptions, 'MakePlot', request.makePlots);
                options = pd_upsert_option(options, 'CancelCallback', request.execution.cancelCallback);
                options = pd_upsert_option(options, 'ProgressCallback', analysisProgress);
                result = vx_chunk_cumulative(filePath, selArgs{:}, options{:});
            case 'massx'
                options = pd_upsert_option(request.analysisOptions, 'MakePlot', request.makePlots);
                options = pd_upsert_option(options, 'CancelCallback', request.execution.cancelCallback);
                options = pd_upsert_option(options, 'ProgressCallback', analysisProgress);
                result = mass_x_cumulative(filePath, selArgs{:}, options{:});
            case 'network2d'
                options = pd_upsert_option(request.analysisOptions, 'MakePlots', request.makePlots);
                options = pd_upsert_option(options, 'CancelCallback', request.execution.cancelCallback);
                options = pd_upsert_option(options, 'ProgressCallback', analysisProgress);
                result = analyze_chunk_network2d(filePath, selArgs{:}, ...
                    'CoordScale', request.chunk.coordScale, options{:});
        end

        pd_execution_checkpoint(request, 1, 'Analysis complete');
        result = pd_add_result_metadata(result, request.analysisType);
        if request.makePlots
            pd_style_result_figures(result, request.plotOptions);
        end
        result.preflight = preflight;
        result.run = struct('startedAt', startedAt, 'durationSeconds', toc(timerValue), ...
            'logFile', logFile, 'request', pd_request_manifest(request));
        result = pd_reduce_result(result, request.resultLevel);
        pd_log(logFile, 'info', 'run', sprintf('Completed %s in %.3f s', ...
            request.analysisType, result.run.durationSeconds));
    catch err
        if strcmp(err.identifier, 'postdata:UserCancelled')
            pd_log(logFile, 'warning', 'run', 'Analysis cancelled by user.');
        else
            pd_log(logFile, 'error', 'run', sprintf('%s | %s', err.identifier, err.message));
        end
        rethrow(err);
    end
end

function forwardAnalysisProgress(request, fraction, message)
    fraction = max(0, min(1, fraction));
    request.execution.progressCallback(0.05 + 0.90 * fraction, message);
end
