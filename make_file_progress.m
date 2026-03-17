function tracker = make_file_progress(progressMode, operationName, filePath)
%MAKE_FILE_PROGRESS Create a progress reporter for file-reading operations.
% MATLAB R2016b compatible.
%
% tracker = make_file_progress(progressMode, operationName, filePath)
%
% Supported progressMode values:
%   'auto'    : use waitbar on desktop MATLAB, otherwise console output
%   'off'     : disable progress output
%   'waitbar' : force waitbar
%   'console' : force command-window progress
%
% Returned struct fields:
%   tracker.mode   : resolved mode after auto/fallback
%   tracker.update : function handle accepting fraction in [0,1]
%   tracker.finish : function handle, emits completion and closes UI
%   tracker.close  : function handle, closes UI without forcing completion

    if nargin < 1 || isempty(progressMode)
        progressMode = 'auto';
    end
    if nargin < 2 || isempty(operationName)
        operationName = 'Reading file';
    end
    if nargin < 3
        filePath = '';
    end

    progressMode = toChar(progressMode);
    operationName = toChar(operationName);
    filePath = toChar(filePath);

    requestedMode = lower(strtrim(progressMode));
    validateMode(requestedMode);

    displayName = shortFileName(filePath);
    if isempty(displayName)
        label = operationName;
    else
        label = sprintf('%s: %s', operationName, displayName);
    end

    resolvedMode = resolveMode(requestedMode);
    wb = [];
    lastWaitbarPct = -1;
    consoleStarted = false;
    consoleNextPct = 5;
    consoleLastPct = 0;

    if strcmp(resolvedMode, 'waitbar')
        try
            wb = waitbar(0, sprintf('%s\n0%%', label), 'Name', label);
            lastWaitbarPct = 0;
        catch
            resolvedMode = 'console';
        end
    end

    if strcmp(resolvedMode, 'console')
        emitConsoleStart();
    end

    tracker = struct();
    tracker.mode = resolvedMode;
    tracker.update = @updateProgress;
    tracker.finish = @finishProgress;
    tracker.close = @closeProgress;

    function updateProgress(fraction)
        fraction = sanitizeFraction(fraction);
        pct = floor(100 * fraction);

        if strcmp(resolvedMode, 'waitbar')
            if isempty(wb) || ~ishghandle(wb)
                return;
            end
            if pct <= lastWaitbarPct && pct < 100
                return;
            end
            try
                waitbar(fraction, wb, sprintf('%s\n%d%%', label, pct));
                lastWaitbarPct = pct;
            catch
                closeWaitbar();
                resolvedMode = 'console';
                emitConsoleStart();
                emitConsoleProgress(pct);
            end
            return;
        end

        if strcmp(resolvedMode, 'console')
            emitConsoleProgress(pct);
        end
    end

    function finishProgress()
        updateProgress(1.0);
        if strcmp(resolvedMode, 'console') && consoleLastPct < 100
            fprintf('%s (100%%)\n', label);
            consoleLastPct = 100;
            consoleNextPct = 105;
        end
        closeProgress();
    end

    function closeProgress()
        if strcmp(resolvedMode, 'waitbar')
            closeWaitbar();
        end
    end

    function emitConsoleStart()
        if consoleStarted
            return;
        end
        fprintf('%s (0%%)\n', label);
        consoleStarted = true;
        consoleLastPct = 0;
        consoleNextPct = 5;
    end

    function emitConsoleProgress(pct)
        if ~consoleStarted
            emitConsoleStart();
        end
        pct = max(0, min(100, pct));
        while consoleNextPct <= pct
            fprintf('%s (%d%%)\n', label, consoleNextPct);
            consoleLastPct = consoleNextPct;
            consoleNextPct = consoleNextPct + 5;
        end
    end

    function closeWaitbar()
        if ~isempty(wb) && ishghandle(wb)
            try
                close(wb);
            catch
                % ignore close failures
            end
        end
        wb = [];
    end
end

function mode = resolveMode(requestedMode)
    if strcmp(requestedMode, 'auto')
        if usejava('desktop') && usejava('awt')
            mode = 'waitbar';
        else
            mode = 'console';
        end
        return;
    end
    mode = requestedMode;
end

function validateMode(mode)
    validModes = {'auto', 'off', 'waitbar', 'console'};
    if ~any(strcmp(mode, validModes))
        error('make_file_progress:BadMode', ...
            'ProgressMode must be auto/off/waitbar/console.');
    end
end

function s = toChar(v)
    if isstring(v)
        s = char(v);
    else
        s = v;
    end
end

function name = shortFileName(filePath)
    if isempty(filePath)
        name = '';
        return;
    end
    [~, base, ext] = fileparts(filePath);
    if isempty(base) && isempty(ext)
        name = filePath;
    else
        name = [base, ext];
    end
end

function fraction = sanitizeFraction(fraction)
    if nargin < 1 || isempty(fraction) || ~isscalar(fraction) || ~isfinite(fraction)
        fraction = 0;
        return;
    end
    fraction = max(0, min(1, fraction));
end
