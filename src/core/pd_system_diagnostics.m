function report = pd_system_diagnostics()
%PD_SYSTEM_DIAGNOSTICS Inspect runtime and installation prerequisites.
% The check is read-only and is safe to call from scripts or the GUI.

    coreDir = fileparts(mfilename('fullpath'));
    sourceDir = fileparts(coreDir);
    projectRoot = fileparts(sourceDir);
    info = pd_version();
    requiredModules = {'analysis','core','io','plot','export'};
    moduleAvailable = false(size(requiredModules));
    for i = 1:numel(requiredModules)
        moduleAvailable(i) = exist(fullfile(sourceDir, requiredModules{i}), ...
            'dir') == 7;
    end

    report = struct();
    report.software = info;
    report.projectRoot = projectRoot;
    report.matlabVersion = version;
    report.matlabRelease = version('-release');
    report.isCompatibleMatlab = ~verLessThan('matlab', '9.1');
    report.desktopAvailable = usejava('desktop') && usejava('awt');
    hasSkeleton = exist('bwskel', 'file') == 2 || exist('bwmorph', 'file') == 2;
    report.imageProcessingToolboxAvailable = ...
        exist('bwdist', 'file') == 2 && hasSkeleton;
    report.requiredModules = requiredModules;
    report.moduleAvailable = moduleAvailable;
    report.sampleDataAvailable = exist(fullfile(projectRoot, 'fixtures', 'sample'), ...
        'dir') == 7;
    report.defaultOutputDirectory = fullfile(projectRoot, 'outputs');
    report.logFile = fullfile(prefdir, 'POST_DATA2', 'logs', 'postdata.log');
    report.warnings = {};
    if ~report.isCompatibleMatlab
        report.warnings{end + 1} = ...
            'MATLAB R2016b (9.1) or newer is required.'; %#ok<AGROW>
    end
    if ~all(moduleAvailable)
        missing = requiredModules(~moduleAvailable);
        report.warnings{end + 1} = ...
            ['Missing source modules: ', strjoin(missing, ', ')]; %#ok<AGROW>
    end
    if ~report.sampleDataAvailable
        report.warnings{end + 1} = ...
            'Sample fixtures are missing; built-in examples are unavailable.'; %#ok<AGROW>
    end
    if ~report.imageProcessingToolboxAvailable
        report.warnings{end + 1} = ...
            ['Image Processing Toolbox is unavailable; optional network ' ...
             'skeleton/thickness metrics will remain disabled.']; %#ok<AGROW>
    end
    report.passed = report.isCompatibleMatlab && all(moduleAvailable);
end
