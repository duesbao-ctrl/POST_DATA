function out = run_analysis(taskType, varargin)
%RUN_ANALYSIS Unified entry for three analysis categories.
%   out = RUN_ANALYSIS('chunk', ...)
%   out = RUN_ANALYSIS('cluster', ...)
%   out = RUN_ANALYSIS('vx', ...)
%   out = RUN_ANALYSIS('network2d', ...)
%
% Common selector params (for all 3 categories):
%   SelectBy: 'Index' | 'TimeStep' | 'Time'
%   Index, TimeStep, Time, SlurmPath, SlurmModuleIndex, ProgressMode
%
% taskType='chunk' required params:
%   ChunkDim : '1d' or '2d'
%   Variable : e.g. c_rho / T / velocity / pressure / density / Sxx ... / vonMisesS
% Optional:
%   ChunkFile (.txt), dV, PlotOptions, DoPlot, CoordScale, CoordRangeX, CoordRangeY
%
% taskType='cluster' optional params:
%   ClusterFile (.txt), ClusterOptions
%
% taskType='vx' optional params:
%   VxFile (.txt), VxOptions
%
% taskType='network2d' optional params:
%   ChunkFile (.txt), NetworkOptions

    p = inputParser;
    p.addRequired('taskType', @isTextScalar);

    p.addParameter('BaseDir', pwd, @isTextScalar);

    p.addParameter('SelectBy', 'Index', @isTextScalar);
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);

    p.addParameter('ChunkDim', '2d', @isTextScalar);
    p.addParameter('ChunkFile', '', @isTextScalar);
    p.addParameter('Variable', 'c_rho', @isTextScalar);
    p.addParameter('dV', [], @isnumeric);
    p.addParameter('DoPlot', true, @islogical);
    p.addParameter('PlotOptions', {}, @iscell);
    p.addParameter('CoordScale', 1, @isnumeric);
    p.addParameter('CoordRangeX', [], @isnumeric);
    p.addParameter('CoordRangeY', [], @isnumeric);

    p.addParameter('ClusterFile', '', @isTextScalar);
    p.addParameter('ClusterOptions', {}, @iscell);

    p.addParameter('VxFile', '', @isTextScalar);
    p.addParameter('VxOptions', {}, @iscell);

    p.addParameter('NetworkOptions', {}, @iscell);

    p.parse(taskType, varargin{:});
    opt = p.Results;
    taskType = toChar(taskType);
    opt.BaseDir = toChar(opt.BaseDir);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);
    opt.ChunkDim = toChar(opt.ChunkDim);
    opt.ChunkFile = toChar(opt.ChunkFile);
    opt.Variable = toChar(opt.Variable);
    opt.ClusterFile = toChar(opt.ClusterFile);
    opt.VxFile = toChar(opt.VxFile);

    t = lower(strtrim(taskType));

    % Auto slurm path for Time mode when missing.
    if strcmpi(opt.SelectBy, 'Time') && isempty(opt.SlurmPath)
        opt.SlurmPath = get_slurm_txt_fullpath(opt.BaseDir);
    end

    selArgs = {'SelectBy', opt.SelectBy, ...
               'Index', opt.Index, ...
               'TimeStep', opt.TimeStep, ...
               'Time', opt.Time, ...
               'SlurmPath', opt.SlurmPath, ...
               'SlurmModuleIndex', opt.SlurmModuleIndex, ...
               'ProgressMode', opt.ProgressMode};

    switch t
        case 'chunk'
            chunkFile = resolveChunkFile(opt.BaseDir, opt.ChunkFile, opt.ChunkDim);
            out = analyze_chunk_field(chunkFile, ...
                'ChunkDim', opt.ChunkDim, ...
                'Variable', opt.Variable, ...
                selArgs{:}, ...
                'dV', opt.dV, ...
                'DoPlot', opt.DoPlot, ...
                'PlotOptions', opt.PlotOptions, ...
                'CoordScale', opt.CoordScale, ...
                'CoordRangeX', opt.CoordRangeX, ...
                'CoordRangeY', opt.CoordRangeY);

        case 'cluster'
            clusterFile = resolveByPattern(opt.BaseDir, opt.ClusterFile, 'cluster_chunk*.txt', true, 'ClusterFile');
            out = cluster_postprocess(clusterFile, selArgs{:}, opt.ClusterOptions{:});

        case 'vx'
            vxFile = resolveByPattern(opt.BaseDir, opt.VxFile, 'vx_chunk*.txt', true, 'VxFile');
            out = vx_chunk_cumulative(vxFile, selArgs{:}, opt.VxOptions{:});

        case 'network2d'
            if ~strcmpi(opt.ChunkDim, '2d')
                error('run_analysis:BadNetworkChunkDim', ...
                    'taskType=''network2d'' requires ChunkDim=''2d''.');
            end
            chunkFile = resolveChunkFile(opt.BaseDir, opt.ChunkFile, '2d');
            out = analyze_chunk_network2d(chunkFile, selArgs{:}, ...
                'CoordScale', opt.CoordScale, ...
                opt.NetworkOptions{:});

        otherwise
            error('run_analysis:BadTaskType', 'taskType must be chunk/cluster/vx/network2d.');
    end
end

function tf = isTextScalar(v)
    tf = ischar(v) || (isstring(v) && isscalar(v));
end

function s = toChar(v)
    if isstring(v)
        s = char(v);
    else
        s = v;
    end
end

function p = resolveChunkFile(baseDir, chunkFile, chunkDim)
    if ~isempty(chunkFile)
        p = resolvePreferredTextFile(baseDir, chunkFile, 'ChunkFile');
        return;
    end

    dimTag = lower(strtrim(chunkDim));
    if strcmp(dimTag, '1d')
        p = resolveByPattern(baseDir, '', 'bin1d*.txt');
    elseif strcmp(dimTag, '2d')
        p = resolveByPattern(baseDir, '', 'bin2d*.txt');
    else
        error('run_analysis:BadChunkDim', 'ChunkDim must be ''1d'' or ''2d''.');
    end
end

function p = resolveByPattern(baseDir, preferred, pattern, requireTxt, argName)
    if nargin < 4
        requireTxt = false;
    end
    if nargin < 5
        argName = 'File';
    end
    if ~isempty(preferred)
        p = resolvePreferredFile(baseDir, preferred, argName, requireTxt);
        return;
    end

    files = dir(fullfile(baseDir, pattern));
    files = files(~[files.isdir]);

    if isempty(files)
        error('run_analysis:NoMatch', 'No file matches %s in %s', pattern, baseDir);
    end
    if numel(files) > 1
        names = {files.name};
        error('run_analysis:MultipleMatches', ...
            'Multiple files match %s in %s: %s', pattern, baseDir, strjoin(names, ', '));
    end

    p = fullfile(baseDir, files(1).name);
end

function p = resolvePreferredTextFile(baseDir, preferred, argName)
    p = resolvePreferredFile(baseDir, preferred, argName, true);
end

function p = resolvePreferredFile(baseDir, preferred, argName, requireTxt)
    if nargin < 4
        requireTxt = false;
    end

    if requireTxt && ~hasTxtSuffix(preferred)
        error('run_analysis:MissingTxtSuffix', ...
            '%s must include the .txt suffix to avoid matching cache files: %s', ...
            argName, preferred);
    end

    if isAbsolutePath(preferred)
        p = preferred;
    else
        p = fullfile(baseDir, preferred);
    end
    if ~exist(p, 'file')
        error('run_analysis:PreferredFileNotFound', 'File not found: %s', p);
    end
end

function tf = hasTxtSuffix(p)
    tf = ~isempty(regexpi(p, '\.txt$', 'once'));
end

function tf = isAbsolutePath(p)
    tf = false;
    if isempty(p)
        return;
    end
    if numel(p) >= 2 && p(2) == ':'
        tf = true;
        return;
    end
    if numel(p) >= 2 && strcmp(p(1:2), '\\')
        tf = true;
    end
end
