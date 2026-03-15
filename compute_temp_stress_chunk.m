function out = compute_temp_stress_chunk(chunkPath, varargin)
%COMPUTE_TEMP_STRESS_CHUNK Compute temperature and Cauchy stress from chunk data.
% MATLAB R2016b compatible.
%
% Required columns (or aliases):
%   Ncount, v_mass, v_mvx, v_mvy, v_mvz, c_ke,
%   v_mvvxx, v_mvvyy, v_mvvzz, v_mvvxy, v_mvvxz, v_mvvyz,
%   c_virial[1..6] (after makeValidName typically c_virial1..c_virial6).
%
% Selection modes:
%   SelectBy = 'Index' | 'TimeStep' | 'Time'
%
% dV:
%   pass 'dV' explicitly as a positive scalar.

    p = inputParser;
    p.addRequired('chunkPath', @isTextScalar);

    p.addParameter('SelectBy', 'Index', @isTextScalar);      % Index | TimeStep | Time
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);

    p.addParameter('dV', [], @isnumeric);

    p.addParameter('kB', 8.617343e-5, @isnumeric);
    p.addParameter('bar2GPa', 1e-4, @isnumeric);
    p.addParameter('mvv2e', 1.0364269e-4, @isnumeric);
    p.addParameter('nktv2p', 1.6021765e6, @isnumeric);

    p.parse(chunkPath, varargin{:});
    opt = p.Results;
    chunkPath = toChar(chunkPath);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex};

    step = read_chunk_step_fast(chunkPath, selectorArgs{:});
    stepIdx = step.stepIndex;
    selectionInfo = buildSelectionInfo(step, opt);
    D = step.data;
    col = step.colIndex;

    if isempty(opt.dV) || ~isscalar(opt.dV) || ~isfinite(opt.dV) || opt.dV <= 0
        error('compute_temp_stress_chunk:BadDV', ...
            'dV is required and must be a positive finite scalar.');
    end
    dV = opt.dV;

    idx.N    = findCol(col, {'Ncount'});
    idx.M    = findCol(col, {'v_mass','mass'});
    idx.mvx  = findCol(col, {'v_mvx','mvx'});
    idx.mvy  = findCol(col, {'v_mvy','mvy'});
    idx.mvz  = findCol(col, {'v_mvz','mvz'});
    idx.KE   = findCol(col, {'c_ke','ke'});

    idx.mvvxx = tryFindCol(col, {'v_mvvxx','mvvxx'});
    idx.mvvyy = tryFindCol(col, {'v_mvvyy','mvvyy'});
    idx.mvvzz = tryFindCol(col, {'v_mvvzz','mvvzz'});
    idx.mvvxy = tryFindCol(col, {'v_mvvxy','mvvxy'});
    idx.mvvxz = tryFindCol(col, {'v_mvvxz','mvvxz'});
    idx.mvvyz = tryFindCol(col, {'v_mvvyz','mvvyz'});

    idx.Wxx = tryFindCol(col, {'c_virial[1]','c_virial_1_','c_virial_1','c_virial1'});
    idx.Wyy = tryFindCol(col, {'c_virial[2]','c_virial_2_','c_virial_2','c_virial2'});
    idx.Wzz = tryFindCol(col, {'c_virial[3]','c_virial_3_','c_virial_3','c_virial3'});
    idx.Wxy = tryFindCol(col, {'c_virial[4]','c_virial_4_','c_virial_4','c_virial4'});
    idx.Wxz = tryFindCol(col, {'c_virial[5]','c_virial_5_','c_virial_5','c_virial5'});
    idx.Wyz = tryFindCol(col, {'c_virial[6]','c_virial_6_','c_virial_6','c_virial6'});

    N = D(:, idx.N);
    M = D(:, idx.M);
    mvx = D(:, idx.mvx);
    mvy = D(:, idx.mvy);
    mvz = D(:, idx.mvz);
    KE  = D(:, idx.KE);

    mvvxx = getOrEmpty(D, idx.mvvxx);
    mvvyy = getOrEmpty(D, idx.mvvyy);
    mvvzz = getOrEmpty(D, idx.mvvzz);
    mvvxy = getOrEmpty(D, idx.mvvxy);
    mvvxz = getOrEmpty(D, idx.mvvxz);
    mvvyz = getOrEmpty(D, idx.mvvyz);

    Wxx = getOrEmpty(D, idx.Wxx);
    Wyy = getOrEmpty(D, idx.Wyy);
    Wzz = getOrEmpty(D, idx.Wzz);
    Wxy = getOrEmpty(D, idx.Wxy);
    Wxz = getOrEmpty(D, idx.Wxz);
    Wyz = getOrEmpty(D, idx.Wyz);

    vx = safeDiv(mvx, M);
    vy = safeDiv(mvy, M);
    vz = safeDiv(mvz, M);

    KE_bulk = 0.5 .* opt.mvv2e .* M .* (vx.^2 + vy.^2 + vz.^2);
    KE_th   = KE + KE_bulk - opt.mvv2e .* (vx.*mvx + vy.*mvy + vz.*mvz);

    T = safeDiv(2 .* KE_th, 3 .* N .* opt.kB);

    Sxx_th = [];
    Syy_th = [];
    Szz_th = [];
    Sxy_th = [];
    Sxz_th = [];
    Syz_th = [];

    sigma_xx = [];
    sigma_yy = [];
    sigma_zz = [];
    sigma_xy = [];
    sigma_xz = [];
    sigma_yz = [];

    computed = struct('xx', false, 'yy', false, 'zz', false, ...
                      'xy', false, 'xz', false, 'yz', false);

    if ~isempty(mvvxx) && ~isempty(Wxx)
        Sxx_th = opt.mvv2e .* (mvvxx + M .* vx.^2 - 2.0 .* vx .* mvx) .* opt.nktv2p;
        sigma_xx = (Wxx - Sxx_th) ./ dV .* opt.bar2GPa;
        computed.xx = true;
    end
    if ~isempty(mvvyy) && ~isempty(Wyy)
        Syy_th = opt.mvv2e .* (mvvyy + M .* vy.^2 - 2.0 .* vy .* mvy) .* opt.nktv2p;
        sigma_yy = (Wyy - Syy_th) ./ dV .* opt.bar2GPa;
        computed.yy = true;
    end
    if ~isempty(mvvzz) && ~isempty(Wzz)
        Szz_th = opt.mvv2e .* (mvvzz + M .* vz.^2 - 2.0 .* vz .* mvz) .* opt.nktv2p;
        sigma_zz = (Wzz - Szz_th) ./ dV .* opt.bar2GPa;
        computed.zz = true;
    end
    if ~isempty(mvvxy) && ~isempty(Wxy)
        Sxy_th = opt.mvv2e .* (mvvxy + M .* vx .* vy - vx .* mvy - vy .* mvx) .* opt.nktv2p;
        sigma_xy = (Wxy - Sxy_th) ./ dV .* opt.bar2GPa;
        computed.xy = true;
    end
    if ~isempty(mvvxz) && ~isempty(Wxz)
        Sxz_th = opt.mvv2e .* (mvvxz + M .* vx .* vz - vx .* mvz - vz .* mvx) .* opt.nktv2p;
        sigma_xz = (Wxz - Sxz_th) ./ dV .* opt.bar2GPa;
        computed.xz = true;
    end
    if ~isempty(mvvyz) && ~isempty(Wyz)
        Syz_th = opt.mvv2e .* (mvvyz + M .* vy .* vz - vy .* mvz - vz .* mvy) .* opt.nktv2p;
        sigma_yz = (Wyz - Syz_th) ./ dV .* opt.bar2GPa;
        computed.yz = true;
    end

    out = struct();
    out.filePath = chunkPath;
    out.selection = selectionInfo;
    out.stepIndex = stepIdx;
    out.timestep = step.timestep;
    out.dV = dV;

    out.N = N;
    out.M = M;
    out.vx = vx;
    out.vy = vy;
    out.vz = vz;
    out.KE = KE;
    out.KE_bulk = KE_bulk;
    out.KE_th = KE_th;
    out.T = T;

    out.S_th = struct('xx', Sxx_th, 'yy', Syy_th, 'zz', Szz_th, ...
                      'xy', Sxy_th, 'xz', Sxz_th, 'yz', Syz_th);
    out.sigma = struct('xx', sigma_xx, 'yy', sigma_yy, 'zz', sigma_zz, ...
                       'xy', sigma_xy, 'xz', sigma_xz, 'yz', sigma_yz);
    out.computed = computed;

    out.table = table(T, 'VariableNames', {'T'});
    if computed.xx, out.table.sigma_xx = sigma_xx; end
    if computed.yy, out.table.sigma_yy = sigma_yy; end
    if computed.zz, out.table.sigma_zz = sigma_zz; end
    if computed.xy, out.table.sigma_xy = sigma_xy; end
    if computed.xz, out.table.sigma_xz = sigma_xz; end
    if computed.yz, out.table.sigma_yz = sigma_yz; end
end

function idx = findCol(colMap, candidates)
    idx = [];
    for i = 1:numel(candidates)
        name = matlab.lang.makeValidName(candidates{i});
        if isfield(colMap, name)
            idx = colMap.(name);
            return;
        end
    end
    error('compute_temp_stress_chunk:MissingColumn', ...
        'Missing required column. Tried: %s', strjoin(candidates, ', '));
end

function idx = tryFindCol(colMap, candidates)
    idx = [];
    for i = 1:numel(candidates)
        name = matlab.lang.makeValidName(candidates{i});
        if isfield(colMap, name)
            idx = colMap.(name);
            return;
        end
    end
end

function v = getOrEmpty(D, idx)
    if isempty(idx)
        v = [];
    else
        v = D(:, idx);
    end
end

function z = safeDiv(a, b)
    z = nan(size(a));
    m = (b ~= 0) & isfinite(b);
    z(m) = a(m) ./ b(m);
end

function info = buildSelectionInfo(step, opt)
    mode = lower(strtrim(opt.SelectBy));
    info = struct();
    info.mode = mode;
    info.stepIndex = step.stepIndex;
    info.mappedTimeStep = step.timestep;

    switch mode
        case 'index'
            info.requestedIndex = opt.Index;

        case 'timestep'
            info.requestedTimeStep = opt.TimeStep;

        case 'time'
            info.requestedTime = opt.Time;

        otherwise
            % read_chunk_step_fast has already validated selector.
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
