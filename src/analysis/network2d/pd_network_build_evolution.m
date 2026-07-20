function output = pd_network_build_evolution(chunkFile, options, snapshotFunction)
%PD_NETWORK_BUILD_EVOLUTION Select, execute, and aggregate network snapshots.

    axisName = resolveAxis(options);
    index = loadIndex(chunkFile, options);
    stepIndexAll = (1:numel(index.timesteps)).';
    timestepAll = index.timesteps(:);
    timeAll = nan(size(timestepAll));
    if strcmp(axisName, 'time')
        timeAll = resolvePhysicalTime(index, timestepAll, options, chunkFile);
        axisValueAll = timeAll;
    elseif ~isempty(options.SlurmPath)
        timeAll = resolvePhysicalTime(index, timestepAll, options, chunkFile);
        axisValueAll = timestepAll;
    elseif strcmp(axisName, 'timestep')
        axisValueAll = timestepAll;
    else
        axisValueAll = stepIndexAll;
    end

    selected = selectSnapshots(axisValueAll, options.EvolutionRange, options.EvolutionStride);
    output = pd_network_empty_evolution();
    output.axis = axisName;
    if isempty(selected), return; end
    output = allocateOutput(output, axisValueAll, stepIndexAll, timestepAll, timeAll, selected);

    count = numel(selected);
    for i = 1:count
        drawnow;
        if options.CancelCallback()
            error('postdata:UserCancelled', 'Evolution analysis cancelled by user.');
        end
        options.ProgressCallback((i - 1) / count, ...
            sprintf('Evolution snapshot %d of %d', i, count));
        step = read_chunk_step_fast(chunkFile, 'SelectBy', 'Index', ...
            'Index', output.stepIndex(i), 'ProgressMode', 'off', ...
            'CancelCallback', options.CancelCallback);
        stats = snapshotFunction(step);
        output = assignSnapshot(output, stats, i);
    end
    options.ProgressCallback(1, 'Evolution analysis complete');
    output = buildTransitions(output);
end

function times = resolvePhysicalTime(index, timesteps, options, chunkFile)
    if isfield(index, 'physicalTimes')
        embedded = index.physicalTimes(:);
        if numel(embedded) == numel(timesteps) && all(isfinite(embedded))
            times = embedded;
            return;
        end
    end
    times = mapPhysicalTime(timesteps, options, chunkFile);
end

function selected = selectSnapshots(axisValues, range, stride)
    mask = true(size(axisValues));
    if ~isempty(range)
        mask = isfinite(axisValues) & axisValues >= range(1) & axisValues <= range(2);
    end
    selected = find(mask);
    selected = selected(1:max(1, round(stride)):end);
end

function output = allocateOutput(output, axisValues, stepIndices, timesteps, times, selected)
    count = numel(selected);
    output.axisValue = axisValues(selected).';
    output.stepIndex = stepIndices(selected).';
    output.timestep = timesteps(selected).';
    output.time = times(selected).';
    output.geometry.phi = nan(1, count);
    output.geometry.interfaceLength = nan(1, count);
    output.geometry.specificInterface = nan(1, count);
    phases = {'pore','matrix'};
    for i = 1:numel(phases)
        phase = phases{i};
        output.topology.(phase).beta0 = nan(1, count);
        output.topology.(phase).beta1 = nan(1, count);
        output.topology.(phase).chi = nan(1, count);
        output.connectivity.(phase).largestFraction = nan(1, count);
        output.connectivity.(phase).percolatesX = false(1, count);
        output.connectivity.(phase).percolatesY = false(1, count);
    end
    output.thickness.matrix.min = nan(1, count);
    output.thickness.matrix.mean = nan(1, count);
    output.thickness.matrix.p1 = nan(1, count);
    output.thickness.matrix.p5 = nan(1, count);
    output.fragmentation.matrix.count = nan(1, count);
    output.fragmentation.matrix.largestFraction = nan(1, count);
end

function output = assignSnapshot(output, stats, index)
    output.geometry.phi(index) = stats.geometry.phi;
    output.geometry.interfaceLength(index) = stats.geometry.interfaceLength;
    output.geometry.specificInterface(index) = stats.geometry.specificInterface;
    phases = {'pore','matrix'};
    for i = 1:numel(phases)
        phase = phases{i};
        output.topology.(phase).beta0(index) = stats.topology.(phase).beta0;
        output.topology.(phase).beta1(index) = stats.topology.(phase).beta1;
        output.topology.(phase).chi(index) = stats.topology.(phase).chi;
        output.connectivity.(phase).largestFraction(index) = stats.connectivity.(phase).largestFraction;
        output.connectivity.(phase).percolatesX(index) = stats.connectivity.(phase).percolatesX;
        output.connectivity.(phase).percolatesY(index) = stats.connectivity.(phase).percolatesY;
    end
    output.thickness.matrix.min(index) = stats.thickness.matrix.min;
    output.thickness.matrix.mean(index) = stats.thickness.matrix.mean;
    output.thickness.matrix.p1(index) = stats.thickness.matrix.p1;
    output.thickness.matrix.p5(index) = stats.thickness.matrix.p5;
    output.fragmentation.matrix.count(index) = stats.fragmentation.matrix.count;
    output.fragmentation.matrix.largestFraction(index) = stats.fragmentation.matrix.largestFraction;
end

function output = buildTransitions(output)
    output.transition.deltaPhi = leadingDiff(output.geometry.phi);
    output.transition.deltaSpecificInterface = leadingDiff(output.geometry.specificInterface);
    phases = {'pore','matrix'};
    for i = 1:numel(phases)
        phase = phases{i};
        output.transition.(phase).deltaBeta0 = leadingDiff(output.topology.(phase).beta0);
        output.transition.(phase).deltaBeta1 = leadingDiff(output.topology.(phase).beta1);
        output.transition.(phase).deltaLargestFraction = ...
            leadingDiff(output.connectivity.(phase).largestFraction);
    end
end

function axisName = resolveAxis(options)
    axisName = options.EvolutionSelectBy;
    if isempty(axisName), axisName = lower(strtrim(options.SelectBy)); end
    if ~any(strcmp(axisName, {'timestep','time'})), axisName = 'index'; end
end

function index = loadIndex(chunkFile, options)
    % Use the in-memory reader contract. Evolution must also work when the
    % input directory is read-only and no sidecar cache can be created.
    step = read_chunk_step_fast(chunkFile, 'SelectBy', 'Index', 'Index', 1, ...
        'ProgressMode', 'off', 'CancelCallback', options.CancelCallback);
    index = step.index;
end

function times = mapPhysicalTime(timesteps, options, chunkFile)
    slurmPath = options.SlurmPath;
    if isempty(slurmPath), slurmPath = get_slurm_txt_fullpath(fileparts(chunkFile)); end
    if isempty(slurmPath) || ~exist(slurmPath, 'file')
        error('analyze_chunk_network2d:MissingEvolutionSlurmPath', ...
            'EvolutionSelectBy=Time requires a valid SlurmPath.');
    end
    slurm = read_slurm_stepcpu(slurmPath);
    moduleIndex = round(options.SlurmModuleIndex);
    if moduleIndex < 1 || moduleIndex > numel(slurm.modules)
        error('analyze_chunk_network2d:BadEvolutionSlurmModule', ...
            'SlurmModuleIndex out of range for evolution time mapping.');
    end
    columns = slurm.modules(moduleIndex).colIndex;
    if ~isfield(columns, 'Time') || ~isfield(columns, 'Step')
        error('analyze_chunk_network2d:BadEvolutionSlurmColumns', ...
            'Slurm block must contain Time and Step columns for evolution mapping.');
    end
    slurmTimes = slurm.modules(moduleIndex).data(:, columns.Time);
    slurmSteps = slurm.modules(moduleIndex).data(:, columns.Step);
    valid = isfinite(slurmTimes) & isfinite(slurmSteps);
    slurmTimes = slurmTimes(valid);
    slurmSteps = slurmSteps(valid);
    times = nan(size(timesteps));
    for i = 1:numel(timesteps)
        [~, nearest] = min(abs(slurmSteps - timesteps(i)));
        times(i) = slurmTimes(nearest);
    end
end

function difference = leadingDiff(values)
    values = values(:).';
    difference = nan(size(values));
    if numel(values) >= 2, difference(2:end) = diff(values); end
end
