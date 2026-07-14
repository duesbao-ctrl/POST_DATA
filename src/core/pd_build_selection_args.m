function args = pd_build_selection_args(request)
%PD_BUILD_SELECTION_ARGS Convert request selection into analyzer arguments.

    selection = request.selection;
    args = {'SelectBy', selection.mode, ...
            'SlurmPath', selection.slurmPath, ...
            'SlurmModuleIndex', selection.slurmModuleIndex, ...
            'ProgressMode', request.progressMode};
    switch lower(strtrim(selection.mode))
        case 'index'
            args = [args, {'Index', selection.value}];
        case 'timestep'
            args = [args, {'TimeStep', selection.value}];
        case 'time'
            args = [args, {'Time', selection.value}];
    end
end
