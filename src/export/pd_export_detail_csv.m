function files = pd_export_detail_csv(result, basePath)
%PD_EXPORT_DETAIL_CSV Export analysis-specific vectors and tables.

    files = {};
    switch result.analysisType
        case 'chunk'
            path = [basePath, '_field.csv'];
            if isempty(result.y)
                pd_write_numeric_csv(path, {'x','value'}, [result.x(:), result.value(:)]);
            else
                pd_write_numeric_csv(path, {'x','y','value'}, ...
                    [result.x(:), result.y(:), result.value(:)]);
            end
            files{end + 1} = path;
        case 'cluster'
            path = [basePath, '_clusters.csv'];
            index = (1:numel(result.diameter)).';
            pd_write_numeric_csv(path, {'index','equivalentDiameter'}, ...
                [index, result.diameter(:)]);
            files{end + 1} = path;
        case 'vx'
            path = [basePath, '_velocity.csv'];
            headers = {'velocity'};
            data = result.velocity(:);
            for i = 1:numel(result.densityVars)
                name = regexprep(result.densityVars{i}, '[^A-Za-z0-9_]', '_');
                headers{end + 1} = ['density_', name]; %#ok<AGROW>
                headers{end + 1} = ['cumulative_', name]; %#ok<AGROW>
                data(:, end + 1) = result.density(:, i); %#ok<AGROW>
                data(:, end + 1) = result.cumulativeDensity(:, i); %#ok<AGROW>
            end
            pd_write_numeric_csv(path, headers, data);
            files{end + 1} = path;
        case 'massx'
            path = [basePath, '_mass_x.csv'];
            headers = {'coordinate'};
            data = result.coordinate(:);
            for i = 1:numel(result.densityVars)
                name = regexprep(result.densityVars{i}, '[^A-Za-z0-9_]', '_');
                if isfield(result, 'localParticleCount')
                    headers{end + 1} = ['particle_count_', name]; %#ok<AGROW>
                    data(:, end + 1) = result.localParticleCount(:, i); %#ok<AGROW>
                end
                headers{end + 1} = ['density_', name]; %#ok<AGROW>
                headers{end + 1} = ['cumulative_', name]; %#ok<AGROW>
                data(:, end + 1) = result.density(:, i); %#ok<AGROW>
                data(:, end + 1) = result.cumulativeDensity(:, i); %#ok<AGROW>
            end
            pd_write_numeric_csv(path, headers, data);
            files{end + 1} = path;
        case 'network2d'
            files = exportNetwork(result, basePath);
    end
end

function files = exportNetwork(result, basePath)
    files = {};
    fraction = double(result.poreMask);
    if isfield(result, 'cutCell') && isfield(result.cutCell, 'pore') && ...
            isfield(result.cutCell.pore, 'fraction')
        fraction = result.cutCell.pore.fraction;
    end
    [xGrid, yGrid] = meshgrid(result.xCenters, result.yCenters);
    path = [basePath, '_phase_grid.csv'];
    pd_write_numeric_csv(path, {'x','y','valid','poreFraction'}, ...
        [xGrid(:), yGrid(:), double(result.validMask(:)), fraction(:)]);
    files{end + 1} = path;

    phases = {'pore','matrix'};
    for i = 1:numel(phases)
        phaseName = phases{i};
        if ~isfield(result, phaseName) || ~isfield(result.(phaseName), 'components')
            continue;
        end
        comp = result.(phaseName).components;
        fields = {'area','equivDiameter','centroidX','centroidY', ...
            'perimeterOpen','interfacePerimeter'};
        n = numel(comp.area);
        data = (1:n).';
        headers = {'component'};
        for j = 1:numel(fields)
            if isfield(comp, fields{j})
                data(:, end + 1) = comp.(fields{j})(:); %#ok<AGROW>
                headers{end + 1} = fields{j}; %#ok<AGROW>
            end
        end
        path = [basePath, '_', phaseName, '_components.csv'];
        pd_write_numeric_csv(path, headers, data);
        files{end + 1} = path; %#ok<AGROW>
    end
end
