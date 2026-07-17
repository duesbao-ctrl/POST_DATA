function views = pd_result_plot_views(result)
%PD_RESULT_PLOT_VIEWS Return available visualization views for a result.

    if ~isstruct(result) || ~isfield(result, 'analysisType')
        error('postdata:plotViews:BadResult', 'Result must contain analysisType.');
    end
    views = emptyViews();
    switch lower(strtrim(pd_to_char(result.analysisType)))
        case 'chunk'
            views = addView(views, 'field', pd_ui_text('Field distribution'));
            views = addView(views, 'histogram', pd_ui_text('Value histogram'));
            if isfield(result, 'dimension') && strcmpi(result.dimension, '2d')
                views = addView(views, 'profile-x', pd_ui_text('Mean profile in X direction'));
                views = addView(views, 'profile-y', pd_ui_text('Mean profile in Y direction'));
            end
        case 'cluster'
            views = addView(views, 'count', pd_ui_text('Particle-size count distribution'));
            views = addView(views, 'probability', pd_ui_text('Particle-size probability distribution'));
            views = addView(views, 'cdf', pd_ui_text('Particle-size cumulative distribution'));
            if isfield(result, 'meanByBin') && ~isempty(result.meanByBin.centers)
                views = addView(views, 'mean', pd_ui_text('Mean particle size by position'));
            end
        case 'vx'
            views = addView(views, 'cumulative', pd_ui_text('Cumulative distribution'));
            views = addView(views, 'differential', pd_ui_text('Differential distribution'));
        case 'massx'
            views = addView(views, 'cumulative', pd_ui_text('Cumulative distribution (monotonic decreasing)'));
            views = addView(views, 'differential', pd_ui_text('Local areal density (non-cumulative)'));
            if isfield(result, 'localParticleCount')
                views = addView(views, 'particle-count', ...
                    pd_ui_text('SPH particle count by X bin'));
            end
        case 'network2d'
            views = addView(views, 'phase', pd_ui_text('Pore phase distribution'));
            if isfield(result, 'pore') && isfield(result.pore, 'labelGrid')
                views = addView(views, 'pore-label', pd_ui_text('Pore connected components'));
            end
            if isfield(result, 'matrix') && isfield(result.matrix, 'labelGrid')
                views = addView(views, 'matrix-label', pd_ui_text('Matrix connected components'));
            end
            if hasComponentDiameter(result, 'pore')
                views = addView(views, 'pore-diameter', pd_ui_text('Pore equivalent diameter'));
            end
            if hasComponentDiameter(result, 'matrix')
                views = addView(views, 'matrix-diameter', pd_ui_text('Matrix equivalent diameter'));
            end
            if hasProfile(result, 'x')
                views = addView(views, 'profile-x', pd_ui_text('Porosity profile in X direction'));
            end
            if hasProfile(result, 'y')
                views = addView(views, 'profile-y', pd_ui_text('Porosity profile in Y direction'));
            end
        otherwise
            error('postdata:plotViews:UnknownType', 'Unsupported analysis type.');
    end
end

function tf = hasComponentDiameter(result, phase)
    tf = isfield(result, phase) && isfield(result.(phase), 'components') && ...
        isfield(result.(phase).components, 'equivDiameter') && ...
        ~isempty(result.(phase).components.equivDiameter);
end

function tf = hasProfile(result, axisName)
    tf = isfield(result, 'profile') && isfield(result.profile, axisName) && ...
        isfield(result.profile.(axisName), 'centers') && ...
        ~isempty(result.profile.(axisName).centers);
end

function views = emptyViews()
    views = struct('Id', {}, 'Label', {});
end

function views = addView(views, id, label)
    views(end + 1) = struct('Id', id, 'Label', label); %#ok<AGROW>
end
