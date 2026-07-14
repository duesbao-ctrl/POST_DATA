function views = pd_result_plot_views(result)
%PD_RESULT_PLOT_VIEWS Return available visualization views for a result.

    if ~isstruct(result) || ~isfield(result, 'analysisType')
        error('postdata:plotViews:BadResult', 'Result must contain analysisType.');
    end
    views = emptyViews();
    switch lower(strtrim(pd_to_char(result.analysisType)))
        case 'chunk'
            views = addView(views, 'field', '场分布');
            views = addView(views, 'histogram', '数值直方图');
            if isfield(result, 'y') && ~isempty(result.y)
                views = addView(views, 'profile-x', 'X 方向均值剖面');
                views = addView(views, 'profile-y', 'Y 方向均值剖面');
            end
        case 'cluster'
            views = addView(views, 'count', '粒径计数分布');
            views = addView(views, 'probability', '粒径概率分布');
            views = addView(views, 'cdf', '粒径累积分布');
            if isfield(result, 'meanByBin') && ~isempty(result.meanByBin.centers)
                views = addView(views, 'mean', '平均粒径位置分布');
            end
        case 'vx'
            views = addView(views, 'cumulative', '累积分布');
            views = addView(views, 'differential', '微分分布');
        case 'massx'
            views = addView(views, 'cumulative', '累计分布（保证单调递减）');
            views = addView(views, 'differential', '局部面密度（非累计）');
        case 'network2d'
            views = addView(views, 'phase', '孔隙相分布');
            if isfield(result, 'pore') && isfield(result.pore, 'labelGrid')
                views = addView(views, 'pore-label', '孔隙连通组分');
            end
            if isfield(result, 'matrix') && isfield(result.matrix, 'labelGrid')
                views = addView(views, 'matrix-label', '基体连通组分');
            end
            if hasComponentDiameter(result, 'pore')
                views = addView(views, 'pore-diameter', '孔隙等效直径');
            end
            if hasComponentDiameter(result, 'matrix')
                views = addView(views, 'matrix-diameter', '基体等效直径');
            end
            if hasProfile(result, 'x')
                views = addView(views, 'profile-x', 'X 方向孔隙率剖面');
            end
            if hasProfile(result, 'y')
                views = addView(views, 'profile-y', 'Y 方向孔隙率剖面');
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
