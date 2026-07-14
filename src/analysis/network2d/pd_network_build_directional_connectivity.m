function connectivity = pd_network_build_directional_connectivity(comp, boundary, phaseArea)
    connectivity = struct();
    connectivity.x = buildConnectivityAxis(comp, boundary, phaseArea, 'x');
    connectivity.y = buildConnectivityAxis(comp, boundary, phaseArea, 'y');
end

function info = buildConnectivityAxis(comp, boundary, phaseArea, axisName)
    if strcmp(axisName, 'x')
        if pd_network_is_periodic_axis(boundary, 'x')
            componentMask = comp.wrapX;
            criterion = 'wraps-periodic-boundary';
        else
            componentMask = comp.percolatesX;
            criterion = 'touches-both-open-boundaries';
        end
    else
        if pd_network_is_periodic_axis(boundary, 'y')
            componentMask = comp.wrapY;
            criterion = 'wraps-periodic-boundary';
        else
            componentMask = comp.percolatesY;
            criterion = 'touches-both-open-boundaries';
        end
    end

    ids = comp.label(componentMask);
    areas = comp.area(componentMask);

    info = struct();
    info.axis = axisName;
    info.criterion = criterion;
    info.isConnected = any(componentMask);
    info.componentIds = ids(:).';
    info.componentCount = numel(ids);
    info.totalConnectedArea = sum(areas);
    info.connectedAreaFraction = safeDivideOrZero(sum(areas), phaseArea);
    if isempty(areas)
        info.largestConnectedArea = 0;
    else
        info.largestConnectedArea = max(areas);
    end
end

function y = safeDivideOrZero(a, b)
    if isempty(b) || ~isfinite(b) || b == 0
        if isempty(a) || ~isfinite(a)
            y = NaN;
        elseif a == 0
            y = 0;
        else
            y = NaN;
        end
    else
        y = a ./ b;
    end
end
