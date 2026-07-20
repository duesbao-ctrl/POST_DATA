function value = pd_network_is_periodic_axis(boundary, axisName)
%PD_NETWORK_IS_PERIODIC_AXIS Test whether one grid axis is periodic.

    if strcmp(axisName, 'x')
        value = strcmp(boundary, 'periodic-x') || strcmp(boundary, 'periodic-xy');
    elseif strcmp(axisName, 'y')
        value = strcmp(boundary, 'periodic-y') || strcmp(boundary, 'periodic-xy');
    else
        error('pd_network_is_periodic_axis:InvalidAxis', ...
            'Axis must be ''x'' or ''y''; received ''%s''.', axisName);
    end
end
