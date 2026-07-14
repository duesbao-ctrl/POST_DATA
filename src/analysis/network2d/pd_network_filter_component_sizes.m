function [diameterOut, areaOut, mask] = ...
        pd_network_filter_component_sizes(diameter, area, diameterRange)
%PD_NETWORK_FILTER_COMPONENT_SIZES Apply the equivalent-diameter range.

    diameter = diameter(:);
    area = area(:);
    mask = isfinite(diameter) & (diameter > 0);
    if ~isempty(diameterRange)
        lo = min(diameterRange(:));
        hi = max(diameterRange(:));
        mask = mask & (diameter >= lo) & (diameter <= hi);
    end
    diameterOut = diameter(mask);
    areaOut = area(mask);
end
