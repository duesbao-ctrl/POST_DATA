function [rows, columns] = pd_subplot_grid(count)
%PD_SUBPLOT_GRID Choose a compact near-square dashboard arrangement.

    if ~(isnumeric(count) && isscalar(count) && isfinite(count) && ...
            count >= 1 && count == round(count))
        error('postdata:BadViewCount', 'View count must be a positive integer.');
    end
    rows = max(1, floor(sqrt(count)));
    columns = ceil(count / rows);
end
