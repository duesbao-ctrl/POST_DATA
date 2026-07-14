function values = pd_stats_quantiles(input, probabilities)
%PD_STATS_QUANTILES Linearly interpolated quantiles without toolbox dependency.

    input = sort(input(isfinite(input)), 'ascend');
    count = numel(input);
    values = nan(size(probabilities));
    if count == 0, return; end
    if count == 1
        values(:) = input(1);
        return;
    end
    for i = 1:numel(probabilities)
        probability = min(max(probabilities(i), 0), 1);
        position = 1 + (count - 1) * probability;
        lowerIndex = floor(position);
        upperIndex = ceil(position);
        if lowerIndex == upperIndex
            values(i) = input(lowerIndex);
        else
            weight = position - lowerIndex;
            values(i) = input(lowerIndex) * (1 - weight) + input(upperIndex) * weight;
        end
    end
end
