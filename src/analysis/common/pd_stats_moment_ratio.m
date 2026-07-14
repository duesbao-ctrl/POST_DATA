function value = pd_stats_moment_ratio(input, numeratorPower, denominatorPower)
%PD_STATS_MOMENT_RATIO Generalized moment-ratio mean used by all analyzers.

    input = input(:);
    input = input(isfinite(input));
    if isempty(input)
        value = NaN;
        return;
    end
    numerator = sum(input .^ numeratorPower);
    denominator = sum(input .^ denominatorPower);
    if ~isfinite(numerator) || ~isfinite(denominator) || denominator == 0
        value = NaN;
    else
        value = numerator / denominator;
    end
end
