function stats = pd_stats_summary(input, numeratorPower, denominatorPower)
%PD_STATS_SUMMARY Shared descriptive statistics with generalized mean.

    if nargin < 2, numeratorPower = 1; end
    if nargin < 3, denominatorPower = 0; end
    input = input(:);
    input = input(isfinite(input));
    stats = struct();
    stats.n = numel(input);
    stats.quantileProb = [0.05, 0.25, 0.50, 0.75, 0.95];
    if isempty(input)
        stats.min = NaN;
        stats.max = NaN;
        stats.mean = NaN;
        stats.std = NaN;
        stats.median = NaN;
        stats.quantileValue = nan(size(stats.quantileProb));
        return;
    end
    stats.min = min(input);
    stats.max = max(input);
    stats.mean = pd_stats_moment_ratio(input, numeratorPower, denominatorPower);
    stats.std = std(input);
    stats.median = median(input);
    stats.quantileValue = pd_stats_quantiles(input, stats.quantileProb);
end
