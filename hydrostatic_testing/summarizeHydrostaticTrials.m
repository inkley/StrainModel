function summary = summarizeHydrostaticTrials(diagnostics)
% SUMMARIZEHYDROSTATICTRIALS Group dimensionless normalized trial slopes.
% SEM is sample SD/sqrt(n), excluding missing same-day references. It does
% not propagate shared-reference uncertainty. Groups with n <= 1 have NaN
% SEM.
scales  = [75; 80; 85; 90; 100]; n = numel(scales);
summary = table(scales,zeros(n,1),zeros(n,1),nan(n,1),nan(n,1), ...
    'VariableNames',{'scale','raw_n','normalized_n','mean','sem'});

for k = 1:n
    rows = diagnostics(diagnostics.scale == scales(k),:);
    values = rows.normalized_slope(isfinite(rows.normalized_slope));
    summary.raw_n(k) = height(rows); summary.normalized_n(k) = numel(values);
    if ~isempty(values), summary.mean(k) = mean(values); end
    if numel(values) > 1, summary.sem(k) = std(values,0)/sqrt(numel(values)); end
end
end
