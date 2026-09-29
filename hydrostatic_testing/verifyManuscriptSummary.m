function comparison = verifyManuscriptSummary(reductionDir)
% VERIFYMANUSCRIPTSUMMARY Check Table II counts and three-decimal reporting.
% Baseline is the September 28 final v13 manuscript, not a new fitted target.

s = readtable(fullfile(reductionDir,'experimental_summary.csv'));
expected = [75 10 7 1.225 .042;80 6 3 1.270 .045;85 6 6 1.327 .082;90 11 10 1.152 .085];
assert(height(s) == 5 && numel(unique(s.scale)) == 5,'Expected five unique groups.');
comparison = array2table(expected,'VariableNames',{'scale','raw_n','normalized_n','mean','sem'});

for k = 1:4
    row = s(s.scale == expected(k,1),:);
    assert(height(row) == 1,'Missing membrane group.');
    assert(isequal([row.raw_n row.normalized_n],expected(k,2:3)), 'Manuscript trial-count mismatch.');
    assert(all(abs([row.mean row.sem]-expected(k,4:5)) < .0005), 'Manuscript mean/SEM rounding mismatch.');
end

bare = s(s.scale == 100,:);
assert(height(bare) == 1 && bare.raw_n == 17 && abs(bare.mean-1) < 1e-10,'Bare-port reference mismatch.');
fprintf('Table II counts and reported means/SEM match the v13 manuscript.\n');
end
