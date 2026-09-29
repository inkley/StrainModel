function renderExperimentalFigures(reductionDir, outputDir)
% RENDEREXPERIMENTALFIGURES Draw Figures 5–6 from saved recording means.
% No raw recordings are opened. Both figures show uncorrected voltages and
% between-trial sample SD. Fits in Figure 6 are not Figure 9 trial slopes.

arguments
    reductionDir (1,1) string
    outputDir (1,1) string = fullfile(reductionDir,'figures')
end

pointsFile = fullfile(reductionDir,'calibration_points.csv');
p = readtable(pointsFile);
assert(all(ismember({'trial','scale','beta','voltage'},p.Properties.VariableNames)), ...
    'Missing calibration point columns.');
scales = [75 80 85 90 100];
assert(isequal(sort(unique(p.scale))',scales),'Expected all five interface conditions.');
for id = unique(p.trial)'
    rows = p(p.trial == id,:);
    assert(height(rows) == 25 && numel(unique(rows.scale)) == 1, ...
        'Trial %d must have 25 angles and a single scale.',id);
    assert(max(abs(sort(rows.beta)-(0:24)'*pi/12)) < 1e-10, ...
        'Trial %d has incomplete or duplicate angles.',id);
end
assert(all(isfinite(p.voltage)),'Nonfinite recording means.');
renderFigure05(pointsFile,outputDir);
renderFigure06(pointsFile,outputDir);
end
