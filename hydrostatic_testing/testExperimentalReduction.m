function testExperimentalReduction
% Synthetic ingestion plus archived-summary regression; no cloud data read.
root        = fileparts(mfilename('fullpath')); addpath(root);
fixture     = tempname; mkdir(fixture);
workbook    = fullfile(fixture,'register.xlsx');

rows = {'trial','date','scale','fluid','resolution','reliability','threshold','note'; ...
    1,'2025-01-01','None','Air',1,'Good','Good',''; ...
    2,'2025-01-01',.85,'Air',1,'Good','Good',''; ...
    3,'2025-01-02',.85,'Air',1,'Good','Good',''; ...
    4,'2025-01-02','None','Air',1,'Good','Bad (A)',''; ...
    5,'2025-01-01',.85,'Mineral Oil',1,'Good','Good',''; ...
    6,'2025-01-01',.85,'Air',1,'Bad','Good',''};

writecell(rows,workbook);

audit = reduceHydrostaticTrials(fixture,workbook,fullfile(fixture,'audit'),AuditOnly=true);
assert(isequal(audit.selection.trial(audit.selection.selected),[1;2;3]));
pressure = 1000*9.81*.040004*sin((0:24)'*pi/12);
slopes = [.001 .0013 .0015 .001]; noise = repmat([-.01 .01],1,5000);

for i = 1:4
    folder = fullfile(fixture,sprintf('Trial %d',i),'Data'); mkdir(folder);
    for j = 1:25
        angle = (j-1)*15;
        if i == 2 && angle == 195, angle = 205; end
        raw = zeros(3,10000); raw(2,:) = .5+slopes(i)*pressure(j)+noise;
        writematrix(raw,fullfile(folder,sprintf('Record %03d.txt',angle)),'Delimiter','tab');
    end
end

r = reduceHydrostaticTrials(fixture,workbook,fullfile(fixture,'out'));
assert(max(abs(r.diagnostics.abs_slope-slopes(1:3)')) < 1e-12);
assert(abs(r.diagnostics.normalized_slope(2)-1.3) < 1e-10);
assert(isnan(r.diagnostics.normalized_slope(3)));
assert(sum(~r.manifest.filename_angle_matches) == 1);
assert(max(abs(r.points.within_record_sd-std(noise,0))) < 1e-12);
assert(all(r.diagnostics.recomputed_threshold_pass));
r2 = reduceHydrostaticTrials(fixture,workbook,fullfile(fixture,'alternate'),Selection='longitudinal_good');
assert(abs(r2.diagnostics.normalized_slope(3)-1.5) < 1e-10);
assert(r2.summary.normalized_n(r2.summary.scale == 85) == 2);
assert(abs(r2.summary.sem(r2.summary.scale == 85)-.1) < 1e-10);

% Reject incomplete recordings rather than silently altering sample counts.
movefile(fullfile(fixture,'Trial 1','Data','Record 000.txt'),fullfile(fixture,'saved.txt'));
failed = false;
try
    reduceHydrostaticTrials(fixture,workbook,fullfile(fixture,'incomplete'));
catch err
    failed = contains(err.message,'expected 25 files');
end
assert(failed);

% Recompute archived slopes from exported recording means (not raw samples).
archive = fullfile(root,'..','Archive','calibration_rebuild','both_good');
if isfolder(archive)
    d = readtable(fullfile(archive,'trial_diagnostics.csv'),'TextType','string');
    p = readtable(fullfile(archive,'calibration_points.csv'));
    expected = readtable(fullfile(archive,'experimental_summary.csv'));
    original = d.abs_slope;
    for i = 1:height(d)
        pt = p(p.trial == d.trial(i),:);
        fit = polyfit(pt.pressure,pt.voltage,1); d.abs_slope(i) = abs(fit(1));
    end
    assert(max(abs(original-d.abs_slope)) < 1e-12);
    for i = 1:height(d)
        bare = d.abs_slope(d.scale == 100 & d.date == d.date(i));
        d.normalized_slope(i) = d.abs_slope(i)/mean(bare);
    end
    actual = summarizeHydrostaticTrials(d);
    assert(isequal(actual{:,1:3},expected{:,1:3}));
    assert(max(abs(actual{:,4:5}-expected{:,4:5}),[],'all') < 1e-11);
    fprintf('Archived 50-trial summary parity passed.\n');
else
    fprintf('SKIPPED archived-summary comparison: local Archive directory unavailable.\n');
end
fprintf('All native MATLAB reduction tests passed. Fixtures: %s\n',fixture);
end
