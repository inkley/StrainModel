function result = reduceHydrostaticTrials(dataRoot, workbook, outputDir, options)
% REDUCEHYDROSTATICTRIALS Reduce hydrostatic recordings using MATLAB only.
% dataRoot contains Trial <ID>/Data; workbook has the eight archived columns
% (trial, date, scale, fluid, resolution, reliability, threshold, note).
% Select the worksheet explicitly if the trial register is not the first.
% Records are read only. Generated CSV files overwrite matching output names.
% No recordings are copied to outputDir. AuditOnly reads the workbook only.

arguments
    dataRoot (1,1) string
    workbook (1,1) string
    outputDir (1,1) string
    options.Selection (1,1) string {mustBeMember(options.Selection,["both_good","longitudinal_good"])} = "both_good"
    options.Sheet = 1
    options.AuditOnly (1,1) logical = false
end

assert(isfile(workbook),'Trial register not found: %s',workbook);
assert(~strcmp(char(java.io.File(dataRoot).getCanonicalPath()),char(java.io.File(outputDir).getCanonicalPath())), ...
    'Use a separate output directory, not the raw-data directory.');
cells = readcell(workbook,'Sheet',options.Sheet);
assert(size(cells,2) >= 8,'Expected eight trial-register columns.');
records = cell(0,8);
for k = 2:size(cells,1)
    row = cells(k,1:8); id = row{1};
    if ~(isnumeric(id) && isscalar(id) && isfinite(id)), continue; end
    assert(id == fix(id),'Trial IDs must be integers.');
    scale = row{3}; pct = -1;
    if isnumeric(scale) && isscalar(scale) && isfinite(scale)
        pct = round(100*scale);
    elseif textValue(scale) == "None"
        pct = 100;
    elseif startsWith(textValue(scale),"75%")
        pct = 75;
    end
    fluid = textValue(row{4}); reliability = textValue(row{6}); threshold = textValue(row{7});
    thresholdOK = threshold == "Good" || (options.Selection == "longitudinal_good" && threshold == "Bad (A)");
    selected = reliability == "Good" && thresholdOK && fluid ~= "Mineral Oil" && ismember(pct,[75 80 85 90 100]);
    records(end+1,:) = {id,dateValue(row{2}),pct,fluid,reliability,threshold,textValue(row{8}),selected}; %#ok<AGROW>
end

selection = cell2table(records,'VariableNames',{'trial','date','scale','fluid','subjective','threshold','note','selected'});
assert(~isempty(selection),'No numeric trial IDs found in worksheet.');
assert(numel(unique(selection.trial)) == height(selection),'Duplicate trial IDs in workbook.');
if ~isfolder(outputDir), mkdir(outputDir); end
    writetable(selection,fullfile(outputDir,'trial_selection.csv'));
    chosen = selection(selection.selected,:);
    assert(~isempty(chosen),'No trials meet the requested selection.');
    assert(all(strlength(chosen.date) > 0),'Selected trials require valid dates for same-day normalization.');
    result = struct('selection',selection);
    fprintf('Selection %s: %d trials.\n',options.Selection,height(chosen));
if options.AuditOnly, return; end
    % Angles follow lexicographic recording order, including both 0 and 360 deg.
    % Filename angle tokens are audited but never used to silently remap data.
    beta = (0:24)'*pi/12;
    pressure = 1000*9.81*0.040004*sin(beta); % Pa, longitudinal port separation
    n = height(chosen); means = zeros(n,25); recordSD = means;
    diagnostics = chosen;
    diagnostics.slope = zeros(n,1); diagnostics.abs_slope = zeros(n,1);
    diagnostics.offset = zeros(n,1); diagnostics.rmse = zeros(n,1);
    diagnostics.return_difference = zeros(n,1);
    diagnostics.recomputed_threshold_pass = false(n,1);
    fileRows = cell(n*25,6);
for i = 1:n
    folder = fullfile(dataRoot,sprintf('Trial %d',chosen.trial(i)),'Data');
    files = dir(fullfile(folder,'*'));
    files = files(~[files.isdir] & ~startsWith({files.name},'.'));
    [~,order] = sort(string({files.name})); files = files(order);
    assert(numel(files) == 25,'Trial %d: expected 25 files; found %d. Check cloud availability.',chosen.trial(i),numel(files));
    for j = 1:25
        file = fullfile(files(j).folder,files(j).name);
        raw = load(file,'-ascii');
        assert(ismatrix(raw) && size(raw,1) >= 3 && size(raw,2) >= 10000, ...
            'Unexpected recording shape: %s',file);
        voltage = raw(2,1:10000);
        assert(all(isfinite(voltage)),'Nonfinite voltage sample: %s',file);
        means(i,j) = mean(voltage); recordSD(i,j) = std(voltage,0);
        token = regexp(files(j).name,'\s(\d{3})(?:\.[^.]*)?$','tokens','once');
        filenameAngle = -1;
        if ~isempty(token), filenameAngle = str2double(token{1}); end
        angle = (j-1)*15;
        fileRows((i-1)*25+j,:) = {chosen.trial(i),angle,string(files(j).name),filenameAngle,filenameAngle == angle,files(j).bytes};
    end
    coeff = polyfit(pressure,means(i,:)',1);
    diagnostics.slope(i) = coeff(1); diagnostics.abs_slope(i) = abs(coeff(1));
    diagnostics.offset(i) = coeff(2);
    diagnostics.rmse(i) = sqrt(mean((means(i,:)'-polyval(coeff,pressure)).^2));
    diagnostics.return_difference(i) = means(i,end)-means(i,1);
    diagnostics.recomputed_threshold_pass(i) = diagnostics.rmse(i) <= .08;
    fprintf('Reduced trial %d (%d/%d).\n',chosen.trial(i),i,n);
end

% Selection follows the archived flags; recomputed failures are reported,
% not silently removed. Figure-generation code must review these diagnostics.
diagnostics.same_day_bare_n = zeros(n,1);
diagnostics.normalized_slope = nan(n,1);
for i = 1:n
    bare = diagnostics.abs_slope(diagnostics.scale == 100 & diagnostics.date == diagnostics.date(i));
    diagnostics.same_day_bare_n(i) = numel(bare);
    if ~isempty(bare)
        assert(mean(bare) > 0,'Zero same-day bare-port reference.');
        diagnostics.normalized_slope(i) = diagnostics.abs_slope(i)/mean(bare);
    end
end

summary = summarizeHydrostaticTrials(diagnostics);
manifest = cell2table(fileRows,'VariableNames',{'trial','angle_deg','filename','filename_angle','filename_angle_matches','bytes'});
% Trial-major order preserves the legacy CSV layout. Within-record SD is
% different from the between-trial SD used to draw Figures 5 and 6.
points = table(repelem(chosen.trial,25),repelem(chosen.scale,25),repmat(beta,n,1), ...
    repmat(pressure,n,1),reshape(means',[],1),reshape(recordSD',[],1), ...
    reshape((means-diagnostics.offset)',[],1),'VariableNames', ...
    {'trial','scale','beta','pressure','voltage','within_record_sd','offset_corrected_voltage'});
writetable(diagnostics,fullfile(outputDir,'trial_diagnostics.csv'));
writetable(manifest,fullfile(outputDir,'file_angle_manifest.csv'));
writetable(summary,fullfile(outputDir,'experimental_summary.csv'));
writetable(points,fullfile(outputDir,'calibration_points.csv'));
result.diagnostics = diagnostics; result.manifest = manifest;
result.summary = summary; result.points = points;
disp(summary);
fprintf('Recomputed RMSE failures: %s\n',mat2str(diagnostics.trial(~diagnostics.recomputed_threshold_pass)'));
fprintf('No same-day reference: %s\n',mat2str(diagnostics.trial(diagnostics.same_day_bare_n == 0)'));
disp(manifest(~manifest.filename_angle_matches,:));
end

function value = textValue(input)
% TEXTVALUE Normalize blank spreadsheet cells without changing flag spelling.
if isempty(input) || (isscalar(input) && ismissing(input))
    value = "";
else
    value = string(input);
end
end

function value = dateValue(input)
% DATEVALUE Canonical date for grouping. Invalid dates are never pooled.
value = "";
if isdatetime(input) && ~isnat(input)
    input.Format = 'yyyy-MM-dd'; value = string(input);
elseif isnumeric(input) && isscalar(input) && isfinite(input)
    date = datetime(input,'ConvertFrom','excel','Format','yyyy-MM-dd'); value = string(date);
elseif (ischar(input) || isstring(input)) && strlength(string(input)) > 0
    try
        date = datetime(input,'Format','yyyy-MM-dd');
        if ~isnat(date), value = string(date); end
    catch
        error('Unrecognized workbook date: %s. Use Excel date cells.',string(input));
    end
end
end
