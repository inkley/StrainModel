function result = hydrostatic_testing_v19(dataRoot, options)
% HYDROSTATIC_TESTING_V19 Rebuild manuscript Figures 5–6 using MATLAB only.
% Supply the external recording directory or set SENSOR_TRIAL_DATA_ROOT.

% For a fast redraw:
% HYDROSTATIC_TESTING_V19(Mode="render",OutputDir=existingDir) Raw recordings,
% workbook, and machine-specific paths stay outside Git.

arguments
    dataRoot (1,1) string = string(getenv('SENSOR_TRIAL_DATA_ROOT'))
    options.Mode (1,1) string {mustBeMember(options.Mode,["reduce","render","audit"])} = "reduce"
    options.OutputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')),'results','both_good')
    options.Workbook (1,1) string = ""
    options.Sheet = 1
end

result = struct();
if options.Mode ~= "render"
    assert(strlength(dataRoot) > 0,'Supply dataRoot or set SENSOR_TRIAL_DATA_ROOT.');
    result = runExperimentalReduction(dataRoot,OutputDir=options.OutputDir, ...
    Workbook = options.Workbook, Sheet = options.Sheet, AuditOnly = options.Mode == "audit");
end

if options.Mode == "audit", return; end
    verifyManuscriptSummary(options.OutputDir);
    renderExperimentalFigures(options.OutputDir);
end
