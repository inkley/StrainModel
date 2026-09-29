function result = runExperimentalReduction(dataRoot, options)
% RUNEXPERIMENTALREDUCTION MATLAB-only entry point; raw data stay external.
% Example: runExperimentalReduction('/absolute/path/1 Sensor Module v1')
% Set AuditOnly=true to inspect workbook selection without loading recordings.

arguments
    dataRoot (1,1) string
    options.OutputDir (1,1) string = ""
    options.Workbook (1,1) string = ""
    options.Selection (1,1) string = "both_good"
    options.Sheet = 1
    options.AuditOnly (1,1) logical = false
end
workbook = options.Workbook;
if strlength(options.OutputDir) == 0
    options.OutputDir = fullfile(fileparts(mfilename('fullpath')),'results',options.Selection);
end
if strlength(workbook) == 0, workbook = fullfile(dataRoot,'SM v1 Test Cases.xlsx'); end
result = reduceHydrostaticTrials(dataRoot,workbook,options.OutputDir, ...
    Selection=options.Selection,Sheet=options.Sheet,AuditOnly=options.AuditOnly);
end
