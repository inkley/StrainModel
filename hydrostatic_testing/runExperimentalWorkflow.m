function result = runExperimentalWorkflow(dataRoot, outputDir)
% RUNEXPERIMENTALWORKFLOW Reduce trials, check Table II, then draw Figs. 5–6.
% Stops if current summaries do not agree with the locked manuscript.
% For redraws alone, call renderExperimentalFigures(outputDir).

arguments
    dataRoot (1,1) string
    outputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')),'results','both_good')
end

result = runExperimentalReduction(dataRoot,OutputDir=outputDir);
verifyManuscriptSummary(outputDir);
renderExperimentalFigures(outputDir);
end
