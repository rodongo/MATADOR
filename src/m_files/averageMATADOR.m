% Brief summary of this function.
% Detailed explanation of this function.
function simResults = averageMATADOR(dModel, files, diffFBS, VRefDis, saveLoc, sampleName)
% averageMATADOR is helper function for MATADOR for performing averaged analysis. 
% It evaluates metabolite consumption using MOMA and rMTA.
%
% This function performs anti-metabolite simulations by closing reactions 
% associated with specific metabolites and comparing metabolic scores.
%
% Usage:
%   simResults = averageMATADOR(dModel, files, diffFBS, VRefDis, saveLoc, sampleName)

    % --- Data Loading and Initialization ---
    cenData = readtable(files, 'FileType', 'spreadsheet');
    
    % Focus on consuming reactions and remove entries with missing data
    cenData = cenData(:, ["Metabolite", "AssociatedConsuming"]); 
    cenData = rmmissing(cenData);
    numMetabs = size(cenData, 1);

    % Initialize output structures/tables if data exists
    if numMetabs > 0
        % Preallocate arrays for performance
        scores = zeros(numMetabs, 2);
        runTimes = zeros(numMetabs, 2);
        rxnIdxList = cell(1, numMetabs);
        
        % Initialize struct array for detailed results
        metabDetails = struct('MetName', [], 'RunTime', []);

        % --- Parallel Processing of Metabolite Simulations ---
        parfor k = 1:numMetabs
            currentMetab = string(cenData.Metabolite(k));
            rawRxnString = string(cenData.AssociatedConsuming(k));
            
            % Clean and split reaction strings
            targetRxns = strtrim(strsplit(rawRxnString, ','));
            
            % Map reactions to model indices and filter out missing reactions
            validRxnLogic = ismember(dModel.rxns, targetRxns);
            rxnIdxList{k} = find(validRxnLogic);
            filteredRxns = dModel.rxns(validRxnLogic);

            % Timing and Simulation
            tSim = tic;
            
            % Simulation 1: MATADOR Comparison (MOMA-based)
            resMOMA = MATADOR(dModel, currentMetab, filteredRxns, diffFBS, VRefDis);
            
            % Simulation 2: rMTA Metabolite
            resRMTA = rMTA_Metabolite(dModel, currentMetab, filteredRxns, diffFBS, VRefDis, 'parameterK', 1);
            
            executionTime = toc(tSim);
            
            % Store results in thread-safe temporary variables
            runTimes(k, :) = [resMOMA.RunTimeMOMA, resRMTA.rMTAtimer];
            scores(k, :) = [resMOMA.aDATScoreMOMA, resRMTA.score_rMTA];
            
            fprintf('\tMetabolite [%s] solved in: %4.2f seconds\n', currentMetab, executionTime);
        end

        % --- Data Formatting and Export ---
        metabTS_1 = table(cenData.Metabolite, scores(:,1), scores(:,2), ...
            'VariableNames', {'Metabolite', 'aDATScore_MOMA', 'aDATScore_rMTA'});
            
        metabTS_Time = table(cenData.Metabolite, runTimes(:,1), runTimes(:,2), ...
            'VariableNames', {'Metabolite', 'Time_MOMA', 'Time_rMTA'});

        % Save to Excel
        outputFile = fullfile(saveLoc, [sampleName, '.xlsx']);
        writetable(metabTS_1, outputFile, 'Sheet', 'TScores');
        writetable(metabTS_Time, outputFile, 'Sheet', 'RunTime');

        % Populate Results Structure
        simResults.aDAT_ScoresTable = metabTS_1;
        simResults.RunTime = metabTS_Time;
        simResults.ConsumingRxnsIndices = rxnIdxList;
    else
        simResults = struct();
        warning('No valid metabolite data found in the provided file.');
    end
end