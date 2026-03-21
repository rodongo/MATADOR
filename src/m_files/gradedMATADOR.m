function simResults = gradedMATADOR(dModel, files, diffFBS, VRefDis, gradeLevels, saveLoc, sampleName)
% gradedMATADOR  Run graded MATADOR sensitivity analysis for anti-metabolite simulations.
%
%   simResults = gradedMATADOR(dModel, files, diffFBS, VRefDis, ...
%                               gradeLevels, saveLoc, sampleName)
%
%   Performs graded metabolic disruption analysis using the MATADOR framework.
%   For each metabolite listed in the input spreadsheet, the function closes
%   all associated consuming reactions at each specified grade level and calls
%   the core scoring routine (mainGradedMATADOR). Results are written to an
%   Excel workbook and returned as a structured output.
%
%   The underlying simulation models an anti-metabolite scenario: all reactions
%   that consume the target metabolite are suppressed, and the metabolite's
%   uptake from the environment is blocked, regardless of its availability.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   dModel      - (struct)  COBRA-compatible metabolic model. Must contain
%                           the field 'rxns' (cell array of reaction IDs).
%
%   files       - (char)    Full file path to the input spreadsheet (.xlsx).
%                           The sheet must contain at minimum two columns:
%                             'Metabolite'          – metabolite name
%                             'AssociatedConsuming'  – comma-separated list
%                                                      of reaction IDs that
%                                                      consume the metabolite
%
%   diffFBS     -   Forward - Backward - Unchanged (+1;0;-1) for each reaction
%                           based on differential gene expression and
%                           reference flux distribution. The output from
%                           diffexprs2rxnFBS function in cobratoolbox
%                           
%
%   VRefDis     - (numeric)  Reference flux distribution used as the baseline
%                           for MATADOR score computation.
%
%   gradeLevels - (numeric) Row or column vector of grade levels (e.g.,
%                           [0, 25, 50, 75]) representing the degree of
%                           anti-metabolite inhibition to simulate (%).
%
%   saveLoc     - (char)    Directory path where the output Excel file will
%                           be written. Must end with a file separator (e.g.,
%                           '/results/').
%
%   sampleName  - (char)    Identifier string prepended to the output filename.
%                           Output file: <saveLoc><sampleName>_SensitivityAnalysis.xlsx
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   simResults  - (struct)  Structure with the following field:
%       .MATADOR_ScoresTable  - (table) Metabolite-by-grade MATADOR scores.
%                               Columns: Metabolite, MATADOR_0, MATADOR_25,
%                               MATADOR_50, MATADOR_75 (or as many columns
%                               as entries in gradeLevels).
%                               Returns empty struct if no valid metabolites
%                               are found in the input file.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   - Reactions listed in 'AssociatedConsuming' that are absent from dModel
%     are silently excluded; only model-present reactions are used.
%   - The sentinel value 65535 may appear in upstream data and should be
%     treated as missing/overflow in downstream analyses.
%   - Timing for each metabolite-grade call is measured via tic/toc and
%     stored in 'time_aDAT' (currently not aggregated; extend as needed).
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
%   gradeLevels = [0, 25, 50, 75];
%   results = gradedMATADOR(myModel, 'metabolites.xlsx', diffFBS, VRef, ...
%                            gradeLevels, './output/', 'SampleA');
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   mainGradedMATADOR  – core MATADOR scoring function (must be on MATLAB path)
%
% -------------------------------------------------------------------------
% AUTHORS / VERSION
% -------------------------------------------------------------------------
%   MATADOR Toolbox  |  See accompanying documentation for citation details.
% =========================================================================

    % ------------------------------------------------------------------
    % 1. LOAD METABOLITE–REACTION MAPPING
    % ------------------------------------------------------------------
    % Read the input spreadsheet and retain only the two columns required
    % for consuming-reaction simulations. Rows with missing values in
    % either column are discarded.
    cenData = readtable(files, 'FileType', 'spreadsheet');
    cenData = cenData(:, ["Metabolite", "AssociatedConsuming"]);
    cenData = rmmissing(cenData);   % remove rows with any missing values

    % ------------------------------------------------------------------
    % 2. RUN GRADED SIMULATIONS (if valid metabolite entries exist)
    % ------------------------------------------------------------------
    if size(cenData, 1) > 0

        % Pre-allocate score matrix: rows = metabolites, cols = grade levels
        metabTS  = zeros(size(cenData, 1), numel(gradeLevels));
        rxnIdx1  = cell(1, size(cenData, 1));   % stores reaction indices per metabolite

        % ---- Outer loop: iterate over metabolites --------------------
        for k = 1:size(cenData, 1)

            % ---- Inner loop: iterate over grade levels ---------------
            for l = 1:numel(gradeLevels)

                % -- Parse reaction list for metabolite k --------------
                % Reactions are stored as a single comma-separated string;
                % split into a cell array and strip any whitespace.
                metabK             = struct();
                metabK.MetabName   = cenData.Metabolite(k);
                metabK.Reaction    = string(cenData.AssociatedConsuming(k));
                metabK.Reaction    = cellstr(strsplit(metabK.Reaction, ','));
                metabK.Reaction    = metabK.Reaction';
                metabK.Reaction    = strrep(metabK.Reaction, ' ', '');

                % -- Map reaction IDs to model indices ------------------
                % Record indices regardless of model membership (for reference).
                rxnIdx1{k} = find(ismember(dModel.rxns, metabK.Reaction));

                % Restrict to reactions that actually exist in the model;
                % reactions absent from dModel are silently skipped.
                metabK.Reaction = dModel.rxns(ismember(dModel.rxns, metabK.Reaction));

                % -- Call core MATADOR scoring routine ------------------
                timerVal  = tic;
                metabKTS  = mainGradedMATADOR(dModel, ...
                                               string(metabK.MetabName), ...
                                               metabK.Reaction, ...
                                               diffFBS, ...
                                               gradeLevels(l), ...
                                               VRefDis);
                time_aDAT = toc(timerVal);   %#ok<NASGU> timing for external logging if needed

                % -- Store MATADOR score for this metabolite × grade ---
                metabTS(k, l) = metabKTS.MATADOR_Score;

            end % grade level loop
        end % metabolite loop

        % ------------------------------------------------------------------
        % 3. ASSEMBLE RESULTS TABLE
        % ------------------------------------------------------------------
        % Package scores into a labelled table with one row per metabolite
        % and one column per grade level.
        metabTS_1 = table(cenData.Metabolite, ...
                          metabTS(:,1), metabTS(:,2), ...
                          metabTS(:,3), metabTS(:,4), ...
                          metabTS(:,5));
        metabTS_1.Properties.VariableNames = ...
            {'Metabolite', 'MATADOR_0', 'MATADOR_25', 'MATADOR_50', 'MATADOR_75', 'MATADOR_100'};

    else
        % No valid rows found in input file; return empty placeholder
        metabTS_1 = struct();
    end

    % ------------------------------------------------------------------
    % 4. WRITE OUTPUT TO EXCEL
    % ------------------------------------------------------------------
    % Save the scores table to the 'TScores' sheet of the output workbook,
    % provided that a non-empty result was produced.
    if ~isempty(metabTS_1)
        outputPath = append(saveLoc, sampleName, '_SensitivityAnalysis.xlsx');
        writetable(metabTS_1, outputPath, ...
                   'FileType', 'spreadsheet', ...
                   'Sheet',    'TScores');
    end

    % ------------------------------------------------------------------
    % 5. PACKAGE AND RETURN RESULTS
    % ------------------------------------------------------------------
    simResults.MATADOR_ScoresTable = metabTS_1;

end % gradedMATADOR