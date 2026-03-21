% Brief summary of this function.
% Detailed explanation of this function.
function personalizedModels = generatePersonalizedModels(model, expressionData)
% generatePersonalizedModels  Build context-specific metabolic models from
%                              omics expression data for healthy and disease cohorts.
%
%   personalizedModels = generatePersonalizedModels(model, expressionData)
%
%   For each omics dataset supplied in expressionData, this function generates
%   sample-level context-specific metabolic models for both a healthy/control
%   group and a disease group using the iMAT algorithm (iMATCustom). Gene
%   expression values are first mapped to reaction activity scores, and
%   dataset-specific quantile thresholds are derived from the pooled cohort
%   to ensure comparable model extraction across conditions.
%
%   The resulting personalised models encode which reactions are active in
%   each individual sample, enabling downstream flux simulations (e.g.,
%   MATADOR) to be grounded in patient- or sample-specific network topology.
%
%   Parallelisation: healthy and disease cohorts are each processed using
%   parfor; a parallel pool should be open before calling this function to
%   take full advantage of multi-core hardware.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   model           - (struct)  COBRA-compatible genome-scale metabolic model.
%                               Must contain at minimum:
%                                 .rxns  – cell array of reaction identifiers
%                                 .genes – cell array of gene identifiers
%                                          (required by mapExpressionToReactions)
%
%   expressionData  - (struct array)  Array of omics dataset structures, one
%                               element per dataset/cohort comparison.
%                               Each element must contain:
%                                 .OmicData  – dataset metadata or identifier
%                                              (passed through to output unchanged)
%                                 .Genes     – (G × 1 cell) gene identifiers
%                                              matching rows of .Control and .Disease
%                                 .Control   – (G × Nc double) expression matrix
%                                              for the healthy/control group;
%                                              G genes × Nc samples
%                                 .Disease   – (G × Nd double) expression matrix
%                                              for the disease group;
%                                              G genes × Nd samples
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   personalizedModels - (struct array)  One element per expressionData entry.
%                        Each element contains:
%       .OmicData    – dataset metadata (passed through from input)
%       .redModelsH  – (1 × Nc cell)    context-specific COBRA model for each
%                                        healthy/control sample, with unused
%                                        genes removed
%       .healthyMod  – (R × Nc logical) binary reaction-presence matrix for
%                                        healthy models; R = numel(model.rxns),
%                                        entry (r,j) = 1 if reaction r is active
%                                        in healthy sample j
%       .redModelsD  – (1 × Nd cell)    context-specific COBRA model for each
%                                        disease sample, with unused genes removed
%       .diseaseMod  – (R × Nd logical) binary reaction-presence matrix for
%                                        disease models (same convention as
%                                        .healthyMod)
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW
% -------------------------------------------------------------------------
%   For each dataset i:
%     1. Pool all samples (control + disease) and compute the 25th and 75th
%        percentile expression thresholds (th_low, th_high) across the
%        entire combined matrix. These thresholds define the iMAT low/high
%        expression categories in a dataset-relative manner.
%     2. For each control sample j (parallel):
%          a. Construct a gene-expression struct for sample j.
%          b. Map gene expression to reaction scores (mapExpressionToReactions).
%          c. Extract a context-specific model (iMATCustom) using th_low, th_high.
%          d. Remove genes no longer present in the extracted model.
%          e. Record reaction membership as a binary vector.
%     3. Repeat Step 2 for each disease sample j (parallel).
%     4. Aggregate reduced models and binary matrices into the output struct.
%
% -------------------------------------------------------------------------
% THRESHOLDING RATIONALE
% -------------------------------------------------------------------------
%   Thresholds are computed from the pooled (control + disease) expression
%   matrix rather than each group independently. This ensures that the iMAT
%   high/low classification is consistent across conditions, preventing
%   artificial differences in model topology caused by group-specific
%   expression scales.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   - parfor loops require the Parallel Computing Toolbox. If unavailable,
%     replace both parfor keywords with for without any other code changes.
%   - iMATCustom is expected to return a COBRA model struct whose .rxns field
%     is a subset of model.rxns. The binary membership vector is derived via
%     ismember, so reaction ordering in model.rxns is preserved.
%   - removeUnusedGenes trims the .genes field to only those genes with at
%     least one associated active reaction, reducing memory footprint.
%   - Expression values in .Control and .Disease should be pre-normalised
%     and on a consistent scale (e.g., TPM, FPKM, or log2-transformed counts)
%     before being passed to this function.
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   mapExpressionToReactions  – COBRA Toolbox utility
%   iMATCustom                – custom iMAT implementation (must be on path)
%   removeUnusedGenes         – COBRA Toolbox utility
%
% -------------------------------------------------------------------------
% AUTHORS / VERSION
% -------------------------------------------------------------------------
%   MATADOR Toolbox  |  See accompanying documentation for citation details.
% =========================================================================

    personalizedModels = struct();

    % ------------------------------------------------------------------
    % OUTER LOOP: iterate over omics datasets
    % ------------------------------------------------------------------
    for i = 1:numel(expressionData)

        % Pass dataset metadata through to output unchanged
        personalizedModels(i).OmicData = expressionData(i).OmicData;

        tData    = expressionData(i);

        % --------------------------------------------------------------
        % 1. COMPUTE POOLED EXPRESSION THRESHOLDS
        % --------------------------------------------------------------
        % Concatenate control and disease matrices column-wise to form the
        % full cohort expression matrix. Derive 25th and 75th percentile
        % thresholds across all values in the pooled matrix.
        %
        % Using pooled thresholds ensures that iMAT classifies genes as
        % low/high expressed on a scale that is comparable across groups,
        % avoiding topology differences driven by group-level scaling alone.
        GeneExprDS = [tData.Control, tData.Disease];
        th         = quantile(GeneExprDS, [0.25, 0.75], 'all');
        th_low     = th(1);   % lower threshold: genes below this are 'low expressed'
        th_high    = th(2);   % upper threshold: genes above this are 'high expressed'

        % --------------------------------------------------------------
        % 2. PROCESS HEALTHY / CONTROL GROUP  (parallelised)
        % --------------------------------------------------------------
        nCtrl       = size(tData.Control, 2);
        redModelsH  = cell(1, nCtrl);
        healthyMod  = zeros(numel(model.rxns), nCtrl);

        parfor j = 1:nCtrl

            % Build per-sample gene expression struct
            exprH      = struct('gene',  tData.Genes, ...
                                'value', tData.Control(:, j));

            % Map gene-level expression to reaction-level scores
            [expH, ~]  = mapExpressionToReactions(model, exprH);

            % Extract context-specific model via iMAT
            hmodelInfo = iMATCustom(model, expH, th_low, th_high);

            % Store reduced model (unused genes removed to save memory)
            redModelsH{j} = removeUnusedGenes(hmodelInfo);

            % Binary reaction-presence vector: 1 if reaction is active in
            % this sample's model, 0 otherwise (preserves model.rxns ordering)
            healthyMod(:, j) = ismember(model.rxns, hmodelInfo.rxns);

        end % control parfor

        % --------------------------------------------------------------
        % 3. PROCESS DISEASE GROUP  (parallelised)
        % --------------------------------------------------------------
        nDis       = size(tData.Disease, 2);
        redModelsD = cell(1, nDis);
        diseaseMod = zeros(numel(model.rxns), nDis);

        parfor j = 1:nDis

            % Build per-sample gene expression struct
            exprD      = struct('gene',  tData.Genes, ...
                                'value', tData.Disease(:, j));

            % Map gene-level expression to reaction-level scores
            [expD, ~]  = mapExpressionToReactions(model, exprD);

            % Extract context-specific model via iMAT
            hmodelInfo = iMATCustom(model, expD, th_low, th_high);

            % Store reduced model (unused genes removed to save memory)
            redModelsD{j} = removeUnusedGenes(hmodelInfo);

            % Binary reaction-presence vector (same convention as healthyMod)
            diseaseMod(:, j) = ismember(model.rxns, hmodelInfo.rxns);

        end % disease parfor

        % --------------------------------------------------------------
        % 4. STORE RESULTS FOR THIS DATASET
        % --------------------------------------------------------------
        personalizedModels(i).redModelsH = redModelsH;
        personalizedModels(i).healthyMod = healthyMod;
        personalizedModels(i).redModelsD = redModelsD;
        personalizedModels(i).diseaseMod = diseaseMod;

    end % dataset loop

end % generatePersonalizedModels