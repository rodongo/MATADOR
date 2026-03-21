function antiMetab = mainGradedMATADOR(model, simMetab, KOrxn, rxnFBS, gradeLevel, VRefDis)
% mainGradedMATADOR  Core MATADOR scoring engine for a single metabolite at one grade level.
%
%   antiMetab = mainGradedMATADOR(model, simMetab, KOrxn, rxnFBS, gradeLevel, VRefDis)
%
%   Performs a metabolite-centric Metabolic Transformation simulation by
%   proportionally restricting the flux bounds of all reactions that consume
%   (or produce) a target metabolite, then solving a Minimisation of Metabolic
%   Adjustment (MOMA) quadratic programme (QP). The resulting flux distribution
%   is compared to a healthy-state reference via the MATADOR scoring function
%   (MATADOR_TS), yielding a scalar score that reflects how well the perturbed network
%   can recover towards the reference phenotype.
%
%   Reversible reactions are treated asymmetrically: only the consuming
%   direction (positive flux) is blocked; the producing direction is scaled by
%   the grade level to preserve thermodynamic feasibility.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   model       - (struct)  COBRA-compatible metabolic model containing:
%                             .rxns   – cell array of reaction identifiers
%                             .S      – stoichiometric matrix (m × n)
%                             .lb     – lower flux bounds (n × 1)
%                             .ub     – upper flux bounds (n × 1)
%                             .rev    – reversibility flags (n × 1, 1 = reversible)
%
%   simMetab    - (string)  Name of the metabolite being simulated.
%                           Used for progress reporting only.
%
%   KOrxn       - (cell or char) Reaction ID(s) to be perturbed. These are
%                           the reactions associated with the target metabolite
%                           (consuming and/or producing). The function
%                           internally converts them to model indices.
%
%   rxnFBS      - Forward - Backward - Unchanged (+1;0;-1) for each reaction
%                           based on differential gene expression and
%                           reference flux distribution. The output from
%                           diffexprs2rxnFBS function in cobratoolbox
%
%   gradeLevel  - (scalar)  Fractional inhibition level in [0, 1], where 0
%                           means complete inhibition and 1 means no inhibition.
%                           Flux bounds are scaled by this factor.
%                           NOTE: values such as 0, 0.25, 0.50, 0.75 passed from
%                           the calling function (gradedMATADOR) represent
%                           percentages and are used as-is here — ensure
%                           consistent scaling in the caller.
%
%   VRefDis     - (numeric) Reference (healthy) flux distribution vector
%                           (n × 1). Used as the MOMA target and as the
%                           baseline for MATADOR score computation.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   antiMetab   - (struct)  Structure with the following fields:
%       .MATADOR_Score  – (scalar) Similarity score between the perturbed
%                          flux state and the healthy reference. Returns
%                          -Inf if the QP is infeasible or the perturbed
%                          solution is numerically degenerate (norm < 1).
%       .fluxMATADOR    – (numeric) Optimal flux vector (n × 1) returned by
%                          the MOMA QP solver. Returns 0 if infeasible.
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW
% -------------------------------------------------------------------------
%   1. Convert reaction IDs (KOrxn) to model indices.
%   2. Scale flux bounds of target reactions by gradeLevel.
%      - Irreversible reactions: both lb and ub scaled.
%      - Reversible reactions:   ub set to 0 (consuming direction blocked);
%                                lb scaled by gradeLevel (producing direction).
%   3. Build a MOMA QP from the perturbed model:
%        minimise  ||v - VRefDis||^2
%        subject to  Sv = 0,  lb ≤ v ≤ ub
%   4. Solve the QP via solveCobraQPCustom.
%   5. Compute the MATADOR score via aTS if a feasible solution is found.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   - A flux solution with norm < 1 is treated as numerically degenerate
%     and assigned a score of -Inf. This threshold assumes typical flux
%     magnitudes on the order of 1e4; adjust if working with rescaled models.
%   - The commented-out Spearman correlation line is retained as a reference
%     to an alternative scoring approach explored during development.
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   buildLPproblemFromModel  – COBRA Toolbox utility
%   solveCobraQPCustom       – custom QP solver wrapper (must be on path)
%   MATADOR_TS               – MATADOR scoring function (must be on path)
%
% -------------------------------------------------------------------------
% AUTHORS / VERSION
% -------------------------------------------------------------------------
%   MATADOR Toolbox  |  See accompanying documentation for citation details.
% =========================================================================

    % ------------------------------------------------------------------
    % 1. INITIALISE OUTPUT AND PROGRESS REPORTING
    % ------------------------------------------------------------------
    fprintf('==============================================================================\n');
    fprintf('======================= Graded MATADOR ======================================\n');
    fprintf('==============================================================================\n');
    fprintf('======== %s metabolite simulation in progress ==============================\n', simMetab);
    fprintf('==============================================================================\n');

    antiMetab = struct();

    % ------------------------------------------------------------------
    % 2. MAP REACTION IDs TO MODEL INDICES AND COMPUTE PERTURBED BOUNDS
    % ------------------------------------------------------------------
    % Convert reaction identifiers to positional indices in model.rxns.
    KOrxn = find(ismember(model.rxns, KOrxn));

    % Scale the original bounds by the grade level (fractional inhibition).
    % A gradeLevel of 0 leaves bounds unchanged; 1 fully blocks the reaction.
    exprLevelUB = model.ub(KOrxn) * gradeLevel;
    exprLevelLB = model.lb(KOrxn) * gradeLevel;

    % ------------------------------------------------------------------
    % 3. ASYMMETRIC TREATMENT OF REVERSIBLE REACTIONS
    % ------------------------------------------------------------------
    % For reversible reactions, block only the consuming direction (ub → 0)
    % while scaling the producing direction (lb) by the grade level.
    % This preserves network feasibility and avoids artificially preventing
    % the cell from synthesising the metabolite.
    for i = 1:length(KOrxn)
        if model.rev(KOrxn(i)) == 1
            exprLevelUB(i) = 0;                       % block consuming flux
            exprLevelLB(i) = -1000 * gradeLevel;      % scale producing flux
        end
    end

    % ------------------------------------------------------------------
    % 4. BUILD AND CONFIGURE THE MOMA QP
    % ------------------------------------------------------------------
    % Construct the base LP structure from the model, then augment it into
    % a QP that minimises the Euclidean distance to the reference distribution:
    %
    %   minimise  ||v - VRefDis||^2  =  v'*F*v  +  c'*v  + const
    %   subject to  A*v = b,  lb ≤ v ≤ ub
    %
    % where F = 2*I  and  c = -2*VRefDis  (the constant is ignored).

    fprintf('========================================================================\n');
    fprintf('====================== MOMA QP in progress =============================\n');
    fprintf('========================================================================\n');

    timerVal = tic;

    QPproblem        = buildLPproblemFromModel(model);
    [~, nRxns]       = size(model.S);
    [~, nVars]       = size(QPproblem.A);

    % Linear term: derivative of ||v - VRef||^2 w.r.t. v gives -2*VRef
    QPproblem.c(1:nRxns)           = -2 * VRefDis;

    % Quadratic term: Hessian of ||v - VRef||^2 is 2*I (sparse for efficiency)
    QPproblem.F                    = sparse(nVars, nVars);
    QPproblem.F(1:nRxns, 1:nRxns) = 2 * speye(nRxns);

    QPproblem.osense = +1;   % minimisation

    % ------------------------------------------------------------------
    % 5. APPLY PERTURBED BOUNDS AND SOLVE
    % ------------------------------------------------------------------
    % Override the bounds for the target reactions with the grade-scaled
    % values computed in Steps 2–3, then solve the QP.
    momaQPProblem           = QPproblem;
    momaQPProblem.ub(KOrxn) = exprLevelUB;
    momaQPProblem.lb(KOrxn) = exprLevelLB;

    MOMAsolution = solveCobraQPCustom(momaQPProblem);

    % ------------------------------------------------------------------
    % 6. COMPUTE MATADOR SCORE FROM QP SOLUTION
    % ------------------------------------------------------------------
    VRefInd        = struct();
    VRefInd.mMTA   = [];
    MATADOR_Score  = [];

    if ~strcmp(MOMAsolution.status, 'INFEASIBLE')

        VRefInd.mMTA = MOMAsolution.x;

        if ~isempty(KOrxn) && norm(VRefInd.mMTA) < 1
            % Degenerate solution: flux norm near zero indicates the
            % perturbation has effectively shut down the network.
            % Assign worst-case score rather than propagate a misleading value.
            MATADOR_Score = -Inf;
        else
            % Compute similarity between perturbed and reference flux states.
            % Inactive reactions (zero flux in both states) are excluded by aTS.
            MATADOR_Score = MATADOR_TS(VRefInd.mMTA, VRefDis, rxnFBS);
        end

    else
        % QP is infeasible under the imposed perturbation; assign worst-case score.
        VRefInd.mMTA  = 0;
        MATADOR_Score = -Inf;
    end

    time_MOMA = toc(timerVal);

    fprintf('================================================================================\n');
    fprintf('\t MOMA QP solved in: %4.2f seconds (%4.2f minutes)\n', time_MOMA, time_MOMA / 60);
    fprintf('================================================================================\n');

    % ------------------------------------------------------------------
    % 7. PACKAGE AND RETURN RESULTS
    % ------------------------------------------------------------------
    antiMetab.MATADOR_Score  = MATADOR_Score;
    antiMetab.fluxMATADOR    = MOMAsolution.x;

end % mainGradedMATADOR