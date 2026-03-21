function antiMetabScore = MATADOR(model, simMetab, KOrxn, rxnFBS, VRefDis)
% MATADOR  Binary metabolite-centric Metabolic Transformation simulation.
%
%   antiMetabScore = MATADOR(model, simMetab, KOrxn, rxnFBS, VRefDis)
%
%   Performs a complete block of all reactions associated with a
%   target metabolite, then solves a Minimisation of Metabolic Adjustment
%   (MOMA) quadratic programme (QP) to obtain the flux distribution that
%   best satisfies steady-state constraints while remaining as close as
%   possible to a healthy-state reference. The resulting flux vector is
%   scored by MATADOR_TS to quantify how well the network can recover
%   towards (or diverge from) the reference phenotype.
%
%   This function implements the full-inhibition (grade = 1) special case
%   of the graded simulation performed by mainGradedMATADOR. All consuming
%   reactions are fully blocked; reversible reactions retain their producing
%   direction at a fixed bound of -1000 (model units).
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
%   KOrxn       - (cell or char) Reaction ID(s) to be fully blocked. These
%                           correspond to the consuming (and/or producing)
%                           reactions of the target metabolite. The function
%                           internally converts them to model indices.
%
%   rxnFBS      - Forward - Backward - Unchanged (+1;0;-1) for each reaction
%                           based on differential gene expression and
%                           reference flux distribution. The output from
%                           diffexprs2rxnFBS function in cobratoolbox
%
%   VRefDis     - (numeric) Reference (healthy) flux distribution vector
%                           (n × 1). Serves as both the MOMA optimisation
%                           target and the baseline for score computation.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   antiMetabScore - (struct)  Structure with the following fields:
%       .MATADOR_Score  – (scalar) MATADOR similarity score between the
%                          fully-perturbed flux state and the healthy
%                          reference. Returns -Inf if the QP is infeasible
%                          or the solution is numerically degenerate
%                          (flux norm < 1).
%       .fluxMOMA       – (numeric) Optimal flux vector (n × 1) from the
%                          MOMA QP. Returns 0 if infeasible.
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW
% -------------------------------------------------------------------------
%   1. Initialise perturbed bounds to zero for all target reactions.
%   2. For reversible reactions: ub → 0 (block consuming direction);
%                                lb → -1000 (allow producing direction).
%   3. Build a MOMA QP:
%        minimise  ||v - VRefDis||^2
%        subject to  Sv = 0,  lb ≤ v ≤ ub
%   4. Solve via solveCobraQPCustom.
%   5. Score the solution with MATADOR_TS if feasible and non-degenerate.
%
% -------------------------------------------------------------------------
% RELATIONSHIP TO GRADED VARIANT
% -------------------------------------------------------------------------
%   MATADOR performs a single, complete inhibition (equivalent to
%   gradeLevel = 1 in mainGradedMATADOR), whereas mainGradedMATADOR sweeps
%   across fractional inhibition levels. Key implementation differences:
%     - Bounds are initialised to zero here (not scaled from model.ub/lb).
%     - Uses buildOptProblemFromModel (vs. buildLPproblemFromModel).
%     - Calls MATADOR_TS (vs. aTS) for scoring.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   - The flux-norm threshold of 1 assumes typical network flux magnitudes
%     on the order of 1e4. Adjust if working with a rescaled model.
%   - The producing-direction bound of -1000 for reversible reactions is a
%     fixed heuristic; consider parameterising if model units differ.
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   buildOptProblemFromModel  – COBRA Toolbox utility
%   solveCobraQPCustom        – custom QP solver wrapper (must be on path)
%   MATADOR_TS                – MATADOR scoring function (must be on path)
%
% -------------------------------------------------------------------------
% AUTHORS / VERSION
% -------------------------------------------------------------------------
%   MATADOR Toolbox  |  See accompanying documentation for citation details.
% =========================================================================

    % ------------------------------------------------------------------
    % 1. INITIALISE OUTPUT AND PROGRESS REPORTING
    % ------------------------------------------------------------------
    fprintf('====================================================\n');
    fprintf('==================== MATADOR =======================\n');
    fprintf('====================================================\n');
    fprintf('===== %s metabolite simulation in progress =========\n', simMetab);
    fprintf('====================================================\n');

    antiMetabScore = struct();

    % ------------------------------------------------------------------
    % 2. INITIALISE PERTURBED BOUNDS AND MAP REACTION IDs TO INDICES
    % ------------------------------------------------------------------
    % Default perturbed bounds are zero for all target reactions,
    % representing complete inhibition of flux in both directions.
    exprLevelUB = zeros(1, length(KOrxn));
    exprLevelLB = zeros(1, length(KOrxn));

    % Convert reaction identifiers to positional indices in model.rxns.
    KOrxn = find(ismember(model.rxns, KOrxn));

    % ------------------------------------------------------------------
    % 3. ASYMMETRIC TREATMENT OF REVERSIBLE REACTIONS
    % ------------------------------------------------------------------
    % For reversible reactions, block only the consuming direction (ub → 0)
    % while permitting flux in the producing direction up to a fixed limit.
    % This prevents the network from becoming infeasible due to thermodynamic
    % constraints while still eliminating net consumption of the metabolite.
    for i = 1:length(KOrxn)
        if model.rev(KOrxn(i)) == 1
            exprLevelUB(i) = 0;       % block consuming (forward) flux
            exprLevelLB(i) = -1000;   % permit producing (reverse) flux
        end
    end

    % ------------------------------------------------------------------
    % 4. BUILD AND CONFIGURE THE MOMA QP
    % ------------------------------------------------------------------
    % Construct the optimisation problem from the model, then configure it
    % as a QP minimising the squared Euclidean distance to the reference:
    %
    %   minimise  ||v - VRefDis||^2  =  v'*F*v  +  c'*v  + const
    %
    % where  F = 2*I  (Hessian)  and  c = -2*VRefDis  (linear term).
    % The constant is independent of v and omitted from the objective.

    fprintf('========================================================================\n');
    fprintf('====================== MOMA QP in progress =============================\n');
    fprintf('========================================================================\n');

    timerVal = tic;

    QPproblem        = buildOptProblemFromModel(model);
    [~, nRxns]       = size(model.S);
    [~, nVars]       = size(QPproblem.A);

    QPproblem.c(1:nRxns)           = -2 * VRefDis;          % linear term
    QPproblem.F                    = sparse(nVars, nVars);   % allocate sparse Hessian
    QPproblem.F(1:nRxns, 1:nRxns) = 2 * speye(nRxns);       % quadratic term
    QPproblem.osense               = +1;                     % minimisation

    % ------------------------------------------------------------------
    % 5. APPLY PERTURBED BOUNDS AND SOLVE
    % ------------------------------------------------------------------
    % Copy the QP, override bounds for the target reactions, and solve.
    momaQPProblem           = QPproblem;
    momaQPProblem.ub(KOrxn) = exprLevelUB;
    momaQPProblem.lb(KOrxn) = exprLevelLB;

    MOMAsolution = solveCobraQPCustom(momaQPProblem);

    % ------------------------------------------------------------------
    % 6. COMPUTE MATADOR SCORE FROM QP SOLUTION
    % ------------------------------------------------------------------
    VRefInd      = struct();
    VRefInd.mMTA = [];
    MATADOR_Score      = [];

    if ~strcmp(MOMAsolution.status, 'INFEASIBLE')

        VRefInd.mMTA = MOMAsolution.x;

        if ~isempty(KOrxn) && norm(VRefInd.mMTA) < 1
            % Degenerate solution: flux norm near zero implies the imposed
            % block has effectively silenced the entire network.
            % Assign worst-case score to flag this condition.
            MATADOR_Score = -Inf;
        else
            % Compute the MATADOR similarity score between the perturbed
            % flux distribution and the healthy reference.
            MATADOR_Score = MATADOR_TS(VRefInd.mMTA, VRefDis, rxnFBS);
        end

    else
        % QP infeasible: perturbation cannot be accommodated by the network.
        VRefInd.mMTA = 0;
        MATADOR_Score      = -Inf;
    end

    time_MOMA = toc(timerVal);

    fprintf('================================================================================\n');
    fprintf('\t MOMA QP solved in: %4.2f seconds (%4.2f minutes)\n', time_MOMA, time_MOMA / 60);
    fprintf('================================================================================\n');

    % ------------------------------------------------------------------
    % 7. PACKAGE AND RETURN RESULTS
    % ------------------------------------------------------------------
    antiMetabScore.MATADOR_Score = MATADOR_Score;
    antiMetabScore.fluxMOMA      = MOMAsolution.x;

end % MATADOR