% MATADOR: Anti-metabolite Computational Pipeline
% This script executes the full MATADOR workflow for sensitivity analysis at various flux grades:
% Environment initialization and toolbox setup.
% Loading expression & DEG data.
% generating contextualized/personalized metabolic models using iMAT.
% Statistical test for reaction activity.
% Sampling metabolic networks and defining set of reactions that must change to restore 'healthy' state.
% Anti-metabolite simulation (Average).
% Requirements:
%     - MATLAB (R2023b v23.2)
%     - Gurobi (v11.0)
%     - COBRATOOLBOX (3.1)
%     - RAVEN Toolbox (v2.8.7)

% PART 1: Setting up environment
clearvars;
addpath('./src/')
addpath('./cobratoolbox-master')
addpath('./RAVEN-main/installation')
addpath('./Human-GEM-main/GPRs')
addpath('./Human-GEM-main/tINIT')

initCobraToolbox
checkInstallation
cd('.')

% PART 2: Load expression and differential gene expression data

% Load expression data
% The processed gene expression data should be in an Excel worksheet with 
% "Controls" samples in sheet 1 and "Disease" samples in sheet 2
expressionData = struct;
fileList1 = dir(fullfile("Expression Data", '*.xlsx'));
for i=1:numel(fileList1)
%     i = 1;
    expressionData(i).OmicData = strrep(fileList1(i).name,'.xlsx','');
    sheets = sheetnames(append(fileList1(i).folder,'/', fileList1(i).name));
    data1 = readtable(append(fileList1(i).folder,'/', fileList1(i).name),'Sheet',sheets{1});
    expressionData(i).Genes = string(data1.GeneID);
    expressionData(i).Control = table2array(data1(:,2:end));
    data1 = readtable(append(fileList1(i).folder,'/', fileList1(i).name),'Sheet',sheets{2});
    expressionData(i).Disease = table2array(data1(:,2:end));
end
clear i data1 genes sheets fileList1
save('expressionData.mat','expressionData')
% Load differential gene expression data
degData_Average = struct;
fileList1 = dir(fullfile("DEG Data/", '*.xlsx'));
for i=1:numel(fileList1)
    degData_Average(i).OmicData = strrep(fileList1(i).name, '_DEG.xlsx','');
    diffExpr = readtable(append(fileList1(i).folder,'/', fileList1(i).name));
    diffExpr = diffExpr(:,[1,2,3]);
    diffExpr.Properties.VariableNames = {'gene','logFC','pval'};
    degData_Average(i).diffExpr= diffExpr;
end
clear fileList1 diffExpr i sheets

% Personalized DEG data
SampleName = readtable('Expression Data/ROSMAP.xlsx','ReadVariableNames',true,'Sheet','AD');
SampleName = SampleName.Properties.VariableNames;
SampleName(1) = [];
degData_Personalized = struct;
controlExp = AnalysisData.Control;
for i=1:size(AnalysisData.Disease,2)
    aveExpr = mean([AnalysisData.Disease(:,i) controlExp],2);
    sdExpr = std([AnalysisData.Disease(:,i) controlExp],0,2);
    zExpr = (AnalysisData.Disease(:,i) - aveExpr)./sdExpr;
    pTwoTailed = 2 * (1 - normcdf(abs(zExpr)));
    degData_Personalized(i).PatientID = SampleName(i);
    diffExpr = table(AnalysisData.Genes,zExpr,pTwoTailed);
    diffExpr.Properties.VariableNames = {'gene','logFC','pval'};
    degData_Personalized(i).diffExpr= diffExpr;
end
clear fileList1 diffExpr i sheets
save('degData_Personalized.mat','degData_Personalized')

% PART 3: Prepare metabolic models

% Load genome scale metabolic model and set constraints
gemPath = './Human-GEM-main/Human-GEM.mat';
ihuman = load(gemPath).ihuman;
curatedModel = ravenCobraWrapper(ihuman);

% Define physiological constraints
% biomass (MAR13082), oxygen (MAR09048), glucose (MAR09034)
constraints = struct('biomass', 0.0001, 'oxygen', -0.01, 'glucose', -0.01);

curatedModel.lb(strcmp('MAR13082', curatedModel.rxns)) = constraints.biomass;
curatedModel.ub(strcmp('MAR09048', curatedModel.rxns)) = constraints.oxygen;
curatedModel.ub(strcmp('MAR09034', curatedModel.rxns)) = constraints.glucose;

% PART 4: Generation of personalized metabolic networks

% Generation of personalized metabolic networks
% This step creates a contexualized metabolic network and binary matrix for
% reactions that remain active "1" and those removed by iMAT "0"
rng(2304.2024,'twister')
% Configure Parallel Pool based on hardware capacity
if isempty(gcp('nocreate')), parpool(18); end
changeCobraSolver('gurobi','all');

personalizedMetNet_iMAT = generatePersonalizedModels(curatedModel, expressionData);

save('personalizedMetNet_iMAT.mat','personalizedMetNet_iMAT')

% PART 5: Perform statistical tests on reaction activity

% Statistical Tests: Fisher's Exact Test
% Fisher's exact test is performed to test reaction activity and condition
% i.e disease. The test results are saved to an Excel file in the given
% director. the table should be in the form of:
% ##############################################################
% ########## Reaction PValue  OR  ActiveCon   ActiveAD #########
% ##############################################################
saveLocation = './';
fisherTestRxn(curatedModel, personalizedMetNet_iMAT, saveLocation)

% PART 6: Extract metabolic targets from summary statistics

% Generation of Metabolite Targets
% This section relies on R scripts. Refer to prepareMetaboliteTargetList.R
% script to generate metabolite targets from reaction summary statistics

% PART 7: Sampling metabolic networks

% Sampling Metabolic Networks
changeCobraSolver('gurobi','all');
% a) Based on average data
averageSampledModels = sampleMetNetworks(curatedModel, personalizedMetNet_iMAT, degData_Average, 'AVERAGE',10);
% b) Based on personalized
personalizedSampledModels = sampleMetNetworks(curatedModel, personalizedMetNet_iMAT, degData_Personalized, 'PERSONALIZED',10);


% PART 8: Run Graded MATADOR

% Configuration
rng(1107.2025, 'twister');
targetPath = 'Metabolite Reports/Metabolite Targets/';
resultsPath = 'Metabolite Reports/Antimetabolite Analysis/';
if ~exist(resultsPath, 'dir'), mkdir(resultsPath); end

% Grades
gradeLevels = [0 0.25 0.50 0.75];
% --- Average Analysis ---
fprintf('Executing Graded Average MATADOR for Sensitivity Analysis...\n');
gradedMATADOR_Sensitivity = struct();
for i = 1:numel(averageSampledModels)
    datasetName = AntiMetCompSampTS(i).OmicData;
    targetFile = fullfile(targetPath, [datasetName, '.xlsx']);
    
    gradedMATADOR_Sensitivity(i).DataSet = datasetName;
    gradedMATADOR_Sensitivity(i).MetaboliteSims = gradedMATADOR(...
        averageSampledModels(i).disModel, targetFile, ...
        averageSampledModels(i).diffFBS, averageSampledModels(i).VRefDis, gradeLevels, ...
        resultsPath, datasetName);
end
save('gradedMATADOR_Sensitivity.mat','gradedMATADOR_Sensitivity')