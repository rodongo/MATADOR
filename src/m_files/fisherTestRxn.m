% Brief summary of this function.
% Detailed explanation of this function.
function fisherTestRxn(model, personalizedModels,saveLocation)
for i = 1:numel(personalizedModels)
    conBIN = personalizedModels(i).healthyMod;
    adBIN = personalizedModels(i).diseaseMod;
    
    fishersPVal = cell(1, size(conBIN, 1));
    fishersOR = cell(1, size(conBIN, 1));
    
    parfor l = 1:size(conBIN, 1)
        statData = [sum(adBIN(l,:)==1), sum(conBIN(l,:)==1); ...
                    sum(adBIN(l,:)==0), sum(conBIN(l,:)==0)];
        [~, p, OR] = fishertest(statData);
        fishersPVal{l} = p;
        fishersOR{l} = OR.OddsRatio;
    end
    
    % Save Results
    statTable = table(model.rxns, [fishersPVal{:}]', [fishersOR{:}]', ...
                round(mean(conBIN, 2)*100), round(mean(adBIN, 2)*100), ...
                'VariableNames', {'Reaction','PValue','OR','ActiveCon','ActiveAD'});
                
    writetable(statTable, sprintf(str(saveLocation),'AMS_iMAT_FisherTest.xlsx'), 'Sheet', personalizedModels(i).OmicData);
end
end