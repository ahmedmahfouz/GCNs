%%% 28 Nov. 2013
%%% create a separate GCN per brain region based on the list of selected
%%% structures and donors

function [donorGCN donorExpMat] = gcnPerRegion(donor,structure,...
    gene,expMatrix,expressingGenesIDs,structures,donors,geneSet,geneListFile,...
    corrType,devPeriods,devPeriodsNames)


%% if geneSet is NOT 'All'
if ~strcmpi(geneSet, 'All')
    geneSetInd = geneSetInBrainSpan(geneListFile,geneSet,gene);
    expressingGenesIDs = expressingGenesIDs(find(ismember(expressingGenesIDs,geneSetInd) == 1));
end

if strcmpi(corrType,'Pearson')
    expMatrix = log2(expMatrix+1);
end

%% construct a donor expression matrix (genes,structures,donors)
donorExpMat.expression = NaN(length(expressingGenesIDs),length(structures),length(donors));
for d = 1 : length(donors)
    % find the ids of the current donors
    donorIDs = find(ismember(donor.name,donors{d}) == 1);
    % find the ids of the selected structures corresponding to the current
    % donor
    allDonorStructures = structure.acronym(donorIDs);
    % rearrange donor structures to match the order of the structuresOfInteres
    [tempSTR ia ib] = intersect(structures,allDonorStructures);
    [sortedIA sortingInd] = sort(ia);
    donorStrIDs = donorIDs(ib(sortingInd));
    % create a donor-specific expression matrix
    donorExpMat.expression(:,sortedIA,d) = expMatrix(expressingGenesIDs, donorStrIDs);
    % retrieve a list of ages corresponding to the list of donors
    donorExpMat.ages(d) = donor.age(donorIDs(1));
    clear donorIDs; clear allDonorStructures; clear tempSTR;
    clear ia; clear ib; clear sortedIA; clear sortingInd; clear donorStrIDs;
end

%% impute missing values
donorExpMat = customImpute(donorExpMat,devPeriods);

%% create a region-specific GCN
figure, hold on
for d = 1 : length(structures)
    RHO(:,d,:) = corr(squeeze(donorExpMat.expression(:,d,:)), 'Type', corrType);
    
    subplot(4,4,d), imagesc(squeeze(RHO(:,d,:))),
    title(structures{d}, 'FontWeight', 'bold');
    set(gca, 'YTickLabels', donorExpMat.ages, 'YTick', [1:length(donorExpMat.ages)], 'FontSize', 5)
    set(gca, 'XTickLabels', donorExpMat.ages, 'XTick', [1:length(donorExpMat.ages)], 'FontSize', 5)
    rotateXLabels(gca(), 45 );
    
    currRHO = tril(squeeze(RHO(:,d,:)), -1);
    currRHO = reshape(currRHO, size(currRHO,1)*size(currRHO,2), 1);
    currRHO(find(currRHO == 0)) = [];
    donorGCN(:,d) = currRHO;
    clear currRHO;
end
hold off


