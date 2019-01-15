%%% 28 Nov. 2013
%%% create a separate GCN per donor brain based on the list of selected
%%% structures and donors for a random list of genes

function [donorGCN devPeriodGCN strPairs donorExpMat] = gcnPerDonor_random(donor,structure,...
    gene,expMatrix,expGeneIds,structures,donors,...
    corrType,devPeriods,devPeriodsNames, randomGenes)

%% extract gene indicies
randGeneInd = randi(length(expGeneIds),randomGenes,1); 
expGeneIds = randGeneInd;

%% log-transform the matrix if using Pearson's Correlation
if strcmpi(corrType,'Pearson')
    expMatrix = log2(expMatrix+1);
end

%% construct a donor expression matrix (genes,structures,donors)
donorExpMat.expression = NaN(length(expGeneIds),length(structures),length(donors));
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
    donorExpMat.expression(:,sortedIA,d) = expMatrix(expGeneIds, donorStrIDs);
    % retrieve a list of ages corresponding to the list of donors
    donorExpMat.ages(d) = donor.age(donorIDs(1));
    clear donorIDs; clear allDonorStructures; clear tempSTR;
    clear ia; clear ib; clear sortedIA; clear sortingInd; clear donorStrIDs;
end

%% impute missing values
donorExpMat = customImpute(donorExpMat,devPeriods);

%% create a donor-specific GCN
% figure, hold on
for d = 1 : length(donors)
    RHO(:,:,d) = corr(donorExpMat.expression(:,:,d), 'Type', corrType);
    
%     subplot(5,6,d), imagesc(RHO(:,:,d)),
%     title(donors{d}, 'FontWeight', 'bold');
%     set(gca, 'YTickLabels', structures, 'YTick', [1:length(structures)], 'FontSize', 5)
%     set(gca, 'XTickLabels', structures, 'XTick', [1:length(structures)], 'FontSize', 5)
%     rotateXLabels(gca(), 45 );
    
    currRHO = tril(RHO(:,:,d), -1);
    currRHO = reshape(currRHO, size(currRHO,1)*size(currRHO,2), 1);
    currRHO(find(currRHO == 0)) = [];
    donorGCN(:,d) = currRHO;
    clear currRHO;
end
% hold off

%% create developmental-period specific GCN
% figure, hold on
for  dp = 1 : length(devPeriods)
    currDonors = devPeriods{dp};
    devPeriodRHO(:,:,dp) = mean(RHO(:,:,currDonors),3);
    
%     subplot(1,7,dp), imagesc(devPeriodRHO(:,:,dp)),
%     title(devPeriodsNames{dp}, 'FontWeight', 'bold');
%     set(gca, 'YTickLabels', structures, 'YTick', [1:length(structures)], 'FontSize', 5)
%     set(gca, 'XTickLabels', structures, 'XTick', [1:length(structures)], 'FontSize', 5)
%     rotateXLabels(gca(), 45 );
    
    devPeriodGCN(:,dp) = mean(donorGCN(:,currDonors)');
end
% hold off

%% generate structure-pair names
count = 0;
for i = 1 : length(structures)-1
    for j = i+1 : length(structures)
        if i ~= j
            count = count + 1;
            strPairs{count} = [structures{i} '_' structures{j}];
        end
    end
end
