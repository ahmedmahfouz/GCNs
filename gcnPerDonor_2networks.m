%%% 14 Oct. 2013
%%% create a separate GCN per donor brain based on the list of selected
%%% structures and donors

function [donorGCN_big6 devPeriodGCN_big6 donorGCN_ncx devPeriodGCN_ncx strPairs_big6 strPairs_ncx donorExpMat] = gcnPerDonor(donor,structure,...
    gene,expMatrix,expressingGenesIDs,structures,donors,geneSet,geneListFile,...
    corrType,devPeriods,devPeriodsNames)

%% remove non expressing genes (genes with no expression value above 5 RPKM)
% % [expressingGenesIDs ~] = ind2sub(size(origExpMatrix), find(origExpMatrix > T));
% % expressingGenesIDs = unique(expressingGenesIDs);
% % 
% % if ~strcmpi(geneSet, 'All')
% %     geneSetInd = geneSetInBrainSpan(geneListFile,geneSet,gene);
% %     expressingGenesIDs = expressingGenesIDs(find(ismember(expressingGenesIDs,geneSetInd) == 1));
% % end

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

%% creat donor GCNs based on  structures (group NCX samples) and NCX samples separately
donorExpMat_big6 = donorExpMat.expression(:,1:5,:); 
donorExpMat_big6(:,end+1,:) = squeeze(mean(donorExpMat.expression(:,6:16,:),2)); 
donorExpMat_ncx = donorExpMat.expression(:,6:16,:);
structures_big6 = [structures(1:5), {'NCX'}];
structures_ncx = structures(6:16);

%% create a donor-specific GCN
for d = 1 : length(donors)
    RHO_big6(:,:,d) = corr(donorExpMat_big6(:,:,d), 'Type', corrType);
    currRHO = tril(RHO_big6(:,:,d), -1);
    currRHO = reshape(currRHO, size(currRHO,1)*size(currRHO,2), 1);
    currRHO(find(currRHO == 0)) = [];
    donorGCN_big6(:,d) = currRHO;
    clear currRHO;
    
    RHO_ncx(:,:,d) = corr(donorExpMat_ncx(:,:,d), 'Type', corrType);
    currRHO = tril(RHO_ncx(:,:,d), -1);
    currRHO = reshape(currRHO, size(currRHO,1)*size(currRHO,2), 1);
    currRHO(find(currRHO == 0)) = [];
    donorGCN_ncx(:,d) = currRHO;
    clear currRHO;
end

%% create developmental-period specific GCN
figure, hold on
for  dp = 1 : length(devPeriods)
    currDonors = devPeriods{dp};
    devPeriodRHO_big6(:,:,dp) = mean(RHO_big6(:,:,currDonors),3);
    subplot(2,7,dp), imagesc(devPeriodRHO_big6(:,:,dp)),
    title(devPeriodsNames{dp}, 'FontWeight', 'bold');
    set(gca, 'YTickLabels', structures_big6, 'YTick', [1:length(structures_big6)], 'FontSize', 5)
    set(gca, 'XTickLabels', structures_big6, 'XTick', [1:length(structures_big6)], 'FontSize', 5)
    rotateXLabels(gca(), 45 );
    devPeriodGCN_big6(:,dp) = mean(donorGCN_big6(:,currDonors)');
    
    devPeriodRHO_ncx(:,:,dp) = mean(RHO_ncx(:,:,currDonors),3);
    subplot(2,7,dp+7), imagesc(devPeriodRHO_ncx(:,:,dp)),
    title(devPeriodsNames{dp}, 'FontWeight', 'bold');
    set(gca, 'YTickLabels', structures_ncx, 'YTick', [1:length(structures_ncx)], 'FontSize', 5)
    set(gca, 'XTickLabels', structures_ncx, 'XTick', [1:length(structures_ncx)], 'FontSize', 5)
    rotateXLabels(gca(), 45 );
%     devPeriodGCN_ncx(:,dp) = mean(donorGCN_ncx(:,currDonors)');
    tempR = tril(devPeriodRHO_ncx(:,:,dp),-1);
    tempR = reshape(tempR, size(tempR,1)*size(tempR,2), 1);
    tempR(find(tempR == 0)) = [];
    devPeriodGCN_ncx(:,dp) = tempR;
end
hold off

%% generate structure-pair names
count = 0;
for i = 1 : length(structures_big6)-1
    for j = i+1 : length(structures_big6)
        if i ~= j
            count = count + 1;
            strPairs_big6{count} = [structures_big6{i} '_' structures_big6{j}];
        end
    end
end

count = 0;
for i = 1 : length(structures_ncx)-1
    for j = i+1 : length(structures_ncx)
        if i ~= j
            count = count + 1;
            strPairs_ncx{count} = [structures_ncx{i} '_' structures_ncx{j}];
        end
    end
end


