%%% 22 April 2013
%%% a function to construct GCNs based on a given set of genes, possible
%%% sets of genes are: 'ALL', 'ASD', 'SZ' and 'ND'

function createGCN(geneSet)

filesDirectory = 'files/';
resultsDirectory = 'results/';

%%% [1] Create a GCNs as the correlation between each pair of structures
%%% for each donor separately based on the expression of the geneSet
load([filesDirectory 'donorsExpMat_5RPKM_processed.mat']);

geneSetFolder = [resultsDirectory geneSet '/'];
if ~exist(geneSetFolder, 'dir')
    mkdir(geneSetFolder);
end

%%% Select the corresponding geneSet
if ~strcmp(geneSet,'ALL')    
    gInd = getGeneInd(geneSet);
    for i = 1 : size(donorsExpMat_5RPKM, 3) 
        geneMat = donorsExpMat_5RPKM(gInd,:,i);  
        RHO = corr(geneMat, 'Type', 'Spearman');
        clear geneMat;
    %     RHO = 0.5*log((1+RHO) ./ (1-RHO));
        RHO(find(isinf(RHO) == 1)) = 1;
        outFolder = [geneSetFolder 'strCorr_donors/'];
        if ~exist(outFolder, 'dir')
            mkdir(outFolder);
        end
        save([outFolder '/donor' num2str(i) '_corrMat.mat'], 'RHO');
    %     csvwrite([outFile '.csv'], RHO, 1, 1);
        clear RHO;
    end
else
    for i = 1 : size(donorsExpMat_5RPKM, 3) 
        geneMat = donorsExpMat_5RPKM(:,:,i);  
        RHO = corr(geneMat, 'Type', 'Spearman');
        clear geneMat;
    %     RHO = 0.5*log((1+RHO) ./ (1-RHO));
        RHO(find(isinf(RHO) == 1)) = 1;
        outFolder = [geneSetFolder 'strCorr_donors/'];
        if ~exist(outFolder, 'dir')
            mkdir(outFolder);
        end
        save([outFolder '/donor' num2str(i) '_corrMat.mat'], 'RHO');
    %     csvwrite([outFile '.csv'], RHO, 1, 1);
        clear RHO;
    end
end

%%% [2] Calculate the average GCN per age Group
group{1} = [1, 2, 3 ,4];
group{2} = [5, 6, 7, 8];
group{3} = [9, 10, 11, 12];
group{4} = [13, 14, 15, 16];
group{5} = [17, 18, 19, 20];
group{6} = [21, 22, 23, 24];
group{7} = [25, 26, 27, 28, 29, 30];
for i = 1 : length(group)
    clear tg;
    tg = group{i}; 
    for j = 1 : length(tg)
        load([geneSetFolder 'strCorr_donors/donor' num2str(tg(j)) '_corrMat.mat']);
        corrVol(:,:,j) = RHO;
        clear RHO;
    end
    corrMat = mean(corrVol, 3);
    clear corrVol;
    outFolder = [geneSetFolder 'strCorr_ageStages/'];
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end
    save([outFolder '/ageStage' num2str(i) '_corrMat.mat'], 'corrMat');
%     csvwrite([outFile '.csv'], corrMat, 1, 1);
    clear corrMat;
end



