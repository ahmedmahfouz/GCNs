%%%--Correlation between the expression pattern of all pairs of structures 
%%% at different age groups for a specific gene list (All, ASDs, SZ, ND)

clear all;

filesDirectory = 'files/';
resultsDirectory = 'results/';

%%% Calculate the structure correlation networks for each donor separately

dName = 'ASDs';
load([filesDirectory dName 'GeneIDs_5RPKM']);
dgIND = dIDs_T;

load([filesDirectory 'donorsExpMat_5RPKM.mat']);
donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));

dGeneMat = donorsExpMat_5RPKM(dgIND,:,:);
clear donorsExpMat_5RPKM;

for i = 1 : size(dGeneMat, 3)
    clear geneMat; 
    geneMat = dGeneMat(:,:,i);  
    RHO = corr(geneMat, 'Type', 'Spearman');
    RHO = 0.5*log((1+RHO) ./ (1-RHO));
    RHO(find(isinf(RHO) == 1)) = 1;
    
    outFile = [resultsDirectory 'strCorr\' dName '\Avg(log2+epsilon)\donor' num2str(i)];
    save(outFile, 'RHO');
    csvwrite([resultsDirectory 'strCorr\' dName '\Avg(log2+epsilon)\donor' num2str(i) '_corrMat.csv'], RHO, 1, 1);
    clear RHO;
end

%%% Calculate the average gene correlation per age Group
group{1} = [1, 2, 3 ,4];
group{2} = [5, 6, 7, 8];
group{3} = [9, 10, 11, 12];
group{4} = [13, 14, 15, 16];
group{5} = [17, 18, 19, 20];
group{6} = [21, 22, 23, 24];
group{7} = [25, 26, 27, 28, 29, 30];

folder = [resultsDirectory 'strCorr\' dName '\Avg(log2+epsilon)\'];
for i = 1 : length(group)
    
    clear tg;
    tg = group{i}; 
    for j = 1 : length(tg)
        load([folder 'donor' num2str(tg(j))]);
        corrVol(:,:,j) = RHO;
        clear RHO;
    end
    corrMat = mean(corrVol, 3);
    clear corrVol;
     
    save([resultsDirectory 'strCorr\' dName '\Avg(log2+epsilon)\ageGroup' num2str(i) '_corrMat.mat'], 'corrMat');
    csvwrite([resultsDirectory 'strCorr\' dName '\Avg(log2+epsilon)\ageGroup' num2str(i) '_corrMat.csv'], corrMat, 1, 1);
    clear corrMat;

end
%%%------------------------------------------------------------------------

%%% Calculate the gene correlation networks for concatenated samples
% load('files\newGeneMat_refined.mat');
% 
% for i = 1 : length(newGeneMat)
%    
%     geneMat = newGeneMat{i}; 
%     
%     for j = 1 : size(geneMat,3)
%         
%         if j == 1
%             G = geneMat(:,:,j);
%         else
%             G = [G geneMat(:,:,j)];
%         end
%         
%     end
%     
%     clear geneMat;
%     G = permute(G, [2 1 3]);
%     [RHO, PVAL] = corr(G, 'Type', 'Spearman');
%     clear G;
%     
%     save(['results\geneCorr\Concat\ageGroup' num2str(i) '_corrMat.mat'], 'RHO');
%     save(['results\geneCorr\Concat_pval\ageGroup' num2str(i) '_corrMat.mat'], 'PVAL');
%     clear RHO; clear PVAL;
%     
% end
%%%------------------------------------------------------------------------



