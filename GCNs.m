%%% 4 April 2013
%%% A pipeline to generate Genomic Connectivity Networks (GCNs) based on the
%%% BrainSpan Atlas (developing human brain atlas)


%% [1] Read the processed Expression Matrix (inconsistent donors removed, 
%%% non-expressing genes are removed, only important structures areretained) 
% filesDirectory = 'files/';
% load([filesDirectory 'donorsExpMat_5RPKM.mat']);
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% save([filesDirectory 'donorsExpMat_5RPKM_processed.mat'], 'donorsExpMat_5RPKM');

%% [2] Create a GCNs as the correlation between each pair of structures
%% for each donor separately based on the expression of ALL genes (13,563)
% filesDirectory = 'files/';
% resultsDirectory = 'results/All/';
% load([filesDirectory 'donorsExpMat_5RPKM_processed.mat']);
% for i = 30 : size(donorsExpMat_5RPKM, 3) 
%     geneMat = donorsExpMat_5RPKM(:,:,i);  
%     RHO = corr(geneMat, 'Type', 'Spearman');
%     clear geneMat;
% %     RHO = 0.5*log((1+RHO) ./ (1-RHO));
%     RHO(find(isinf(RHO) == 1)) = 1;
%     outFile = [resultsDirectory 'strCorr_donors/donor' num2str(i) '_corrMat'];
%     save([outFile '.mat'], 'RHO');
% %     csvwrite([outFile '.csv'], RHO, 1, 1);
%     clear RHO;
% end

%% [3] Calculate the average GCN per age Group
% filesDirectory = 'files/';
% resultsDirectory = 'results/';
% group{1} = [1, 2, 3 ,4];
% group{2} = [5, 6, 7, 8];
% group{3} = [9, 10, 11, 12];
% group{4} = [13, 14, 15, 16];
% group{5} = [17, 18, 19, 20];
% group{6} = [21, 22, 23, 24];
% group{7} = [25, 26, 27, 28, 29, 30];
% for i = 1 : length(group)
%     clear tg;
%     tg = group{i}; 
%     for j = 1 : length(tg)
%         load([resultsDirectory 'strCorr_donors\donor' num2str(tg(j)) '_corrMat.mat']);
%         corrVol(:,:,j) = RHO;
%         clear RHO;
%     end
%     corrMat = mean(corrVol, 3);
%     clear corrVol;
%     outFile = [resultsDirectory 'strCorr_ageStages\ageStage' num2str(i) '_corrMat'];
%     save([outFile '.mat'], 'corrMat');
% %     csvwrite([outFile '.csv'], corrMat, 1, 1);
%     clear corrMat;
% end

%% [4] Cluster regions based on their pair-wise correlation change over
%%% time
% filesDirectory = 'files/';
% resultsDirectory = 'results/';
% setName = 'ASD';
% for i = 1 : 30
%     load([resultsDirectory setName '/' 'strCorr_donors/donor' num2str(i) '_corrMat']);
%     RHO = tril(RHO,-1);
%     RHO = reshape(RHO, size(RHO,1)*size(RHO,2), 1);
%     RHO(find(RHO == 0)) = [];
%     allDonorsCorrMat(:,i) = RHO;
%     clear RHO;
% end
% load(['files/strucs.mat']);
% strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24, 26];
% for i = 1 : length(strucIndInc)
%     S{i} = strucs{strucIndInc(i)};
% end
% count = 0;
% for i = 2 : length(S)
%     for j = i-1 : length(S)
%         if i ~= j
%             count = count + 1;
%             temp{count} = [num2str(i) '_' num2str(j)];
%             strPair{count} = [S{i} '_' S{j}];
%         end
%     end
% end
% clear strucIndInc; clear strucs;
% % figure, imagesc(allDonorsCorrMat)
% cT = 0.85; % cutoff threshold
% corrType = 'correlation';
% linkType = 'average'; % linkage type
% % corrDist = 1 - allDonorsCorrMat;
% % clusterTree = linkage(corrDist, linkType);
% % clusters = cluster(clusterTree, 'cutoff', cT, 'criterion', 'distance');
% % f1 = figure;
% % [H, Leafs, P] = dendrogram(clusterTree, 0, 'colorthreshold', cT, 'orientation', 'left');
% C = clustergram(allDonorsCorrMat, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', ...
%     corrType, 'Linkage', linkType, 'Dendrogram', cT, 'Colormap', 'redbluecmap', ...
%     'RowLabels', strPair);
% h = plot(C);
% saveas(h, ['results\' setName '\region-pairs clustering_' setName '_donors.fig']);

%% [5] Calculate topological measures of all the donor networks
resultsDirectory = 'C:/Users/amahfouz/Documents/MATLAB/Results/GenomicConnectivityNetworks/All/strCorr_donors/';
nDonors = 30;
for sub = 1 : nDonors
    load([resultsDirectory 'donor' num2str(sub) '_corrMat.mat']);
    
    outFile = ['donor' num2str(sub) '_netmeasures'];
    NetworkMeasures(RHO, outFile, 1, 'b')
end
%% [6] Calculate topological measures of all the group networks
% resultsDirectory = 'results/All/strCorr_ageStages/';
% nGroups = 7;
% for sub = 1 : nGroups
%     load([resultsDirectory 'ageStage' num2str(sub) '_corrMat.mat']);
%     
%     outFile = ['ageStage' num2str(sub) '_netmeasures'];
%     NetworkMeasures(corrMat, outFile, 1, 'b')
% end
%% [7] combine network measures in one matrix/file
filesDirectory = 'results/NetworkMeasures_donors/';
% % loop on all groups (All, SZ, ASD and ND)
groups = {'All', 'ASD', 'ND', 'SZ'};
for g = 1 : 1
    groupDir = [filesDirectory groups{g} '/'];
    % loop on all the subjects
    subjectsFiles = dir([groupDir, '\*.mat']);
    for s = 1 : length(subjectsFiles)
        myStr = load([groupDir subjectsFiles(s).name]);
        tempName = fieldnames(myStr);
        subjectStruct = getfield(myStr, tempName{1});
        % loop on all measures
        measures = fieldnames(subjectStruct);
        for m = 1 : 15 % there are 15 non empty measures
            measuresMatrix(s,:,m) = getfield(subjectStruct, measures{m});
        end
    end
    save([filesDirectory groups{g} '_measurements Matrix.mat'], 'measuresMatrix');
    clear measuresMatrix;
end
% save the results to excel (donors)
for g = 1 : 1
    load([filesDirectory groups{g} '_measurements Matrix.mat']);
    oF = ['measuresMatrix_' groups{g} '.xlsx'];
    for m = 2 : 15
        xlswrite(oF, 1:30, m-1, 'B1');
        xlswrite(oF, measuresMatrix(1,:,1)', m-1, 'A2');
        xlswrite(oF, measuresMatrix(:,:,m)', m-1, 'B2');
    end
    xlsheets(measures(2:15)', [pwd '/' oF]);
end





