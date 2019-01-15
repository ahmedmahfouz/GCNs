%%% Analyze gene pairs corr change over time

clear all;

D = 'ASDs';

%%% create a geneList correlation matrix
% load(['files/' D 'GeneIDs_5RPKM.mat']);
% for i = 1 : 7
%     load(['results/geneCorr/Avg(log2+epsilon)/ageGroup' num2str(i) '_corrMat.mat']);
%     glCorr(:,:,i) = corrMat(dIDs_T, dIDs_T);
%     clear corrMat; 
% end
% save(['results\corrChange\' D 'Corr.mat'], 'glCorr');

%%% load the gene list correlation matrix
load(['results\corrChange\' D 'Corr.mat']);
t = ones(size(glCorr,1), size(glCorr,2)); 
t = triu(t); t = t * 10000;
p = find(t ~= 10000); clear t;
for i = 1 : 7    
    clear corrDist;
    corrDist(:,:) = glCorr(:,:,i);
    modCorr(:,i) = corrDist(p);
end

%%% (1) remove gene-pairs that never show corr >= T or <= -T
temp = modCorr;
T = 0.8;
for k = 1 : size(temp,1)
    x1 = find(temp(k,:) >= T);
    spC(k) = length(x1);
    x2 = find(temp(k,:) <= (-1*T));
    snC(k) = length(x2);
end
spcID = find(spC > 0);
sncID = find(snC > 0);
sConn = temp([spcID sncID], :);
newP = p([spcID sncID]);

%%% (2) plot all gene-pair connections
% f = figure('Visible', 'off');
% plot(modCorr'), title(['allConnections_' D])
% saveas(f, ['results\corrChange\allConnections_' D '.fig']);
% saveas(f, ['results\corrChange\allConnections_' D '.jpg']);

%%% (3) plot gene-pairs with strong connections
% f = figure;%('Visible', 'off');
% plot(sConn'), title(['strongConnections_' D])
% legend
% saveas(f, ['results\corrChange\strongConnections' D '.fig']);
% saveas(f, ['results\corrChange\strongConnections' D '.jpg']);

%%% (4) cluster gene-pairs based on their correlation pattern
% cT = 0.85; % cutoff threshold
% linkType = 'average'; % linkage type
% corrType = 'correlation';
% corrDist = pdist(sConn, corrType);
% clusterTree = linkage(corrDist, linkType);
% clusters = cluster(clusterTree, 'cutoff', cT, 'criterion', 'distance');
% % f1 = figure;
% % [H, Leafs, P] = dendrogram(clusterTree, 1167, 'colorthreshold', cT, 'orientation', 'left');
% C = clustergram(sConn, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', ...
%     corrType, 'Linkage', linkType, 'Dendrogram', cT, 'Colormap', 'redbluecmap');
% % h = plot(C);
% % saveas(h, ['results\corrChange\gene-pairs clustering_' D '.fig']);
% % saveas(h, ['results\corrChange\gene-pairs clustering_' D '.jpg']);
% % save(['results\corrChange\clusters_' D '.mat'], 'clusters');

%%% (5) return names of genes in each cluster
% load(['results\corrChange\clusters_' D '.mat']);
% [num txt] = xlsread(['Data\Donors\Book1.xls']);
% dgNames = txt(2:end,1);
% dgEn = txt(2:end,2);
% clear num; clear txt;
% gnf = ['results\corrChange\clusters_' D '_Genes.xls'];
% for i = 1 : length(unique(clusters))
%     clustInd = find(clusters == i);
%     [rG cG] = ind2sub([size(glCorr,1), size(glCorr,2)], newP(clustInd));
%     clustGeneInd = union(unique(rG), unique(cG));
%     clustGeneNames{i} = dgNames(clustGeneInd);
%     clustGeneEn{i} = dgEn(clustGeneInd);
%     
%     colName = ['cluster' num2str(i)];
%     xlswrite(gnf, {colName}, i, 'A1');
%     xlswrite(gnf, clustGeneNames{i}, i, 'A2');
%     xlswrite(gnf, clustGeneEn{i}, i, 'B2');
% end

%%% (6) corr-change pattern of each cluster
load(['results\corrChange\clusters_' D '.mat']);
for i = 3 : length(unique(clusters))
    clustInd = find(clusters == i);
    cC = mean(sConn(clustInd,:), 1);
    f1 = figure;%('Visible', 'Off');
    hold on
    plot(sConn(clustInd,:)', 'linewidth', 0.25, 'Color', [0.8 0.8 0.8]);
    plot(cC, 'linewidth', 2); grid on
    title(['Cluster#' num2str(i) ' - ' num2str(length(clustInd)) ' Genes'], 'fontweight', 'bold');
    set(gca, 'XTick', 0:8);
    xlabel('Age Stages', 'fontweight', 'bold');
    ylabel('Median. Exp.', 'fontweight', 'bold');
%     axis([0.8 7.2 0 1])
    hold off
    saveas(f1, ['results\corrChange\clusters_' D '_cluster' num2str(i) '_pc1.jpg']);
    saveas(f1, ['results\corrChange\clusters_' D '_cluster' num2str(i) '_pc1.fig']);
    close(f1);
end


%%% (7) cell-type enrichment
load(['results\corrChange\clusters_' D '.mat']);

load('files/genesStatus_5RPKM.mat');
load('files/entrezIDs.mat');
load('files/gNames.mat');
gIND = find(genesStatus_5RPKM == 1);
gNames_5RPKM = gNames(gIND);
entrezIDs_5RPKM = entrezIDs(gIND);
[num txt] = xlsread(['Data\Donors\' D '.xls']);
dgNames = txt(:,1);
clear num; clear txt;
for i = 1 : length(unique(clusters))
    clustInd = find(clusters == i);
    [rG cG] = ind2sub([size(glCorr,1), size(glCorr,2)], newP(clustInd));
    clustGeneInd = union(unique(rG), unique(cG));
    [c, ia, ib] = intersect(dgNames(clustGeneInd), gNames_5RPKM);
    clustGeneEIDs{i} = entrezIDs_5RPKM(ib);
end
clear genesStatus_5RPKM; clear entrezIDs; clear gNames; clear gNames_5RPKM;
clear entrezIDs_5RPKM;
ctNameList = {'Neurons', 'Oligodendrytes', 'Astrocytes'};
outF = ['results\corrChange\'];
enrichAnalyze_cellType_corrChange(clustGeneEIDs, ctNameList, D, outF);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (3) SZ
% load('results\corrChange\szCorr.mat');
% t = ones(size(szCorr,1), size(szCorr,2)); 
% t = triu(t); t = t * 10000;
% p = find(t == 10000); clear t;
% for i = 1 : 7    
%     
%     clear corrDist;
%     corrDist(:,:) = szCorr(:,:,i);
%     corrDist(p) = [];
%     modCorr(:,i) = corrDist';
%     
% end
% % f = figure('Visible', 'off');
% % plot(modCorr'), title('allConnections_SZ')
% % saveas(f, 'results\corrChange\allConnections_SZ.fig');
% % saveas(f, 'results\corrChange\allConnections_SZ.jpg');
% 
% % toHist = reshape(modCorr, 1, size(modCorr,1)*size(modCorr,2));
% % figure, hist(toHist,1000);
% 
% % corrVar = var(modCorr');
% % figure, hist(corrVar,1000);
% 
% temp = modCorr;
% T = 0.8;
% for k = 1 : size(temp,1)
%     x1 = find(temp(k,:) >= T);
%     spC(k) = length(x1);
%     x2 = find(temp(k,:) <= (-1*T));
%     snC(k) = length(x2);
% end
% spcID = find(spC > 0);
% sncID = find(snC > 0);
% sConn = temp([spcID sncID], :);
% % f = figure;%('Visible', 'off');
% % plot(sConn'), title('strongConnections_SZ')
% % legend
% % saveas(f, 'results\corrChange\strongConnectionsSZ.fig');
% % saveas(f, 'results\corrChange\strongConnectionsSZ.jpg');
% C = clustergram(sConn, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', 'spearman', ...
%         'Linkage', 'average');
