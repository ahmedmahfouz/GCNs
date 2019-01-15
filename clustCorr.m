%%% 25 Nov 2013
%%% Cluster structural correlations 

function clusters = clustCorr(gcn, expMat, strPairs, cT, linkType)

corrDist = pdist(gcn);
clusterTree = linkage(corrDist, linkType);
clusters = cluster(clusterTree, 'cutoff', cT, 'criterion', 'distance');
figure; hold on
subplot(1,3,3)
[H, Leafs, P] = dendrogram(clusterTree, 0, 'colorthreshold', cT, ...
    'orientation', 'right'); set(H,'LineWidth',2); axis off
subplot(1,3,[1,2]), imagesc(gcn(P(end:-1:1),:)), colormap('redbluecmap')
set(gca, 'XTickLabels', expMat.ages, 'XTick', [1:length(expMat.ages)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45);
set(gca, 'YTickLabels', strPairs(P(end:-1:1)), 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
hold off