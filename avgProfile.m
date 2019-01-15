%%% 27 Nov 2013
%%% A function to plot the average profile of given clusters

function avgProfile(clusters, data, rowLabels, colLables, colGroups)

for period = 1 : length(colGroups)
    groupDonors = colGroups{period};
    groupedData(:,period) = mean(data(:,groupDonors)')';
end

COLOR = jet(4);
figure, hold on
for c = 1 : length(unique(clusters))
    currClusterData = groupedData(find(clusters == c), :);
    currClusterRows = rowLabels(find(clusters == c));
    
    L{c} = ['cluster ' num2str(c) ' - ' num2str(size(currClusterData,1)) ' region pairs'];
    plot(mean(currClusterData), 'linewidth', 3, 'color', COLOR(c,:))
%     title(['cluster ' num2str(c) ' - ' num2str(size(currClusterData,1)) ' region pairs'],...
%         'FontWeight', 'bold', 'FontSize', 15)
end
set(gca, 'XTickLabels', colLables, 'XTick', [1:length(colLables)], 'xlim', [0 length(colLables)+1], 'FontWeight', 'bold', 'FontSize', 15)
% rotateXLabels(gca(), 45 );
ylim([0.5 1]);
ylabel('Correlation', 'FontWeight', 'bold', 'FontSize', 15);
xlabel('Developmental Period', 'FontWeight', 'bold', 'FontSize', 15);
legend(L), grid on;
hold off

