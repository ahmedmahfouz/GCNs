%%% 14 Oct. 2013
%%% A pipeline to generate Genomic Connectivity Networks (GCNs) based on the
%%% BrainSpan Atlas (developing human brain atlas)

%% define some variables
atlasVersion = 3; %either 3 or 10
dataDirectory = ['C:/Users/amahfouz/Documents/MATLAB/Data/GenomicConnectivityNetworks/files/genes_matrix_csv_genCodeV' ...
    num2str(atlasVersion) '/'];
expThresh = 5;% define the expression threshold (RPKM)
varThresh = 1;% define the variance threshold 
corrType = 'Spearman';% define the correlation type
% define the file with the gene lists (ASD,SZ,ND)
geneListFile = 'C:/Users/amahfouz/Documents/MATLAB/Data/GenomicConnectivityNetworks/GeneLists.xls';
% % ids of donors with no missing data (based on the selected 16 structures)
% completeDonor = [12837,12365,12288,12296,12890,12297,12830,12979,12298,12841,...
%     12831,12832,13057,12300,12290,12302,12303,12291];
% extract only the 16 structures of interest
structuresToInclude = {'AMY','HIP','STR','MD','CBC','DFC','VFC','MFC','OFC',...
    'IPC','STC','ITC','A1C','M1C','S1C','V1C'};
% extract donors with enough samples
if atlasVersion == 3
    donorsToInclude = {'H376.IIIB.50','H376.IIIB.51','H376.IIIB.52','H376.IIIB.53',...
        'H376.IV.50','H376.IV.51','H376.IV.53','H376.IV.54','H376.IX.50','H376.IX.51',...
        'H376.IX.52','H376.VI.50','H376.VI.52','H376.VII.50','H376.VIII.50',...
        'H376.VIII.51','H376.VIII.52','H376.VIII.53','H376.VIII.54','H376.X.50',...
        'H376.X.51','H376.X.52','H376.X.53','H376.XI.60','H376.XI.50','H376.XI.52','H376.XI.53',...
        'H376.XI.54','H376.XI.55','H376.XI.56'};
    devPeriods = {[1, 2, 3 ,4], [5, 6, 7, 8], [9, 10, 11], [12, 13, 14, 15], ...
        [16, 17, 18, 19, 20], [21, 22, 23, 24], [25, 26, 27, 28, 29, 30]};
elseif atlasVersion == 10
    donorsToInclude = {'H376.IIIB.50','H376.IIIB.51','H376.IIIB.52','H376.IIIB.53',...
        'H376.IV.53','H376.IV.54','H376.IV.50','H376.V.53',...
        'H376.VI.50','H376.VI.52','H376.VIII.51',...
        'H376.VIII.53','H376.VIII.54','H376.VIII.52','H376.VIII.50',...
        'H376.IX.51','H376.IX.52','H376.IX.50','H376.X.51',...
        'H376.X.50','H376.X.53','H376.X.52','H376.XI.60',...
        'H376.XI.50','H376.XI.52','H376.XI.53','H376.XI.54','H376.XI.56'};
    devPeriods = {[1, 2, 3 ,4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15], ...
        [16, 17, 18, 19], [20, 21, 22, 23], [24, 25, 26, 27, 28]};
end
devPeriodsNames = {'Early 2nd Trimester', 'Late 2nd trimester', 'Infancy', ...
    'Early Childhood', 'Late Childhood', 'Adolescence', 'Adulthood'};

devPeriods3 = {1:8, 9:20, 21:30};
devPeriodsNames3 = {'Prenatal', 'Childhood', 'Adulthood'};

%% Read original data from the BrainSpan Atlas (V3) and arrange it matlab files
arrangeBrainSpan(dataDirectory);
load([dataDirectory 'donor.mat']);
load([dataDirectory 'structure.mat']);
load([dataDirectory 'gene.mat']);
load([dataDirectory 'origExpMatrix.mat']);

%% filter the data (5 RPKM in at least one sample)
% select the donors and structures of interest
donorsSelect = find(ismember(donor.name, donorsToInclude));
structuresSelect = find(ismember(structure.acronym, structuresToInclude));
includedSamples = intersect(donorsSelect,structuresSelect);
data = origExpMatrix(:,includedSamples);
% construct a donor expression matrix
donorExpMat.expression = NaN(size(data,1),length(structuresToInclude),length(donorsToInclude));
for d = 1 : length(donorsToInclude)
    % find the ids of the current donors
    donorIDs = find(ismember(donor.name,donorsToInclude{d}) == 1);
    % find the ids of the selected structures corresponding to the current
    % donor
    allDonorStructures = structure.acronym(donorIDs);
    % rearrange donor structures to match the order of the structuresOfInteres
    [tempSTR ia ib] = intersect(structuresToInclude,allDonorStructures);
    [sortedIA sortingInd] = sort(ia);
    donorStrIDs = donorIDs(ib(sortingInd));
    % create a donor-specific expression matrix
    donorExpMat.expression(:,sortedIA,d) = origExpMatrix(:, donorStrIDs);
    % retrieve a list of ages corresponding to the list of donors
    donorExpMat.ages(d) = donor.age(donorIDs(1));
    clear donorIDs; clear allDonorStructures; clear tempSTR;
    clear ia; clear ib; clear sortedIA; clear sortingInd; clear donorStrIDs;
end
% impute missing values
donorExpMat = customImpute(donorExpMat,devPeriods);
% remove genes that have no expression above 5 RPKM in one sample at least
data = reshape(donorExpMat.expression, size(donorExpMat.expression,1), size(donorExpMat.expression,2)*size(donorExpMat.expression,3));
[row col] = find(data >=5);
row = unique(row);

% data_max = max(data');
% expressingGenes = find(data_max ~= 0);
% figure, hist(log2(data_max(expressingGenes)), 100); grid on;
% title('Distribution of the Maximum Expression per Gene', 'FontWeight', 'bold', 'FontSize', 20)
% set(gca, 'xlim', [-15 20], 'ylim', [0 2500], 'FontWeight', 'bold', 'FontSize', 15)
% xlabel('log_2(maximum expression)', 'FontWeight', 'bold', 'FontSize', 20);
% ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)
% get the indicies of expressing genes (with no outliers) and clear
% variables
expressingGenesIDs = row;

%% filter the data (1 RPKM in 80% 0f samples)
% select the donors and structures of interest
donorsSelect = find(ismember(donor.name, donorsToInclude));
structuresSelect = find(ismember(structure.acronym, structuresToInclude));
includedSamples = intersect(donorsSelect,structuresSelect);
data = origExpMatrix(:,includedSamples);
% remove genes that have no expression at all (zero everywhere)
data_max = max(data');
expressingGenes = find(data_max ~= 0);
figure, hist(log2(data_max(expressingGenes)), 100); grid on;
title('Distribution of the Maximum Expression per Gene', 'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'xlim', [-15 20], 'ylim', [0 2500], 'FontWeight', 'bold', 'FontSize', 15)
xlabel('log_2(maximum expression)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)
% remove outliers: genes that have an extremely high expression (>2000 RPKM)
[sortedMaxVals sortingInd] = sort(data_max);
outliers = sortingInd(sortedMaxVals > 2000);
outliersInd = find(ismember(expressingGenes,outliers));
expressingGenes_noOutliers = expressingGenes;
expressingGenes_noOutliers(outliersInd) = [];
figure, hist(log2(data_max(expressingGenes_noOutliers)), 100); grid on;
title({'Distribution of the Maximum Expression per Gene', ...
    ['oultiers removed (' num2str(length(outliersInd)) ' genes with maximum expression > 2000)']}, ...
    'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'xlim', [-15 20], 'ylim', [0 2500], 'FontWeight', 'bold', 'FontSize', 15)
xlabel('log_2(maximum expression)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)
% keep only genes that have expression of 1 RPKM in 80% of the samples
data_thresh = zeros(size(data));
data_thresh(data>1) = 1;
notinManySamples = find(sum(data_thresh') < (0.8*size(data,2))); 
notInManySamples_ind = find(ismember(expressingGenes_noOutliers,notinManySamples));
expressingGenes_noOutliers_inManySamples =  expressingGenes_noOutliers;
expressingGenes_noOutliers_inManySamples(notInManySamples_ind) = [];
figure, hist(log2(data_max(expressingGenes_noOutliers_inManySamples)), 100); grid on;
title({'Distribution of the Maximum Expression per Gene', ...
    ['genes with expression of 1 RPKM in 80% of samples (' num2str(length(expressingGenes_noOutliers_inManySamples)) ' genes)']}, ...
    'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'xlim', [-15 20], 'FontWeight', 'bold', 'FontSize', 15)
xlabel('log_2(maximum expression)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)
% get the indicies of expressing genes (with no outliers) and clear
% variables
expressingGenesIDs = expressingGenes_noOutliers_inManySamples;
clear donorsSelect; clear structuresSelect;
clear data; clear data_max; clear expressingGenes; clear sortedMaxVals;
clear sortingInd; clear outliers; clear outliersInd; clear expressingGenes_noOutliers;
clear data_thresh; clear notinManySamples; clear notInManySamples_ind;
clear expressingGenes_noOutliers_inManySamples;

%% Look for missing data
uniqueDonors = unique(donor.name);
dataAnalysisMat = zeros(length(uniqueDonors), length(structuresToInclude));
for d = 1 : length(uniqueDonors)
    donorIDs = find(ismember(donor.name,uniqueDonors(d)) == 1);
    uniqueAges(d) = donor.age(donorIDs(1));
    allDonorStructures = structure.acronym(donorIDs);
    [tempSTR ia ib] = intersect(structuresToInclude,allDonorStructures);
    dataAnalysisMat(d,ia) = 1;
end

%% analyze the expression distribution
dataArr = reshape(origExpMatrix,1,size(origExpMatrix,1)*size(origExpMatrix,2));
dataArr(dataArr == 0) =[];
figure, hist(log2(dataArr), 100), grid on
title('Gene Expression Distribution', 'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'xlim', [-20 20], 'FontWeight', 'bold', 'FontSize', 15)
xlabel('log_2(expression)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)

%% analyze the variance of genes across time
% remove the genes with no expression data (zeros every where)
[Mask1, filteredData] = genevarfilter(origExpMatrix, 'Percentile', 0);
figure, hist(log2(var(filteredData')), 100), grid on
title('Gene Variance Distribution', 'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'xlim', [-30 30], 'FontWeight', 'bold', 'FontSize', 15)
xlabel('log_2(variance)', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)
% analyze the variance of genes with an expression value of at least 'expThresh'
[expressingGenesIDs ~] = ind2sub(size(origExpMatrix), find(origExpMatrix > expThresh));
expressingGenesIDs = unique(expressingGenesIDs);
[Mask_exp, FData_exp] = genevarfilter(origExpMatrix(expressingGenesIDs,:), 'Percentile', 0);
figure, hist(log2(var(FData_exp')), 100), grid on
title(['Gene Variance Distribution of genes with expression > ' num2str(expThresh)], 'FontWeight', 'bold', 'FontSize', 20)
xlabel('log_2(variance)', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)
% remove genes with variance < 1
[Mask_var, FData_var] = genevarfilter(origExpMatrix, 'AbsValue', 1);
figure, hist(log2(var(FData_var')), 100), grid on
title('Gene Variance Distribution of genes with variance > 1', 'FontWeight', 'bold', 'FontSize', 20)
xlabel('log_2(variance)', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 20)

%% Filter genes with low variance
% if geneSet is 'All'
% [maskVar filteredExpMat] = genevarfilter(origExpMatrix, 'AbsValue', 1);
% [maskVar filteredExpMat] = genevarfilter(origExpMatrix, 'Percentile', 0);
% expressingGenesIDs = find(maskVar == 1);

%% Normalize across genes (rows)
origExpMatrix = (origExpMatrix - repmat(nanmean(origExpMatrix')', 1, size(origExpMatrix,2))) ./ repmat(nanstd(origExpMatrix')', 1, size(origExpMatrix,2));

%% Analyze a network of random genes
randomGenesNumber = 380;
[donorGCN_rand devPeriodGCN_rand strPairs donorExpMat_rand] = gcnPerDonor_random(...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    corrType,devPeriods,devPeriodsNames,randomGenesNumber);

%% Analyze one region across time
geneSet = 'All';% define the set of genes (All, ASD, SZ, ND)
[donorGCN_All donorExpMat_All] = gcnPerRegion( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

%% Analyze one gene across region pairs
geneSet = 'All';% define the set of genes (All, ASD, SZ, ND)
[geneGCN] = gcnPerGene( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

%% create genomic connectivity matrix per donor for the list of all genes and asd gene lists
geneSet = 'All';% define the set of genes (All, ASD, SZ, ND)
[donorGCN_All_sq donorGCN_All devPeriodGCN_All strPairs donorExpMat_All] = gcnPerDonor( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

geneSet = 'ASD';% define the set of genes (All, ASD, SZ, ND)
[donorGCN_ASD_sq donorGCN_ASD devPeriodGCN_ASD strPairs donorExpMat_ASD] = gcnPerDonor( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods3,devPeriodsNames3);

%% SPIE FIGURE
ageStage1_expMat = mean(donorExpMat_All.expression(:,:,1:4),3);
ageStage2_expMat = mean(donorExpMat_All.expression(:,:,5:8),3);
ageStage7_expMat = mean(donorExpMat_All.expression(:,:,25:end),3);

figure, imagesc(zscore(ageStage1_expMat,[],2)), colormap(redbluecmap)
axis off
figure, imagesc(zscore(ageStage2_expMat,[],2)), colormap(redbluecmap)
axis off
figure, imagesc(zscore(ageStage7_expMat,[],2)), colormap(redbluecmap)
axis off


%% Analyze the significance of the correlation between region pairs
nPerm = 10000;
for p = 1 : nPerm
    randomGenesNumber = 455;
    [donorGCN_rand(:,:,p) devPeriodGCN_rand(:,:,p) strPairs donorExpMat_rand] = gcnPerDonor_random(...
        donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
        structuresToInclude,donorsToInclude, ...
        corrType,devPeriods3,devPeriodsNames3,randomGenesNumber);
end

%% compare CBC correlation between ASD and random
load('devPeriodGCN_rand_V3_3periods.mat');
C6 = distinguishable_colors(6);
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
x = 100;
y = 100;
W = 1300;
H = 700;
F1 = 30;
LF = 20;
LW = 5;

for S = 1 : 5
    inputStr = structuresToInclude{S};
    K = strfind(strPairs, inputStr);
    Index = find(not(cellfun('isempty', K)));
    f1 = figure, 
%     set(f1, 'Position', [x y W H])
    hold on
    for i = 1 : 4
        currLabel = strPairs{Index(i)};
        yL = strfind(currLabel, inputStr);
        if yL == 1
            L{i} = currLabel(length(inputStr)+2:end);
        else
            L{i} = currLabel(1:strfind(currLabel, '_')-1);
        end
        cInd = find(ismember(structuresToInclude(1:5),L{i}) == 1);
        H(i) = plot(devPeriodGCN_ASD(Index(i),:)', 'LineWidth', LW, 'Color', C6(cInd,:))
        errorbar(mean(devPeriodGCN_rand(Index(i),:,:),3),std(devPeriodGCN_rand(Index(i),:,:),[],3),...
            'LineWidth', LW, 'LineStyle', '--', 'Color', C6(cInd,:));
        grid on;
    end
    
    H(i+1) = errorbar(mean(devPeriodGCN_ASD(Index(5:end),:)),...
        std(devPeriodGCN_ASD(Index(i+1:end),:)),...
        'LineWidth', LW, 'Color', C6(end,:))
    ncxRand_mean = mean(devPeriodGCN_rand(Index(i+1:end),:,:),1);
    ncxRand_std = std(devPeriodGCN_rand(Index(i+1:end),:,:),[],1);
    errorbar(mean(ncxRand_mean,3),std(ncxRand_std,[],3),...
            'LineWidth', LW, 'LineStyle', '--', 'Color', C6(end,:));
    L{i+1} = 'NCX';
    
    set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
        'xlim',[0 length(devPeriods3)+1],'FontSize',F1,'FontWeight','bold')
    ylim([0.5 1]);
    ylabel('Correlation','FontSize',F1,'FontWeight','bold');
    legend(H,L,'FontWeight','bold','FontSize',LF)
    hold off
    set(gcf,'NextPlot','add');
    axes;
    h = title(inputStr, 'FontWeight', 'bold', 'FontSize', F1);
    set(gca,'Visible','off');
    set(h,'Visible','on');
end

%% compare CBC correlation between ASD and random
% load('devPeriodGCN_rand.mat');
for S = 1 : 5
    inputStr = structuresToInclude{S};
    K = strfind(strPairs, inputStr);
    Index = find(not(cellfun('isempty', K)));
    figure, 
    for i = 1 : length(Index)
        subplot(3,5,i),hold on
        plot(devPeriodGCN_ASD(Index(i),:)', 'LineWidth', 2, 'Color', 'r')
        errorbar(mean(devPeriodGCN_rand(Index(i),:,:),3),std(devPeriodGCN_rand(Index(i),:,:),[],3),...
            'LineWidth', 2);
        grid on;
        currLabel = strPairs{Index(i)};
        yL = strfind(currLabel, inputStr);
        if yL == 1
            ylabel(currLabel(length(inputStr)+2:end), 'FontWeight', 'bold');
        else
            ylabel(currLabel(1:3), 'FontWeight', 'bold');
        end
        set(gca,'XTickLabels',devPeriods3,'XTick',[1:length(devPeriods3)],'xlim',[0 length(devPeriods3)+1])
        ylim([0.5 1]);
        hold off
    end
    set(gcf,'NextPlot','add');
    axes;
    h = title(inputStr, 'FontWeight', 'bold', 'FontSize', 20);
    set(gca,'Visible','off');
    set(h,'Visible','on');
end

%% analyze the difference between the 'All' and 'ASD' networks (donors)
% plot the structural correlations of All and ASD
figure, hold on
subplot(1,2,1), imagesc(donorGCN_All), colormap('redbluecmap')
title('Correlation Change All', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', donorExpMat_All.ages, 'XTick', [1:length(donorsToInclude)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
subplot(1,2,2), imagesc(donorGCN_ASD), colormap('redbluecmap')
title('Correlation Change ASD', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', donorExpMat_All.ages, 'XTick', [1:length(donorsToInclude)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
hold off
% plot the difference in structural correlations between All and ASD
figure, imagesc(donorGCN_All-donorGCN_ASD), colormap('redbluecmap')
title('Differences in Structural Correlations between All and ASD', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', donorExpMat_All.ages, 'XTick', [1:length(donorExpMat_All.ages)], 'FontWeight', 'bold', 'FontSize', 20)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)

%% analyze the difference between the 'All' and 'ASD' networks (developmental periods)
% plot the structural correlations of All and ASD
figure, hold on
subplot(1,2,1), imagesc(devPeriodGCN_All), colormap('redbluecmap')
title('Correlation Change All', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', devPeriodsNames, 'XTick', [1:length(devPeriodsNames)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
xlabel('Developmental Periods', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Region Pairs', 'FontWeight', 'bold', 'FontSize', 15);
subplot(1,2,2), imagesc(devPeriodGCN_ASD), colormap('redbluecmap')
title('Correlation Change ASD', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', devPeriodsNames, 'XTick', [1:length(devPeriodsNames)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
xlabel('Developmental Periods', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Region Pairs', 'FontWeight', 'bold', 'FontSize', 15);
hold off
% plot the difference in structural correlations between All and ASD
figure, imagesc(devPeriodGCN_All-devPeriodGCN_ASD), colormap('redbluecmap')
title('Differences in Structural Correlations between All and ASD', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', devPeriodsNames, 'XTick', [1:length(devPeriodsNames)], 'FontWeight', 'bold', 'FontSize', 20)
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
xlabel('Developmental Periods', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Region Pairs', 'FontWeight', 'bold', 'FontSize', 15);

%% cluster the correlation changes over time
cT = 0.4; % ASD = 0.45, All = 0.35
linkType = 'average'; % linkage type
regPairClusters = clustCorr(donorGCN_All, donorExpMat_All, strPairs, cT, linkType);

%% cluster the correlation changes over time
cT = 0.35; % ASD = 0.45, All = 0.35
linkType = 'average'; % linkage type
regPairClusters = clustCorr(donorGCN_All, donorExpMat_All, strPairs, cT, linkType);

%% cluster the differences in corr change between All and ASD
cT = 0.15; % cutoff threshold
linkType = 'average'; % linkage type
clustCorr(donorGCN_All-donorGCN_ASD, donorExpMat_ASD, strPairs, cT, linkType);

%% plot the average profiles of clusters
avgProfile(regPairClusters, donorGCN_ASD, strPairs, devPeriodsNames, devPeriods);

%% create genomic connectivity matrix per donor for the list of cell-type specific genes 
geneSet = 'Neurons';
[donorGCN_N devPeriodGCN_N strPairs donorExpMat_N] = gcnPerDonor( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

geneSet = 'Oligodendrocytes';
[donorGCN_O devPeriodGCN_O strPairs donorExpMat_O] = gcnPerDonor( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);
 
geneSet = 'Astrocytes';
[donorGCN_A devPeriodGCN_A strPairs donorExpMat_A] = gcnPerDonor( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

%% Analyze the difeerences between Neurons/Oligodendrocytes and Astrocytes networks
% plot the structural correlations of N/O/A
figure, hold on
subplot(1,3,1), imagesc(donorGCN_A), colormap('redbluecmap')
title('Correlation Change Neurons', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', donorExpMat_A.ages, 'XTick', [1:length(donorsToInclude)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
subplot(1,3,2), imagesc(donorGCN_O), colormap('redbluecmap')
title('Correlation Change Oligodendrocytes', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', donorExpMat_O.ages, 'XTick', [1:length(donorsToInclude)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
subplot(1,3,3), imagesc(donorGCN_A), colormap('redbluecmap')
title('Correlation Change Astrocytes', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabels', donorExpMat_A.ages, 'XTick', [1:length(donorsToInclude)], 'FontWeight', 'bold', 'FontSize', 10)
rotateXLabels(gca(), 45 );
set(gca, 'YTickLabels', strPairs, 'YTick', [1:length(strPairs)], 'FontWeight', 'bold', 'FontSize', 5)
hold off

%% Calculate network topological measures for each donor's GCN
% add the path of the BCT toolbox
bctToolBox = 'C:/Users/amahfouz/Documents/MATLAB/Libraries/2012-12-04 BCT';
bgl = 'C:/Users/amahfouz/Documents/MATLAB/Libraries/matlab_bgl';
addpath(bctToolBox, bgl);
% define the output directory
resultsDirectory = 'C:/Users/amahfouz/Documents/MATLAB/Results/GenomicConnectivityNetworks/DonorsGCNs_NetMeasures/';
% define the output file
geneSet = 'ASD';
oF = [resultsDirectory  'measuresMatrix' geneSet '.xlsx'];
% call the gcnMeasure function
gcnMeasures(dataDirectory,resultsDirectory,donorsToInclude,donorGCN_ASD_sq,oF);

%% create genomic connectivity matrix per donor for the list of all genes and asd gene lists
% for big 6 structures and ncx separately
geneSet = 'All';% define the set of genes (All, ASD, SZ, ND)
[donorGCN_big6 devPeriodGCN_big6 donorGCN_ncx devPeriodGCN_ncx strPairs_big6 strPairs_ncx donorExpMat] = gcnPerDonor_2networks( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

geneSet = 'ASD';% define the set of genes (All, ASD, SZ, ND)
[donorGCN_big6 devPeriodGCN_big6 donorGCN_ncx devPeriodGCN_ncx strPairs_big6 strPairs_ncx donorExpMat] = gcnPerDonor_2networks( ...
    donor,structure,gene,origExpMatrix, expressingGenesIDs, ...
    structuresToInclude,donorsToInclude, ...
    geneSet,geneListFile,corrType,devPeriods,devPeriodsNames);

%% cluster the correlations for big6 and ncx
cT = 0.5; % ASD = 0.45, All = 0.35
linkType = 'average'; % linkage type
regPairClusters = clustCorr(donorGCN_big6, donorExpMat, strPairs_big6, cT, linkType);

%% prepare the network measures figures for SPIE 2014
measuresFile = 'C:\Users\amahfouz\Documents\MATLAB\Data\GenomicConnectivityNetworks\analysis_updated.xlsx';
measures_ALL = xlsread(measuresFile,2,'B2:AE14');
measures_ASD = xlsread(measuresFile,2,'B17:AE29');
[num measures] = xlsread(measuresFile,2,'A2:A14');
clear num;
for P = 1 : length(devPeriods3)
    currDonors = devPeriods3{P};
    meanMeasures_ALL(:,P) = mean(measures_ALL(:,currDonors)');
    standardError_ALL(:,P) = (std(measures_ALL(:,currDonors)')) / sqrt(length(currDonors));
    meanMeasures_ASD(:,P) = mean(measures_ASD(:,currDonors)');
    standardError_ASD(:,P) = (std(measures_ASD(:,currDonors)')) / sqrt(length(currDonors));
    [h,p] = ttest(measures_ALL(:,currDonors)',measures_ASD(:,currDonors)',0.05,'both');
    signTest(:,P) = p';
end

%% plot Assortativity (SPIE poster)
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
x = 100;
y = 100;
W = 1300;
H = 700;

F1 = 30;
LF = 20;
LW = 5;

f1 = figure, 
set(f1, 'Position', [x y W H])
hold on
errorbar(meanMeasures_ALL(1,:),standardError_ALL(1,:),...
    'LineWidth',LW,'Color','blue');
errorbar(meanMeasures_ASD(1,:),standardError_ASD(1,:),...
    'LineWidth',LW,'Color','red');
grid on
legend({'AllGenes-GCNs', 'ASD-GCNs'},'FontWeight','bold','FontSize',LF);
ylabel({'Assortativity',''},'FontWeight','bold','FontSize',F1);
set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
    'xlim',[0 length(devPeriods3)+1],'FontWeight', 'bold', 'FontSize', F1)
ylim([-0.2 0.5]);
text(2.97,0.44,'*','FontWeight','bold','FontSize',40,'Color','black')
hold off
saveas(f1, 'Assortativity.png');
% plot Modularity
f1 = figure, 
set(f1, 'Position', [x y W H])
hold on
errorbar(meanMeasures_ALL(3,:),standardError_ALL(3,:),...
    'LineWidth',LW,'Color','blue');
errorbar(meanMeasures_ASD(3,:),standardError_ASD(3,:),...
    'LineWidth',LW,'Color','red');
grid on
legend({'AllGenes-GCNs', 'ASD-GCNs'},'FontWeight','bold','FontSize',LF);
ylabel({'Modularity',''},'FontWeight','bold','FontSize',F1);
set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
    'xlim',[0 length(devPeriods3)+1],'FontWeight', 'bold', 'FontSize', F1)
ylim([0 0.2]);
hold off
saveas(f1, 'Modularity.png');
% plot Clustering Coefficient
f1 = figure, 
set(f1, 'Position', [x y W H])
hold on
errorbar(meanMeasures_ALL(5,:),standardError_ALL(5,:),...
    'LineWidth',LW,'Color','blue');
errorbar(meanMeasures_ASD(5,:),standardError_ASD(5,:),...
    'LineWidth',LW,'Color','red');
grid on
legend({'AllGenes-GCNs', 'ASD-GCNs'},'FontWeight','bold','FontSize',LF);
ylabel({'Clustering Coefficient',''},'FontWeight','bold','FontSize',F1);
set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
    'xlim',[0 length(devPeriods3)+1],'FontWeight', 'bold', 'FontSize', F1)
ylim([0.5 0.8]);
hold off
saveas(f1, 'ClusteringCoefficient.png');
% plot Path Length
f1 = figure, 
set(f1, 'Position', [x y W H])
hold on
errorbar(meanMeasures_ALL(7,:),standardError_ALL(7,:),...
    'LineWidth',LW,'Color','blue');
errorbar(meanMeasures_ASD(7,:),standardError_ASD(7,:),...
    'LineWidth',LW,'Color','red');
grid on
legend({'AllGenes-GCNs', 'ASD-GCNs'},'FontWeight','bold','FontSize',LF);
ylabel({'Path Length',''},'FontWeight','bold','FontSize',F1);
set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
    'xlim',[0 length(devPeriods3)+1],'FontWeight', 'bold', 'FontSize', F1)
ylim([1.5 2]);
text(2.97,1.84,'*','FontWeight','bold','FontSize',40,'Color','black')
hold off
saveas(f1, 'PathLength.png');
% plot Small-Worledness
f1 = figure, 
set(f1, 'Position', [x y W H])
hold on
errorbar(meanMeasures_ALL(9,:),standardError_ALL(9,:),...
    'LineWidth',LW,'Color','blue');
errorbar(meanMeasures_ASD(9,:),standardError_ASD(9,:),...
    'LineWidth',LW,'Color','red');
grid on
legend({'AllGenes-GCNs', 'ASD-GCNs'},'FontWeight','bold','FontSize',LF);
ylabel({'Small-Worledness',''},'FontWeight','bold','FontSize',F1);
set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
    'xlim',[0 length(devPeriods3)+1],'FontWeight', 'bold', 'FontSize',F1)
ylim([0.8 1.1]);
text(2.97,1.02,'*','FontWeight','bold','FontSize',40,'Color','black')
hold off
saveas(f1, 'SmallWorldness.png');
% plot Efficiency
f1 = figure, 
set(f1, 'Position', [x y W H])
hold on
errorbar(meanMeasures_ALL(10,:),standardError_ALL(10,:),...
    'LineWidth',LW,'Color','blue');
errorbar(meanMeasures_ASD(10,:),standardError_ASD(10,:),...
    'LineWidth',LW,'Color','red');
grid on
legend({'AllGenes-GCNs', 'ASD-GCNs'},'FontWeight','bold','FontSize',LF);
ylabel({'Efficiency',''},'FontWeight','bold','FontSize',F1);
set(gca,'XTickLabels',devPeriodsNames3,'XTick',[1:length(devPeriods3)],...
    'xlim',[0 length(devPeriods3)+1],'FontWeight', 'bold', 'FontSize', F1)
ylim([0.68 0.74]);
text(2.97,0.727,'*','FontWeight','bold','FontSize',40,'Color','black')
hold off
saveas(f1, 'Efficiency.png');
