%%% 14 Oct. 2013
%%% read the original data of the BrainSpan Atlas (RNA-SEQ) and arrange in
%%% usable matlab files

function arrangeBrainSpan(dataDir);

%% read column lables from file (columns_metadata.csv)
if ~exist([dataDir 'donor.mat'], 'file') && ~exist([dataDir 'structure.mat'], 'file')
    [num txt] = xlsread([dataDir 'columns_metadata.xlsx']);
    donor.id = num(:,2);
    donor.name = txt(2:end,3);
    donor.age = txt(2:end,4);
    donor.gender = txt(2:end,5);
    structure.id = num(:,6);
    structure.acronym = txt(2:end,7);
    structure.name = txt(2:end,8);
    save([dataDir 'donor.mat'], 'donor');
    save([dataDir 'structure.mat'], 'structure');
else
    display('"donor" & "structure" files exist already...')
end

%% read the row labels from file (rows_metadata.csv)
if ~exist([dataDir 'gene.mat'], 'file')
    [num txt] = xlsread([dataDir 'rows_metadata.xlsx']);
    gene.ensemblID = txt(2:end,3);
    gene.symbol = txt(2:end,4);
    gene.entrezID = num(:,5);
    save([dataDir 'gene.mat'], 'gene');
else
    display('"gene" file exist already...')
end

%% construct the original expression matrix (rows:Genes; columns:samples)
if ~exist([dataDir 'origExpMatrix.mat'], 'file')
    origExpMatrix = csvread([dataDir 'expression_matrix.csv']);
    origExpMatrix(:,1) = [];%remove the first column (row numbers)
    save([dataDir 'origExpMatrix.mat'], 'origExpMatrix');
else
    display('"expression matrix" file exist already...')
end




