%%% 22 April 2013
%%% get indicies of a geneSet within the BrainSpan Atlas

function geneSetInd = geneSetInBrainSpan(geneListFile,geneSet,gene)

% read disease genes entrezIDs & names
[num txt] = xlsread(geneListFile, geneSet);

if ~isempty(num)
    deIDs = num;
    % find the indecies of genes based on gene entrez_id
    entrez_idIND = find(ismember(gene.entrezID,deIDs) == 1);
else 
    entrez_idIND = [];
end

deNames = txt(2:end, 1);
% find the indecies of genes based on gene symbols
namesIND = find(ismember(gene.symbol,deNames) == 1);

% return the union of the two sets
geneSetInd = union(namesIND, entrez_idIND);


