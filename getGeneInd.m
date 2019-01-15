%%% 22 April 2013
%%% get indicies of a geneSet

function geneInd = getGeneInd(geneSet)

filesDirectory = 'files/';
resultsDirectory = 'results/';

load([filesDirectory 'genesStatus_5RPKM.mat']);
load([filesDirectory 'gNames.mat']);
load([filesDirectory 'strucs.mat']);
load([filesDirectory 'entrezIDs.mat']);
load([filesDirectory 'enIDs.mat']);

gIND = find(genesStatus_5RPKM == 1);
gN_T = gNames(gIND);
eIDs_T = entrezIDs(gIND);
enIDs_T = enIDs(gIND);
clear gNames; clear entrezIDs; clear enIDs;

% read disease genes entrezIDs & names
[num txt] = xlsread('files\DataFiles\GeneLists.xls', geneSet);
deIDs = num;
deNames = txt(2:end, 1);
clear num; clear txt;

c1 = 0; c2 = 0;
for i = 1 : length(deIDs)
    clear tem1; clear eID; clear temp2; clear n;
    % check with entrezID
    deID = deIDs(i);
    temp1 = find(eIDs_T == deID);
    if length(temp1) ~= 0
        c1 = c1 + 1;
        gInd(c1) = temp1(1);
        deID_Ind(c1) = i;
    end
    % check with gene name
    dn = deNames(i); 
    temp2 = strmatch(dn, gN_T, 'exact');
    if length(temp2) ~= 0
        c2 = c2 + 1;
        gN(c2) = temp2(1);
        dNames_Ind(c2) = i;
    end
end

% check the differences between the entrezIDs and gene names
geneInd = union(gInd, gN);





