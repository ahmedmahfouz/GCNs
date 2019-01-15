function subGeneMat(dName, tRPKM, iExpMat, gStatusFile)

load('files\gNames.mat');
load('files\strucs.mat');
load('files\entrezIDs.mat');
load('files\enIDs.mat');

load(iExpMat);
expMat = donorsExpMat_5RPKM;
clear donorsExpMat_5RPKM;

load(gStatusFile);
statArr = genesStatus_5RPKM;
clear genesStatus_5RPKM;

gIND = find(statArr == 1);
gN_T = gNames(gIND);
eIDs_T = entrezIDs(gIND);
enIDs_T = enIDs(gIND);
clear gNames; clear entrezIDs; clear enIDs;

% read disease genes entrezIDs & names
[num txt] = xlsread('files\DataFiles\GeneLists.xls', dName);
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
dIDs_T = union(gInd, gN);

% temp = union(deID_Ind, dNames_Ind);
dNames_T = gN_T(dIDs_T);
deIDs_T = eIDs_T(dIDs_T);
denIDs_T = enIDs_T(dIDs_T);

% fname = ['Data\Donors\' dName '_TEMP.xls'];
% xlswrite(fname, dNames_T, 1, 'A2');
% xlswrite(fname, denIDs_T, 1, 'B2');

% create & save ASD geneList
% list of included structures + the entrezID header
strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24, 26];
for i = 1 : length(strucIndInc)
    S{i+1} = strucs{strucIndInc(i)};
end
S{1} = 'Gene';

fname = ['Data\Donors\' dName '_TEMP.xls'];

for i = 1 : size(expMat, 3)
        
    xlswrite(fname, log2(expMat(dIDs_T,:,i)+2), i, 'B2');
    xlswrite(fname, deIDs_T, i, 'A2');
    xlswrite(fname, S, i, 'A1');
    
end

dGeneMat = expMat(dIDs_T, :, :);
save(['files\' dName 'GeneIDs_' tRPKM '.mat'], 'dIDs_T');
save(['files\' dName 'GeneMat_' tRPKM '.mat'], 'dGeneMat');
%--------------------------------------------------------------------------


