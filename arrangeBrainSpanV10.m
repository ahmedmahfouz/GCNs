%%% 16 Oct. 2013
%%% read the original data of the BrainSpan Atlas V10 (RNA-SEQ) and arrange in
%%% usable matlab files

function arrangeBrainSpanV10(dataDirV10,dataDirV3);

%% load the BrainSpan V3 donor and structure info
load([dataDirV3 'donor.mat']);
load([dataDirV3 'structure.mat']);

%% read the subject's names
[num txt] = xlsread([dataDirV10 'BrainSpan_RNASeq_Specimen_IDs.xlsx']);
donorNames = txt(2:end,1);
externalNames = txt(2:end,2);

%% read the donors files (separate file for each structure)
if ~exist([dataDirV10 'donorV10.mat'], 'file') && ~exist([dataDirV10 'structureV10.mat'], 'file')
    sampleCount = 0;
    donorFiles = dir(dataDirV10);
    for f = 3 : length(donorFiles)-2
        f
        if isdir([dataDirV10 '/' donorFiles(f+2).name])
            donorStrFiles = dir([dataDirV10 '/' donorFiles(f+2).name]);
            %find the index of the current donor in the file mapping
            %epcimen to external names
            ind = find(strcmpi(externalNames, donorFiles(f+2).name) == 1);
            %find the index of the current donor in the donor list (V3)
            donorIND = unique(find(strcmpi(donor.name,donorNames(ind)) == 1));
            if isempty(donorIND)
                display('donor not found')
                break 
            end
            %fill in the donor structure of the V10 data
            donorV10.name = donor.name(donorIND);
            donorV10.id = donor.id(donorIND);
            donorV10.age = donor.age(donorIND);
            donorV10.gender = donor.gender(donorIND);
            clear donorIND; clear ind;
            %loop on the files (structures) of the current donor
            for strFile = 1 : length(donorStrFiles)-2
                tempStrName = donorStrFiles(strFile+2).name;
                tempPos = strfind(tempStrName,'_');
                currStrName = tempStrName(tempPos(2)+1:tempPos(3)-1);
                % read the tab delimited file of the current structure
                tempS = tdfread([dataDirV10 '/' donorFiles(f+2).name '/' donorStrFiles(strFile+2).name]);
                % retrieve gene data from one file (gene order is the same in all files)
                if ~exist('gene', 'var')
                    % entrez_id
                    geneV10.entrezID = cellstr(tempS.EntrezID);
                    % the gencode column has the ensemble_id|gene_symbol (separator position is 16)
                    gencode = tempS.GencodeID;
                    geneV10.ensebleID = cellstr(gencode(:,1:15));
                    geneV10.symbol = cellstr(gencode(:,17:end));
                    save([dataDirV10 'geneV10.mat'], 'geneV10');
                    clear gencode;
                end
                sampleCount = sampleCount + 1;
                structureV10.acronym = currStrName;
                origExpMatrixV10(:,sampleCount) = tempS.RPKM;
                clear tempS; clear currStrName; clear tempStrName; clear tempPos;
            end
            clear donorStrFiles;
        end
        clear donorStrFiles;
    end
else
    display('"donor" & "structure" files exist already...')
end




