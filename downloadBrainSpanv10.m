%%% 15 Oct. 2013
%%% donload the BrainSpan v10 data

%% define the data directory
dataDir = 'C:/Users/amahfouz/Documents/MATLAB/Data/BrainSpanV10_15Oct2013/'
donorDataDirectory = 'C:/Users/amahfouz/Documents/MATLAB/Data/GenomicConnectivityNetworks/files/genes_matrix_csv_14Oct2013/';
%% define the url of the data (HTTP)
dataURL = 'http://download.alleninstitute.org/brainspan/RNASeq_Gencode_v10/';

%% read the subject's names
[num txt] = xlsread([dataDir 'BrainSpan_RNASeq_Specimen_IDs.xlsx']);
donorNames = txt(2:end,1);
folderNames = txt(2:end,2);

% get the donors data and the structures data
load([donorDataDirectory 'donor.mat']);
load([donorDataDirectory 'structure.mat']);

for f = 1 : length(folderNames)
    % create a donor folder for the output
    outDir = [dataDir folderNames{f} '/'];
    if ~exist(outDir)
        mkdir(outDir);
    end
    % get the indicies of the current donor samples
    donorIND = find(ismember(donor.name, donorNames{f}) == 1);
    % get the list of structures available for the current donor
    donorStructures = structure.acronym(donorIND);
    % retrieve the donors data
    for strs = 1 : length(donorStructures)
        currFolderName = folderNames{f};
        % for files with "R"
        try
            currFileName = [currFolderName(1:3) '_' currFolderName(4:end) '_' donorStructures{strs} '_R_RNASeq_Gene.txt'];
            currURL = [dataURL folderNames{f} '/' currFileName];
            tempData = urlread(currURL);
            clear tempdata;
        % for files with "L"
        catch err
            currFileName = [currFolderName(1:3) '_' currFolderName(4:end) '_' donorStructures{strs} '_L_RNASeq_Gene.txt'];
            currURL = [dataURL folderNames{f} '/' currFileName];
        end
        % define the output file
        outFile = [outDir currFileName];
        try 
            currFile = urlwrite(currURL, outFile);
        catch err
            display([currFileName ' > not found'])
        end
    end
end

