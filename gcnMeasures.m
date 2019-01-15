%%% 15 Oct. 2013
%%% given a correlation matrix (GCN), calculate topological network
%%% measures

function gcnMeasures(dataDirectory,resultsDirectory,completeDonor,donorGCN,oF)

% check if the results directory exists, and create it if not
if ~exist(resultsDirectory, 'dir')
    mkdir(resultsDirectory);
end
% load the donors' info
load([dataDirectory 'donor.mat']);
% for each donor's GCN
for sub = 1 : size(donorGCN,3)
    subIDX = find(strcmpi(donor.name,completeDonor{sub})==1);
    donorID(sub) = unique(donor.id(subIDX));
    donorAge(sub) = unique(donor.age(subIDX));
    donorGender(sub) = unique(donor.gender(subIDX));
    donorName(sub) = unique(donor.name(subIDX));
    outFile = 'X';
    currGCN = double(single(donorGCN(:,:,sub)));
    [meas(sub) adjMat(:,:,sub)] = NetworkMeasures(currGCN, outFile, 1, 'b');
    % combine measures in one big matrix (rows: donorBrains; columns:costs; depth:measures)
    measures = fieldnames(meas(sub));
    for m = 1 : 15 % there are 15 non empty measures
        measuresMatrix(sub,:,m) = getfield(meas(sub), measures{m});
    end
end
% save the output to an excel sheet
for m = 2 : 15
    % write the column headers
    xlswrite(oF, donorName, m-1, 'B1');
    xlswrite(oF, donorID, m-1, 'B2');
    xlswrite(oF, donorAge, m-1, 'B3');
    xlswrite(oF, donorGender, m-1, 'B4');
    xlswrite(oF, {'cost'}, m-1, 'A4');
    % write the cost values to the first column
    xlswrite(oF, measuresMatrix(1,:,1)', m-1, 'A5');
    % write the current measure values 
    xlswrite(oF, measuresMatrix(:,:,m)', m-1, 'B5');
end
xlsheets(measures(2:15)', oF);



