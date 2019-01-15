%%% 2 July 2013
%%% combine network measures in one matrix/file

filesDirectory = 'results/NetworkMeasures_donors/';
% % loop on all groups (All, SZ, ASD and ND)
groups = {'All', 'ASD', 'ND', 'SZ'};
for g = 1 : 1
    groupDir = [filesDirectory groups{g} '/'];
    % loop on all the subjects
    subjectsFiles = dir([groupDir, '\*.mat']);
    for s = 1 : length(subjectsFiles)
        myStr = load([groupDir subjectsFiles(s).name]);
        tempName = fieldnames(myStr);
        subjectStruct = getfield(myStr, tempName{1});
        % loop on all measures
        measures = fieldnames(subjectStruct);
        for m = 1 : 15 % there are 15 non empty measures
            measuresMatrix(s,:,m) = getfield(subjectStruct, measures{m});
        end
    end
    save([filesDirectory groups{g} '_measurements Matrix.mat'], 'measuresMatrix');
    clear measuresMatrix;
end

% save the results to excel (donors)
for g = 1 : 1
    load([filesDirectory groups{g} '_measurements Matrix.mat']);
    oF = ['measuresMatrix_' groups{g} '.xlsx'];
    for m = 2 : 15
        xlswrite(oF, 1:30, m-1, 'B1');
        xlswrite(oF, measuresMatrix(1,:,1)', m-1, 'A2');
        xlswrite(oF, measuresMatrix(:,:,m)', m-1, 'B2');
    end
    xlsheets(measures(2:15)', [pwd '/' oF]);
end


