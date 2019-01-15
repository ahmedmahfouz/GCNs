%%% 20 Nov 2013
%%% impute missing values using a custom scheme based on similar donors

function returnedMatrix = customImpute(inExpMat,devPeriods)

if size(inExpMat.expression,1) > 1
    checkingMatrix = squeeze(sum(inExpMat.expression));
else
    checkingMatrix = squeeze(inExpMat.expression);
end

for don =  1 : size(inExpMat.expression,3)
    % check if the current sample is missing
    missingSamples = find(isnan(checkingMatrix(:,don)) == 1);
    if ~isempty(missingSamples)
        % check how many donors have the same age as the current donor
        donorNumber = find(ismember(inExpMat.ages,inExpMat.ages(don)) == 1);
        if length(donorNumber) > 1 % get the average expression of donors of the same age
            inExpMat.expression(:,missingSamples,don) = nanmean(inExpMat.expression(:,missingSamples,donorNumber),3);
        else % get the average expression of donors within the same age period
            % get which age period is the donor in
            for period = 1 : length(devPeriods)
                P(period) = length(find(devPeriods{period} == don));
            end
            donorsToUse = devPeriods{find(P == 1)}; 
            inExpMat.expression(:,missingSamples,don) = nanmean(inExpMat.expression(:,missingSamples,donorsToUse),3);
        end
    end
end

returnedMatrix = inExpMat;
