%%% A function to Normalize the RPKM expression values across different
%%% samples - see equation (5) in [Anders S & Huber W, Genome Biology 2010] 

function expMat = normalizeExpMat(expMat)

temp = reshape(expMat, size(expMat,1), size(expMat,2)*size(expMat,3));
temp = temp .^ (1 / size(temp,2));
denom = (prod(temp,2)); 

z = find(denom ~= 0);

% count = 0;
for i = 1 : size(expMat,2)
    for j = 1 : size(expMat,3)
        
        s = median(expMat(z,i,j) ./ denom(z));
        expMat(:,i,j) = expMat(:,i,j) / s;
%         count = count + 1;
%         toPlot(count)=s;
        
    end
end

% figure, plot(toPlot);



