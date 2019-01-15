%%% Ahmed: edited to stop saving the results and plotting the graphs

%19 Sept 2010 - Petra Vertes
%This Function uses a non-standard thresholding method, which avoids
%problems with disconnected networks at low cost. It creates the MST of the 
%network and then grows the network edge by edge according to weight in Correlation matrix 
%CAREFUL: Algorithm constructs undirected binary nets 


%arguments of the function:
%Co: the correlation matrix that you want to study
%ext: a string used to name the output file which stores the results
%step: the number of edges to add at each 'step' before recalculating the
%network measures. This will set the density of datapoints on the curves in
%the final result
%c1: a string used to specify the color of the curves plotted at the end

%example use:
%NetworkMeasures(A,'myresults',100,'r');


function [s AdjMat] = NetworkMeasures(Co, ext, step, c1)

%Declare the variables to store all measures that will be used
s.cost=[]; s.k=[]; s.a=[]; s.arand=[]; s.M=[]; s.Mrand=[];
s.C=[]; s.Crand=[]; s.L=[]; s.Lrand=[]; s.Sigma=[]; 
s.E=[]; s.Erand=[]; s.CE=[]; s.CErand=[];
s.Diam=[]; s.Diamrand=[]; s.Bass=[]; s.Bassrand=[];
A=[]; R=[];


%Take absolute value of Correlations and set diagonal to ones:
n=size(Co,1);
Co=abs(Co);      %%%%%%%%%%%%%%%%%%%%%%%%% TAKING ABS VALUE 
Co(1:n+1:n*n)=1; %%%%%%%%%%%%%%%%%%%%%%%%% ONES ON DIAGONAL 

%Create MST (the minimum spanning tree of the network  
MST=kruskal_mst(sparse(sqrt(2*(1-Co))));

%Order C according to decreasing wieghts in the correlation matrix
Co=triu(Co,1);
ind = find(Co+triu(ones(n,n),1)); %%%TRICK: necessary in case there are zeros in the matrix 
Clist = Co(ind);
Cnonz = length(Clist);
[ClistSort, IX] = sort(Clist,'descend');
[row col]=ind2sub([n,n],ind(IX));
dd= length(Clist);

%Store Initial MST in the adjacency matrix A that defines the network
A=full(MST);
[i,j]=find(MST);
for m=1:length(i)
    A(i(m),j(m))= 1; %Co(i(m),j(m));  %(NOT) WEIGHTED VERSION
    A(j(m),i(m))= 1; %Co(i(m),j(m));  %(NOT) WEIGHTED VERISON
end

%find corresponding random matrix R
R=randmio_und_connected(A, 5);

%Start Growing the network according to weights in Co matrix and record Network Measures
%after each edge addition

%Initially, with just the MST: set counters and calculate cost and all measures
t=1;
enum=n-1;
g = 1;
s.cost(g)=enum/(n*(n-1));
gmeasure;

%Now add edges in correct order until all possible edges exist
while (enum < 0.95*n*(n-1)/2)
%     enum
    % if edge wasn't initially included in MST
    if A(row(t),col(t)) == 0
        %add edge
        A(row(t),col(t)) = 1; %Co(row(t),col(t)); %NOT WEIGHTED VERSION
        A(col(t),row(t)) = 1; %Co(row(t),col(t)); %NOT WEIGHTED VERSION       
        enum=enum+1;
        if mod(enum, step) == 0
            %find corresponding R matrix 
            R = randmio_und_connected(A, 5);
            %Increment counter
            g = g + 1;
            %calculate cost
            s.cost(g)=2*enum/(n*(n-1));
            if s.cost(g) >= 0.5 && s.cost(g) <= 0.51
                AdjMat = A;
            end
            %Call function that calculates all measures
            gmeasure; %%THIS FUNCTION CALCULATES THE MEASURES WE WANT
        end
    end
    t=t+1;
end


% %Transfer the structure containing the measures into a correctly named
% %variable for saving.
% eval(['Meas_' ext '= s;']);

% %Save the structure in a .mat file
% save(['Meas_' ext '.mat'], ['Meas_' ext]);

% %plot what you want
% gplot; %THIS FUNCTION PLOTS ALL THE MEASURES AS A FUNCTION OF COST
