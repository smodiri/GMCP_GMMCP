function [nodes, cost, timeSpent] = gmcpsolver( ...
    netCostMat, ...
    nClusterNodes, ...
    nClusters, ...
    isFindMax, ...
    method, ...
    showDebugInfo, ...
    kCliques)

%GMCPSOVER Solves the GMCP and GMMCP problems.
% Args:
%   netCostMat: is the weights of the GMCP graph! Please keep the nodes 
%       belonging to each cluster together. 
%
%   nClusterNodes: Number of Nodes in each Cluster. It can be a single number
%       (For Graphs with the same number of nClusterNodes in each Cluster) or can
%       be a vector of numbers which shows number of members for each
%       cluster.
%
%   nClusters: Number of Clusters
%
%   isFindMax: 1 if finding Maximum, 0 if finding Minimum.
%
%   method:
%       0 if Linear Programmin
%       1 if Binary Linear Programming
%       2 (Default-Fast-Exact) if mixed Linear Programming
%
%   kCliques: Number of cliques in the solution. In case of GMCP, this
%       kCliques = 1 and for GMMCP kCliques should be selected
%       appropreately.

%   Pleaes cite:
%   1- Shayan Modiri Assari, Amir Roshan Zamir, and Mubarak Shah. "Video classification using semantic concept co-occurrences." Computer Vision and Pattern Recognition (CVPR), 2014 IEEE Conference on
%   2- Afshin Dehghan, Shayan Modiri Assari, and Mubarak Shah. "GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking." CVPR. Vol. 1. 2015.

%% Constant Variables
% A cost considered too big to be ever selected (negative of this value
% when dealing with the scores). Pleae do not make this value too big to
% avoid computational errors. You can normalize the input costs/scores
% instead.
INF_IN = 1e6;

%% Setting the default values
if ~exist('method', 'var')
    method = 2;
end

if ~exist('showDebugInfo','var')
    showDebugInfo = false;
end

if ~exist('kCliques', 'var') || isempty(kCliques)
    kCliques = 1;
end

if length(nClusterNodes) == 1
    nClusterNodes = nClusterNodes*ones(1,nClusters);
end

if kCliques > min(nClusterNodes)
    warning('Number of Cliques can not be more than number of nodes in every cluster');
    keyboard;
end

if length(nClusterNodes) ~= nClusters
    warning('Length of nClusterNodes needs to be consistant with the provided nClusters');
    keyboard;
end

%% Adding Cplex Binary Path
% Please comment this section if you are calling this function multiple
% times and call this once outside this function.
if isunix
    cplex_path = '~/Tools/IBM_Cplex/cplex/matlab';
elseif ismac
    cplex_path = '/opt/IBM/ILOG/CPLEX_Studio1262/cplex/matlab';
elseif ispc
    cplex_path = ...
        '~\Program Files\IBM\ILOG\CPLEX_Studio1263\cplex\matlab\x64_win64';
else
    cplex_path = '';
    disp('Can not recognize the platform.');
end
addpath(genpath(cplex_path));

if showDebugInfo
    disp('GMCP LP Solver has started making the matrixes');
end

startT = tic;

%% Pre-processing the graph

% The net cost matrix has to be symmetrical.
netCostMat = (netCostMat + netCostMat')/2;

nVertices = sum(nClusterNodes);

% edgeLookup(i, j) contains the unique edge number between node i and
% node j.
% edgeLookup(i, j) == edgeLookupd(j, i)
edgeLookup = zeros(nVertices,nVertices);

% vertexLookup(i, j) contains the unique vertex number of node j in
% cluster i.
vertexLookup = zeros(nClusters, max(nClusterNodes));

nEdges = 0;
clusterLabel = [];
counterV = 0;
for i = 1:nClusters
    vertexLookup(i,1:nClusterNodes(i)) = counterV+1:counterV+nClusterNodes(i);
    counterV = counterV + nClusterNodes(i);
    clusterLabel = [clusterLabel;i*ones(nClusterNodes(i),1)];
    nEdges = nEdges + (nVertices-nClusterNodes(i))*nClusterNodes(i);
end
nEdges = nEdges/2;

if nVertices ~= size(netCostMat,1)
    warning('Size of netCostMat has to be the same as number of nodes.');
    keyboard;
    netCostMat = netCostMat(1:nVertices,1:nVertices);
end

%% Optimization Constraints
fE = zeros(nEdges + nVertices, 1);
x0 = zeros(nEdges+ nVertices, 1);
ctype = char([ones(1,nEdges)*('C'*1), ones(1,nVertices)*('B'*1)]);

AeqR2 = 0;
Aeq2 = zeros(nVertices, nEdges + nVertices);
beq2 = zeros(nVertices, 1);

edgeCounter = 0;
for i = 1 : nVertices
    iCluster = clusterLabel(i);
    AeqR2 = AeqR2 + 1;
    Aeq2(AeqR2, nEdges + i) = nClusters-1;
    for j = 1:nVertices
        jCluster = clusterLabel(j);
        if iCluster == jCluster
            continue;
        end
        if i<j
            edgeCounter = edgeCounter + 1;
            edgeLookup(i,j) = edgeCounter;
            edgeLookup(j,i) = edgeCounter;
            edgeList(edgeCounter,1) = i;
            edgeList(edgeCounter,2) = j;
            Aeq2(AeqR2, edgeLookup(i,j)) = -1;
            fE(edgeCounter) = netCostMat(i,j);
            if isnan(fE(edgeCounter))
                fE(edgeCounter) = INF_IN;
            end
        else
            Aeq2(AeqR2, edgeLookup(j,i)) = -1;
        end
    end
end

AeqR3 = 0;
Aeq3 = zeros(nClusters, nEdges + nVertices);

count = 0;
for C = 1:nClusters
    AeqR3 = AeqR3 + 1;
    Aeq3(AeqR3, nEdges+ count+1:nEdges+ count+nClusterNodes(C)) = 1;
    count = count + nClusterNodes(C);
end
beq3 = ones(AeqR3, 1) * kCliques;

AR1 = 0;
A1 = sparse(nClusters*(nClusters-1)*(nClusters-2)/2*max(nClusterNodes)^3, nEdges + nVertices);

AeqR4 = 0;
Aeq4 = zeros(sum(nClusterNodes)*(nClusters-1), nEdges + nVertices);
beq4 = zeros(size(Aeq4,1), 1);

for i = 1 : nClusters
    for j = 1 : nClusters
        if i ~= j
            for k = j+1 : nClusters
                if i ~= k
                    for ii = 1 : nClusterNodes(i) % Nodes in cluster i
                        for jj = 1 : nClusterNodes(j) % Nodes in cluster j
                            for kk = 1 : nClusterNodes(k) % Nodes in cluster k
                                AR1 = AR1 + 1;
                                A1(AR1, edgeLookup(vertexLookup(i,ii),vertexLookup(j,jj))) = 1;
                                A1(AR1, edgeLookup(vertexLookup(i,ii),vertexLookup(k,kk))) = 1;
                                A1(AR1, edgeLookup(vertexLookup(j,jj),vertexLookup(k,kk))) = -1;
                            end
                        end
                    end
                end
            end

            for ii = 1 : nClusterNodes(i)-1
                AeqR4 = AeqR4 + 1;
                Aeq4(AeqR4, edgeLookup(vertexLookup(i,ii), vertexLookup(j,1:nClusterNodes(j)))) = 1;
                Aeq4(AeqR4, nEdges + vertexLookup(i,ii)) = -1;
            end
        end
    end
end
A1 = A1(1:AR1,:);
b1 = ones(AR1, 1);

Aeq = [Aeq2;Aeq3;Aeq4];
beq = [beq2;beq3;beq4];

A = A1;
b = b1;

if showDebugInfo
    disp('Optimization has started');
end

try
    options = cplexoptimset();

    if isFindMax
        fE = - fE;
    end
    if showDebugInfo
        disp('everything is ready now, Optimization has started');
    end
    if method==1
        [x, fval, exitflag, output] = cplexbilp(fE, A, b, Aeq, beq, ...
            x0, options);
    elseif method ==0
        devideBy = 5;
        [x, fval, exitflag, output] = cplexlp(fE*50, A, b, Aeq, beq/devideBy, ...
            zeros(size(fE)),ones(size(fE))/devideBy,[ ], options);
    elseif method ==2
        [x, fval, exitflag, output] = cplexmilp(fE, A, b, Aeq, beq, [], [], ...
            [], zeros(size(fE)), ones(size(fE)), ctype, x0, options);
    elseif method == 3
        %It is already calculated before...
    end
catch
    error('If you do not have cplex on your machine, you might use the default matlab solver. This solver is much slower and has a limited number of variables.');
end

if any(abs(x - round(x)) > 10e-8)
    warning('There are some variables which differ at least 10e-8 from the integer numbers. Be careful about it, I will round the solution');
    keyboard;
end
x = round(x);

xOnlyNodes = x(nEdges + 1 : end);
nodesIdx = sort(find(xOnlyNodes == 1));

nodes = zeros(kCliques, nClusters);
for i = 1:kCliques
    nodes(i,1) = nodesIdx(i);
    neighborings = edgeLookup(nodesIdx(i), :);
    neighborings = nonzeros(neighborings);
    nodes(i,2:end) = edgeList(neighborings(x(neighborings) > 0.9)', 2);
end

if isFindMax
    cost = -fval;
else
    cost = fval;
end

timeSpent = toc(startT);

% _____________________________________
% Copyright (c) 2015, Shayan Modiri Assari
% 
% Permission to use, copy, modify, and distribute this software for research
% purposes with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
