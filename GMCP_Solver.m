function [nodes, cost, timeSpent] = GMCP_Solver(Net_Cost_Mat, NN, NC, ...
    isFindMax, method, showDebugInfo, Kcliques)

% Variable List:
% Net_Cost_Mat is the weights of the GMCP graph! Please keep the nodes 
%   belonging to each cluster together. 
%
% NN: Number of Nodes in each Cluster. It can be a single number
%   (For Graphs with the same number of NNs in each Cluster) or can be a
%   vector of numbers which shows number of members for each Cluster
%
% NC: Number of Clusters
%
% isFindMax: 1 if finding Maximum, 0 if finding Minimum
%
% method: 0 if Linear Programmin
%         1 if Binary Linear Programming
%     *** 2 (Default-Fast-Exact) if mixed Linear Programming
%
% Kcliques: Number of cliques in the solution. In case of GMCP, this
%   Kcliques = 1 and for GMMCP Kcliques should be selected appropriately.
%
% Pleaes cite:
% 1- Shayan Modiri Assari, Amir Roshan Zamir, and Mubarak Shah. "Video classification using semantic concept co-occurrences." Computer Vision and Pattern Recognition (CVPR), 2014 IEEE Conference on
% 2- Afshin Dehghan, Shayan Modiri Assari, and Mubarak Shah. "GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking." CVPR. Vol. 1. 2015.

infin = 1e6;

Net_Cost_Mat = (Net_Cost_Mat + Net_Cost_Mat')/2;

if ~exist('method', 'var')
    method = 2;
end

if ~exist('showDebugInfo','var')
    showDebugInfo = 0;
end

if ~exist('Kcliques', 'var') || isempty(Kcliques)
    Kcliques = 1;
end

if length(NN) == 1
    NN = NN*ones(1,NC);
end

if length(NN) ~= NC
    warning('Length of NN needs to be consistant with the provided NC');
    keyboard;
end

if showDebugInfo
    disp('GMCP LP Solver has started making the matrixes');
end

startT = tic;

VN = sum(NN);
NNcumsum = cumsum(NN);

Enum = zeros(VN,VN);
Vnum = zeros(NC, max(NN));

EN = 0;
ClusterLabel = [];
counterV = 0;
for i = 1:NC
    Vnum(i,1:NN(i)) = counterV+1:counterV+NN(i);
    counterV = counterV + NN(i);
    ClusterLabel = [ClusterLabel;i*ones(NN(i),1)];
    EN = EN + (VN-NC-(NN(i)-1))*(NN(i)-1);
end
EN = EN/2;

if VN ~= size(Net_Cost_Mat,1)
    warning('Size of Net_Cost_Mat has to be the same as number of nodes.');
    keyboard;
    Net_Cost_Mat = Net_Cost_Mat(1:VN,1:VN);
end


E = 0;
fE = zeros(EN + NC, 1);
x0 = zeros(EN+ NC, 1);
ctype = char([ones(1,EN)*('B'*1), ones(1,NC)*('C'*1)]);

AeqR1 = 0;
Aeq1 = zeros(NC, EN + NC);
beq1 = ones(NC, 1)*Kcliques*(NC-1);

for i = 1:VN
    Ci = ClusterLabel(i);
    if (i==NNcumsum(Ci))
        Aeq1(Ci, EN+Ci) = 1;
        continue;
    end
    for j = 1:VN
        
        Cj = ClusterLabel(j);
        if (j==NNcumsum(Cj))
            continue;
        end
        
        if Ci == Cj
            continue;
        end
        if i<j
            E = E + 1;
            Enum(i,j) = E;
            Enum(j,i) = E;
            Elist(E,1) = i;
            Elist(E,2) = j;
            fE(E) = Net_Cost_Mat(i,j);
            if isnan(fE(E))
                fE(E) = infin;
            end
        else
        end
        Aeq1(Ci,Enum(i,j)) = 1;
    end
end


AR1 = 0;
A1 = sparse(NC*(NC-1)*(NC-2)/2*max(NN)^3, EN + NC);

AR2 = 0;
A2 = zeros(VN*(NC-1), EN+NC);
b2 = ones(VN*(NC-1),1);


for i = 1:NC
    for j = 1:NC
        if i ~= j
            for k = j+1:NC
                if i ~= k
                    for ii = 1:NN(i)-1
                        for jj = 1:NN(j)-1
                            for kk = 1:NN(k)-1
                                AR1 = AR1 + 1;
                                A1(AR1, Enum(Vnum(i,ii),Vnum(j,jj))) = 1;
                                A1(AR1, Enum(Vnum(i,ii),Vnum(k,kk))) = 1;
                                A1(AR1, Enum(Vnum(j,jj),Vnum(k,kk))) = -1;
                            end
                        end
                    end
                end
            end
            
            for ii = 1:NN(i)-1
                AR2 = AR2 + 1;
                A2(AR2, Enum(Vnum(i,ii), Vnum(j,1:NN(j)-1))) = 1;
            end
        end
    end
end
A1 = A1(1:AR1,:);
b1 = ones(AR1, 1);

Aeq = [Aeq1];
beq = [beq1];

A = [A1;sparse(A2)];
b = [b1;b2];

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
            [], zeros(size(fE)), [], ctype, x0, options);
    elseif method == 3
        %It is already calculated before...
    end
catch
    warning('If you do not have cplex on your machine, you might use the default matlab solver. This solver is much slower and has a limited number of variables.');
    keyboard;
end

if any(abs(x-round(x))>10e-8)
    warning('There are some variables which differ at least 10e-12 from the integer numbers. Be careful about it, I will round the solution');
    keyboard;
end
x = round(x);

x2 = x(1:EN);
x2 = sort(find(x2==1));
x2 = Elist(x2,:);
x2 = unique(x2(:));

counter = 0;
nodes = zeros(Kcliques, NC);
while length(x2)
    counter = counter + 1;
    nodes(counter,ClusterLabel(x2(1))) = x2(1);
    tmp = Enum(x2(1),:);
    tmp = nonzeros(tmp);
    tmp = Elist(tmp(x(tmp)>0.9)',2);
    nodes(counter,ClusterLabel(tmp)) = tmp;
    x2 = x2(~ismember(x2,nodes(counter,:)));
end

Net_Cost_Mat = Net_Cost_Mat + dummyWeight;
Net_Cost_Mat(isnan(Net_Cost_Mat)) = 0;

x2 = find(~ismember(1:VN, NNcumsum));
x2 = x2(~ismember(x2,nonzeros(nodes(:))));
[~,sInd] = sort(sum(Net_Cost_Mat(x2,:),2));
x2 = x2(sInd);

for i = counter + 1 : min(Kcliques, counter + length(x2))
    nodes(i,ClusterLabel(x2(i-counter))) = x2(i-counter);
end

for i = 1:NC
    zeroInd = (nodes(:,i) == 0);
    nodes(zeroInd, i) = NNcumsum(i);
end

cost = 0;
for i = 1:size(nodes,1)
    tmp = Net_Cost_Mat(nodes(i,:),nodes(i,:));
    tmp = tmp .* (1-eye(size(tmp)));
    cost = cost + sum(sum(tmp));
end
cost = cost / 2;

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
