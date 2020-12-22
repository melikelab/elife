function [KLDivergence] = KLEstimator1D(histogramP,histogramQ,minDist)
%KLEstimator1D Estimates the KL Divergence between two sets of IID random
%   variables from two distinct distributions
%   Inputs: histogramP: "True" set of random samples
%           histogramQ: "Null" set of random samples
%           minDist: minimum distance between samples.  Samples that are
%           closer than that are aggregated to a single data point with a higher weight
% Author: PKR, UPENN January 2020
% based off the Paper: Kullback-Leiber Divergence Estimation of Continuous
% Distributions by Fernando Perez-Cruz

%% TO DO, AGGREGATE POINTS THAT ARE TOO CLOSE TOGETHER.

% sort variables in ascending order as column vectors
histogramP = sort(histogramP(:));
histogramQ = sort(histogramQ(:));
% aggregate samples that are too close together (double counts, etc)
[adjHistP,weightsP] = weightedData(histogramP, minDist);
[adjHistQ,weightsQ] = weightedData(histogramQ, minDist);
% number of samples from the adjusted "true" distribution, basis of integrator
N = length(adjHistP);
% get boundary points for the coefficients
xMin = min(min(adjHistP),min(adjHistQ));
xMax = max(max(adjHistP),max(adjHistQ));
xSpacer = (xMax-xMin)/min(N,length(histogramQ));
x0 = xMin-xSpacer;
xF = xMax+xSpacer;
% get coefficients from linear interpolation of the empirical CDFs
[Ap,Bp] = linCDFCoeff(adjHistP,weightsP,x0,xF);
[Aq,Bq] = linCDFCoeff(adjHistQ,weightsQ,x0,xF);

diffP = diff(adjHistP);
delta = min(diffP)-eps;
% get euler derivative samples from histogram P
Px = Ap(1:N).*adjHistP+Bp(1:N);
Pm = Ap(1:N).*(adjHistP-delta)+Bp(1:N);
% get corresponding histQ coefficients for each of the sampled histP points
qxi = zeros(N,1);
qmi = zeros(N,1);
extendedHistQ = [adjHistQ;xF];
for ii = 1:N
   qxi(ii) = find(extendedHistQ>=adjHistP(ii),1);
   qmi(ii) = find(extendedHistQ>=(adjHistP(ii)-delta),1);
end
Qx = Aq(qxi).*adjHistP + Bq(qxi);
Qm = Aq(qmi).*(adjHistP-delta) + Bq(qmi);
% estimate and return the KL Divergence
DeltaP = Px-Pm;
DeltaQ = Qx-Qm;
KLDivergence = 1/N * sum(log(DeltaP)-log(DeltaQ)) - 1;
end

function [adjHist, weights] = weightedData(data, minDist)
weights = 0*data+1;
adjHist = data;
cc = 1;
for ii = 2:length(adjHist)
    if (data(ii)-adjHist(cc)) < minDist
        % combine the two data points, avg. position and shift the histogram
        adjHist(cc) = (adjHist(cc)+data(ii))/2;
        weights(cc) = weights(cc)+1; % stay at current node
    else
        cc=cc+1;
        adjHist(cc) = data(ii);
    end
end
% remove unused positions
adjHist(cc+1:end) = [];
weights(cc+1:end) = [];
end

function Nodes = empCDFNodes(weights)
N = sum(weights);
aN = length(weights);
Nodes = zeros(aN,1);
Nodes(1) = weights(1)-0.5;
for ii = 2:aN
    Nodes(ii) = Nodes(ii-1) + weights(ii);
end
Nodes = Nodes/N;
end

function [A,B] = linCDFCoeff(data,weights,x0,xF)
N = length(data);
Nodes = empCDFNodes(weights);
A = zeros(N+1,1);
B = zeros(N+1,1);
AltData = [x0;data;xF];
AltNodes = [eps;Nodes;1-eps];
for ii = 1:N+1
   NegDiff = 1/(AltData(ii)-AltData(ii+1));
   A(ii) = NegDiff*(AltNodes(ii)-AltNodes(ii+1));
   B(ii) = NegDiff*(-AltData(ii+1)*AltNodes(ii) ...
       + AltData(ii)*AltNodes(ii+1));
end
end

