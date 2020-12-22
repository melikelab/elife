function [KLDivergence] = KLEstimator1DHybrid(histogramP,logPDFuncQ,minDist)
%KLEstimator1D Estimates the KL Divergence between a set of IID random
%   variables from a distinct distribution and a theoretical PDF function
%   of another distribution (the "Null Hypothesis")
%   Inputs: histogramP: "True" set of random samples
%           logPDfunc: function of the log PDF of the null hypothesis
%           minDist: minimum distance between samples.  Samples that are
%           closer than that are aggregated to a single data point with a higher weight
% Author: PKR, UPENN January 2020
% based off the Paper: Kullback-Leiber Divergence Estimation of Continuous
% Distributions by Fernando Perez-Cruz, modified for a theoretical Q

% sort variables in ascending order as column vectors
histogramP = sort(histogramP(:));
weightsP = 0*histogramP+1;
adjHistP = histogramP;
cc = 1;
for ii = 2:length(adjHistP)
    if (histogramP(ii)-adjHistP(cc)) < minDist
        % combine the two data points, avg. position and shift the histogram
        adjHistP(cc) = (adjHistP(cc)+histogramP(ii))/2;
        weightsP(cc) = weightsP(cc)+1; % stay at current node
    else
        cc=cc+1;
        adjHistP(cc) = histogramP(ii);
    end
end
% remove unused positions
adjHistP(cc+1:end) = [];
weightsP(cc+1:end) = [];
% number of samples from the "true" distribution, basis of integrator
N = length(adjHistP);
% get boundary points for the coefficients
xMin = min(adjHistP);
xMax = max(adjHistP);
xSpacer = (xMax-xMin)/N;
x0 = xMin-xSpacer;
xF = xMax+xSpacer;
% get coefficients from linear interpolation of the empirical CDFs
[Ap,Bp] = linCDFCoeff(adjHistP,weightsP,x0,xF);

diffP = diff(adjHistP);
delta = min(diffP)-eps;
% get euler derivative samples from histogram P
Px = Ap(1:N).*adjHistP+Bp(1:N);
Pm = Ap(1:N).*(adjHistP-delta)+Bp(1:N);
% Evaluate Delta at corresponding P points
logDeltaQ = logPDFuncQ(adjHistP-delta/2)+log(delta);
% estimate and return the KL Divergence
DeltaP = Px-Pm;
KLDivergence = 1/N * sum(log(DeltaP)-logDeltaQ) - 0.5772;
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