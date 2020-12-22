% KLEstimatorTest.m
% script to test the KLEstimator1D function on procedurally generated
% histograms given a known divergence score

% Author PKR Upenn January 2020
mu = 0;
sigma = 1;
N = 20050;
M = 20000;
repeats = 100;
errorplot = zeros(repeats,1);
for ii = 1:repeats
   histP = randn(N,1)*sigma+mu;
   histQ = randn(M,1)*sigma+mu;
   histP = sort(histP);
   histQ = sort(histQ);
   KLest = KLEstimator1D(histP,histQ,1e-8);
   errorplot(ii) = 0-KLest;
end

plot(errorplot)