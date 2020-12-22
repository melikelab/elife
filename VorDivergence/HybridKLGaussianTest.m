% HybridKLGaussianTest.m
% script to test the hybrid KL divergence estimator

% define parameters
mu = 0;
sigma = 1;
N = 20050;
%repeats = 100;
%errorplot = zeros(repeats,1);

% define function for analytical gaussian
GaussFunc = @(p,X) 1/sqrt(2*pi)/p(2) * exp(-(X-p(1)).^2/2/p(2)^2);
LogGFunc = @(p,X) -log(p(2))-0.5*log(2*pi) - (X-p(1)).^2/(2*p(2)^2);
%GaussCDFunc = @(p,X) 0.5*(1+erf((X-p(1))/(p(2)*sqrt(2))));
GaussCDist = @(X) LogGFunc([0,1],X);

Sample = mu + sigma*randn(N,1);

KLScore=KLEstimator1DHybrid(Sample,GaussCDist,1e-8)






