% VoronoiEntropyTest.m
% script to test the Voronoi Entropy of Shreyasi's data
% Author: PKR, UPenn, November 2018 


% load Shreya's data
File1 = 'O:\peter\storm1_0589-2086_ch1_iterVorSeg_VornClust_Th0,016.bin';
File2 = 'O:\peter\storm6_07046-06518_allChs_iterVorSeg_VornClust_Th0,016.bin';

LL1 = Insight3( File1 );
LL2 = Insight3( File2 );

xy1 = LL1.getXYcorr;
xy2 = LL2.getXYcorr;

Areas1 = zeros(length(xy1),1);
Areas2 = zeros(length(xy2),1);

% get area distributions for each set of points
% Data 1
DT1 = delaunayTriangulation(xy1);
[V1,C1] = voronoiDiagram(DT1);
VorDat1 = {V1, C1};

for pt=1:length(xy1)
    xt = V1(C1{pt},1);
    yt = V1(C1{pt},2);
    Areas1(pt) = abs(sum( (xt([2:end 1]) - xt).*(yt([2:end 1]) + yt))*0.5);
end
% Data 2
DT2 = delaunayTriangulation(xy2);
[V2,C2] = voronoiDiagram(DT2);
VorDat2 = {V2, C2};

for pt=1:length(xy2)
    xt = V2(C2{pt},1);
    yt = V2(C2{pt},2);
    Areas2(pt) = abs(sum( (xt([2:end 1]) - xt).*(yt([2:end 1]) + yt))*0.5);
end

% set a cut off to remove large areas for now
cutOff = 0.1;
modAreas1 = Areas1(Areas1<cutOff);
modAreas2 = Areas2(Areas2<cutOff);

% Get statistics on each distribution
meanArea1 = mean(modAreas1);
meanArea2 = mean(modAreas2);

% Figure out optimal bin sizes for each data set
[Mbest1,MAP1] = binEstimator(modAreas1);
[Mbest2,MAP2] = binEstimator(modAreas2);

h1 = histogram(modAreas1,Mbest1,'Normalization','Probability');
figure;
h2 = histogram(modAreas2,Mbest2,'Normalization','Probability');
%h1.BinWidth

% define likelihood function
Lfunc = @(a,b,c,X) a*b^(c/a) / gamma(c/a) * X.^(c-1) .* exp(-b*X.^a);

% determine coefficients
a = 1.07950;
b = 3.03226;
c = 3.31122;

% define negative log likelihood function for losses
nLLfunc = @(a,b,c,X) -(c/a)*log(b) - log(a) + gammaln(c/a) - (c-1)*log(X) + b*X.^a;

L1 = nLLfunc(a,b,c,modAreas1/meanArea1);
L2 = nLLfunc(a,b,c,modAreas2/meanArea2);

%EntropyVal = log(a*gamma(b/c) 

% xx = 1e-2:0.01:10;
% FuncVal1 = LLfunc(a,b,c,xx/meanArea1);
% FuncVal2 = LLfunc(a,b,c,xx/meanArea2);