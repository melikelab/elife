% NaiveKLTest.m
% script to perform the KL divergence test using a naive comparison to the
% generalized gamma distribution
% Author: PKR, UPENN, November 2019

% initial functions and parameters to define
% parameters from Tanemura 2003 paper
a_uni = 1.07950;
b_uni = 3.03226;
c_uni = 3.31122;
uniCoeff = [a_uni, b_uni, c_uni];
% function for generalized gamma distribution
LFunc = @(p,X) p(1)*p(2)^(p(3)/p(1)) / gamma(p(3)/p(1)) * X.^(p(3)-1) .* exp(-p(2)*X.^p(1));
% function for the negative log likelihood of the generalized gamma dist
nLLfunc = @(p,X) -(p(3)/p(1))*log(p(2)) - log(p(1)) + gammaln(p(3)/p(1)) ...
    - (p(3)-1)*log(X) + p(2)*X.^p(1);
% entropy function of gamma dist
EntropyVal = @(p) log(p(2)^(-1/p(1))*gamma(p(3)/p(1))/p(1)) + p(3)/p(1) + (1/p(1) - p(3)/p(1))*psi(p(3)/p(1));
% KL-Divergence between two generalize gamma distributions
KLFunc = @(p1,p2) log(p1(1)*p2(2)^(-p2(3)/p2(1)) * gamma(p2(3)/p2(1)) / ...
    (p2(1)*p1(2)^(-p1(3)/p1(1))*gamma(p1(3)/p1(1)))) + ...
    (psi(p1(3)/p1(1))/p1(1) - 1/p1(1)*log(p1(2)) )*(p1(3)-p2(3)) ...
    + gamma( (p1(3)+p2(1))/p1(1)) / gamma( p1(3)/p1(1)) * (p1(2)^(-1/p1(1))/p2(2)^(-1/p2(1)))^(p2(1)) ...
    - p1(3)/p1(1);

% try the Weibull++ parameterization with log likelihood
nWeib = @(p,X) log(X) + log(p(2)) - log(abs(p(3))) + gammaln(1/(p(3)^2)) ...
    - ( ( (p(3)*log(X)-p(1))/p(2) + log(1/(p(3)^2)) - exp(p(3)*(log(X)-p(1))/p(2)) )/p(3)^2);

uniWeib = [-log(b_uni)/a_uni + 1/a_uni * log(c_uni/a_uni), ...
    1/sqrt(c_uni*a_uni), sqrt(a_uni/c_uni)];

%% data selection of PostVA files
[fileNames,path] = uigetfile('*.mat',  'All Files (*.*)','MultiSelect','on');

numFiles = length(fileNames);
estCoeff = cell(numFiles,1);
stdCoeff = cell(numFiles,1);

normedAreas = cell(numFiles,1);

% loop over each file name
for ii = 1:numFiles
   % load file name 
    tempFile = load(fullfile(path,fileNames{ii}));
    
    normedAreas{ii} = [];
    for kk = 1:length(tempFile.VoronoiAreas)
        normedAreas{ii} = [normedAreas{ii}; ...
            tempFile.VoronoiAreas{kk}/mean(tempFile.VoronoiAreas{kk})];
    end
    
    % perform MLE estimation of Areas Distribution
    %nLLsum = @(p) sum(nLLfunc(p,normedAreas));
    nWeibSum = @(p) sum(nWeib(p,normedAreas{ii}));
    estCoeff{ii} = fminsearch(nWeibSum,uniWeib);
    
    % export coefficients to shorter variables for typing out conversion
    tmu = estCoeff{ii}(1);
    tsi = estCoeff{ii}(2);
    tla = estCoeff{ii}(3);
    
    stdCoeff{ii} = [tla/tsi, 1/tla^2*exp(-tla*tmu/tsi), 1/tsi/tla];
    
    % plot a histogram of fit to data
    figure;
    higram(ii) = histogram(normedAreas{ii},200,'Normalization','pdf','EdgeAlpha',0);
    hold on
    funcVal{ii} = LFunc(stdCoeff{ii},higram(ii).BinEdges);
    plot(higram(ii).BinEdges,funcVal{ii},'r','LineWidth',3);
    title('Histogrammed Areas of Data and MLE fit','FontSize',14)
    xlabel('Reduced Areas','FontSize',14)
    ylabel('Probability Density','FontSize',14)
    hold off
end

% pull out KL divergences
Divergence = zeros(numFiles,1);
for ii = 1:numFiles
    Divergence(ii) = KLFunc(stdCoeff{ii},uniCoeff);
end

