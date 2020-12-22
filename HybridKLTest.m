% HybridKLTest.m
% script to perform the KL divergence test using an analytical Q and an
% empirical P, where Q is described by the generalized gamma distribution
% Author: PKR, UPENN, September 2020

% initial functions and parameters to define
% parameters from Tanemura 2003 paper
a_uni = 1.07950;
b_uni = 3.03226;
c_uni = 3.31122;
uniCoeff = [a_uni, b_uni, c_uni];

% log PDF of generalized gamma distribution in log area space
PDFgamlog = @(p,y) p(3)/p(1) * log(p(2)) + log(p(1)) - gammaln(p(3)/p(1)) ...
    + p(3)*y - p(2)*exp(y*p(1));
UNIgamlog = @(y) PDFgamlog(uniCoeff,y);

%% data selection of PostVA files
[fileNames,path] = uigetfile('*.mat',  'All Files (*.*)','MultiSelect','on');

numFiles = length(fileNames);
estCoeff = cell(numFiles,1);
stdCoeff = cell(numFiles,1);

normedAreas = cell(numFiles,1);
logAreas = cell(numFiles,1);
Divergence = zeros(numFiles,1);

Colors = [0.49 0.18 0.56; 0.93 0.69 0.13; 0.85 0.32 0.09 ;0 0.48 0.74]; %(Purple, Yellow, Brown)

% loop over each file name
for ii = 1:numFiles
   % load file name 
    tempFile = load(fullfile(path,fileNames{ii}));
    
    normedAreas{ii} = [];
    logAreas{ii} = [];
    for kk = 1:length(tempFile.VoronoiAreas)
        normedAreas{ii} = [normedAreas{ii}; ...
            tempFile.VoronoiAreas{kk}/mean(tempFile.VoronoiAreas{kk})];
    end
    logAreas{ii} = log(normedAreas{ii});
    
    % perform KL divergence of the histograms
    Divergence(ii) = KLEstimator1DHybrid(logAreas{ii},UNIgamlog,1e-10)
end
% for ii = 1:numFiles
% histogram(log(normedAreas{ii}),100,'normalization', 'probability');
% hold on
% %histogram(log(normedAreas{2}),100,'normalization', 'probability');
% end
% hold off

% histogram(log(normedAreas{1}),100,'normalization', 'pdf','FaceColor','[0.7 0.6940 0.1250]','EdgeColor','[0 0 0]');
% hold on
% histogram(log(normedAreas{2}),100,'normalization', 'pdf','FaceColor','[0.2 0.69 0.5]','EdgeColor','[0 0 0]');
% hold off

histogram(log(normedAreas{1}),100,'normalization', 'pdf','FaceColor','[.3 .3 .3]','EdgeColor','[.3 .3 .3]');
hold on
%histogram(log(normedAreas{2}),100,'normalization', 'pdf','FaceColor','[.7 .7 .7]','EdgeColor','[.7 .7 .7]');
hold off

% histogram(log(normedAreas{1}),100,'normalization', 'pdf','FaceColor','[0.2 0.69 0.5]','EdgeColor','[[0.2 0.69 0.5]');
% hold on
% histogram(log(normedAreas{2}),100,'normalization', 'pdf','FaceColor','[.1 .35 .25]','EdgeColor','[.1 .35 .25]');
% hold off

%Divergence12 = KLEstimator1D(logAreas{1},logAreas{2},1e-10)

