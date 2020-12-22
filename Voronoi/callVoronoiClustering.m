% basic script to perform Voronoi Tesselation and cluster identification

clear
clc
% close all

% input parameters
% nanometers per pixel
pix2nm = 157;
% pixel size of the mask and for for cluster plotting 
p = 20;
% Cluster identification can be manual or automatic
% clusterID = 'automatic'; % performs Monte Carlo simulation
clusterID = 'manual'; 

% number of Monte Carlo iterations
iter = 5; % for automatic ID only
% percent confidence bounds for Monte Carlo output
signif = 99; % for automatic ID only
% minimum number of localizations allowed in each cluster
minNLoc = 5;
% true/false value for plotting the localizations, their Voronoi
% Tesselation and the identified clusters
showIM = true;

% starting place for selecting Localization list
startpath = 'V:\Jason_O\Software testing\ClusterViSu'; 



%%
%%%%%% begin algorithm %%%%%% 
[fname,fpath] = uigetfile(fullfile(startpath, '*.bin'));

% load localization list
LL = Insight3(fullfile(fpath,fname));
% get coordinates in units of pixels
xy = LL.getXYcorr;
% %%
switch clusterID
    case 'automatic'
        % send localizations for tesellation and perform monte carlo simulations
        % to identify the threshold value above which Voronoi polygon area is
        % indistinguishable from randomly distributed localizations. Polygons with
        % area below this threshold could be considered clusters
        [Histograms, intersection, xy, neighborList, maskBW, repidx] = VoronoiMonteCarlo_JO(xy,iter,signif);
        thresh = intersection(1);
        % plot the histograms of Voronoi Areas
        plotVoronoiMCdat(Histograms, thresh, signif);
        
        % send localizations for clustering
        if showIM
            [cluster, fig_h] = clusterVoronoi(xy,thresh,neighborList,minNLoc,[pix2nm,p],showIM,maskBW,LL.filename);
        else
            cluster = clusterVoronoi(xy,thresh,neighborList,minNLoc);
        end
        
    case 'manual'
        cluster = clusterVoronoi(xy,[],[],minNLoc,[pix2nm,p],showIM);
end
