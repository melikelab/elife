% basic script to perform Voronoi Tesselation and cluster identification

clear
clc
% close all

% input parameters
% nanometers per pixel
pix2nm = 157;
% pixel size of the mask and for for cluster plotting 
finalpix = 20;
pixVals = [pix2nm,finalpix];
% Cluster identification can be manual or automatic
clusterID = 'automatic'; % performs Monte Carlo simulation
% clusterID = 'manual'; 

% number of Monte Carlo iterations
iter = 5; % for automatic ID only
% percent confidence bounds for Monte Carlo output
signif = 99.5; % for automatic ID only
% minimum number of localizations allowed in each cluster
minNLoc = 5;
% true/false value for plotting the localizations, their Voronoi
% Tesselation and the identified clusters
showIM = true;

% starting place for selecting Localization list
startpath = 'V:\Jason_O\Software testing\ClusterViSu'; 



% %%
%%%%%% begin algorithm %%%%%% 
[fname,fpath] = uigetfile(fullfile(startpath, '*.bin'));

% %% load localization list
LL = Insight3(fullfile(fpath,fname));
% get coordinates in units of pixels
xy = LL.getXYcorr;
%%
switch clusterID
    case 'automatic'
        %%
        xyorig = xy;
timer1 = tic;
        % Use all localizations for Monte Carlo simulation of voronoi areas
        % from randomized localizations
%         [xy,DT,VorDat,repidx,neighborList] = VoronoiAreas(xy,true);
        [Histograms, intersection, xy, neighborList, mask, repidx, Varea_rnd] = ...
            VoronoiMonteCarlo_JO_beta(xyorig,iter,signif,pixVals);
t1 = toc(timer1)
        % save the output
        savefile = [LL.filename(1:end-4) '_iterVorSegData_beta.mat'];
        save(savefile,'Histograms', 'intersection', 'xy', 'neighborList',...
            'mask', 'repidx', 'Varea_rnd')
        %%
        unithresh = 2*((finalpix/pix2nm)^2)*mask.Area/size(xy,1);
        maxloop = 8;
        Allthresholds = iterativeVoronoiSegmentation(xy,Varea_rnd,maxloop,signif,unithresh);
        plotVoronoiMCdat(Histograms, Allthresholds(:,1), signif);
        set(gca,'XLim',[0 unithresh*1.6])
        %%

timer2 = tic;
        [Histograms, intersection, xy, neighborList, mask, repidx, Varea_rnd] = ...
            VoronoiMonteCarlo_JO(xyorig,iter,signif,pixVals);
t2 = toc(timer2)

% [t1 timer1;t2 timer2]
        % save the output
        savefile = [LL.filename(1:end-4) '_iterVorSegDataNew.mat'];
        save(savefile,'Histograms', 'intersection', 'xy', 'neighborList',...
            'mask', 'repidx', 'Varea_rnd','pix2nm','finalpix','pixVals','signif')
        %

        % Determine the threshold for area/localization assuming a uniform
        % distribution, a.la. SR-Tesseler
        unithresh = 2*((finalpix/pix2nm)^2)*mask.Area/size(xy,1);
%         [mask1.BW,mask1.Area] = Locs2Mask( xy(:,1:2), finalpix, pix2nm );
%         mask1.Area = ((finalpix/pix2nm)^2)*mask1.Area; % put back to units of pixels
%         unithresh=2*mask1.Area/size(xy,1);

        % Apply univorm threshold to establish the base object for sequential
        % cluster segmentation
        
        % send localizations with small Voronoi areas for seq seg
        maxloop = 8;
        % Apply univorm threshold to establish the base object for sequential
        % cluster segmentation
        % improperly, for program testing
%         MCthresholds = iterativeVoronoiSegmentation(xy( xy(:,3)<=unithresh, :),Varea_rnd,maxloop,signif,false);
%         Allthresholds = [unithresh;MCthresholds];
        % properly applying threshold
        Allthresholds = iterativeVoronoiSegmentation(xy,Varea_rnd,maxloop,signif,unithresh,true);
        % Don't apply any initial threshold
%         MCthresholds = iterativeVoronoiSegmentation(xy,Varea_rnd,maxloop,signif,true);
%         Allthresholds = [unithresh;MCthresholds];
        
        % save threshold info
        save(savefile,'Allthresholds','-append')

        % plot Voronoi distributions with obtained thresholds
        plotVoronoiMCdat(Histograms, Allthresholds(:,1), signif);
        set(gca,'XLim',[0 unithresh*1.6])
%%
%         % write a new molecule list assigningxy=xyorig; localizations to channels 0-9
%         % according to their Voronoi area        
        writeIterativeVoronoiSegmentedLL(LL,xy,Allthresholds)
        
        
        %% send localizations for clustering
        cluster = cell(size(Allthresholds,1),1);
        for th = 1:length(Allthresholds)
            if showIM
                [cluster{th}, fig_h] = ...
                    clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc,...
                    [pix2nm,finalpix],showIM,mask.BW,LL.filename);
            else
                cluster{th} = ...
                    clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc);
            end
            %% plot the cluster statistics
            figAnnotation = [LL.filename( find(LL.filename==filesep,1,'last')+1:end-4) ', thresh = ' num2str(Allthresholds(th,1),'%.3g') '    '];
            VClustFig = plotVoronoiClusterStats2(...
                cluster{th}.nLocs, ...
                cluster{th}.areas*(pix2nm^2), ...
                cluster{th}.NND(:,1)*pix2nm, ...
                figAnnotation);
            %% save the clusters in a localization list for visualization
            writeVoronoiClusteredLL(LL,cluster{th},Allthresholds(th,1))
        end
        %%
        if size(Allthresholds,1) > 1
            plotVoronoiArea_nLocs(cluster,Allthresholds,LL.filename( find(LL.filename==filesep,1,'last')+1:end-4))
        end
        
    case 'manual'
        cluster = clusterVoronoi(xy,[],[],minNLoc,[pix2nm,finalpix],showIM);
end
