clear
% clc
saveOutFile = 'VoronoiAreasPX2.csv';
[pixVals,pix2nm,finalpix,clusterID,iter,signif,minNLoc,maxloop,showIM,files] = ...
    setVoronoiClusterParams();
%     %%
nfiles = size(files.data,1);
% If doing non auto thresholding, ask user for threshold(s)
switch clusterID
    case 'manual'
        options.WindowStyle = 'normal';
        manualThresh = inputdlg([...
            'Set the threshold(s) [pix^2] you desire for cluster identification.  '...
            'Can be an array of thresholds.'],...
            'Manual Voronoi Threshold',[1 68],{''},options);
        clear options
        if ~isempty( manualThresh )
            if strcmp(manualThresh{1},''), manualThresh = [];
            else
                manualThresh = str2num( manualThresh{1} );
            end
        end
    case 'cdf'
        options.WindowStule = 'normal';
        cdfThresh = inputdlg([...
            'Set the cdf threshold(s) [% of kept areas] for cluster identification. ', ...
            'Can be an array of thresholds.'],...
            'CDF Voronoi Threshold',[1 68],{''},options);
        clear options
        if ~isempty(cdfThresh)
            if strcmp(cdfThresh{1},''), cdfThresh = [];
            else
                cdfThresh = str2num( cdfThresh{1} );
            end
        end
    case 'pval'
        options.WindowStule = 'normal';
        pvThresh = inputdlg([...
            'Set the pvalues [% of kept uniform areas] for cluster identification. ', ...
            'Can be an array of thresholds.'],...
            'PValue Voronoi Threshold',[1 68],{''},options);
        clear options
        if ~isempty(pvThresh)
            if strcmp(pvThresh{1},''), pvThresh = [];
            else
                pvThresh = str2num( pvThresh{1} );
            end
        end
end
% create an outervariable to store third column of xy variable
xythirddata = cell(nfiles,2);
for ii = 1:nfiles
    xythirddata{ii,1} = files.data{ii,1};
end

%%%%%% begin %%%%%%
for f = 1:nfiles
    % load localization list
    LL = Insight3( fullfile( files.data{f,2},files.data{f,1} ) );
    % get coordinates in units of pixels
    xy = LL.getXYcorr;
    if size(xy,1)<3
        continue
    end
    % make directory for saving data
    [fpath,fname,~] = fileparts( LL.filename );
    spath = [fpath filesep 'test_new'];
    if ~exist(spath,'dir')
        mkdir(spath)
    end
    badIndex = [];
%     if size(xy,1)>= 3  
    switch clusterID
        case 'automatic'           
            timer1 = tic;
            [Histograms, ~, xy, neighborList, mask, DT, VorDat, repidx, Varea_rnd] = ...
                VoronoiMonteCarlo_JO(xy,iter,signif,pixVals);
            t1 = toc(timer1);
            % save the output
            savefile = fullfile(spath, [fname '_iterVorSegData.mat']);
            save(savefile,...
                'pix2nm','finalpix','pixVals','clusterID','iter','signif','minNLoc','maxloop','showIM',... inputs
                'Histograms', 'xy', 'neighborList', 'mask', 'DT', 'VorDat', 'repidx', 'Varea_rnd',... outputs
                'savefile')
            %% Determine the threshold for area/localization assuming a uniform
            % distribution, a.la. SR-Tesseler
            unithresh = 2*((finalpix/pix2nm)^2)*mask.Area/size(xy,1);
            % send localizations with small Voronoi areas for seq seg
            % Apply uniform threshold to establish the base object for sequential
            % cluster segmentation
            Allthresholds = iterativeVoronoiSegmentation(xy,Varea_rnd,maxloop,signif,unithresh,showIM);
            % save threshold info
            save(savefile,'Allthresholds','unithresh','-append')
            % plot Voronoi distributions with obtained thresholds
            plotVoronoiMCdat(Histograms, Allthresholds(:,1), signif);
            set(gca,'XLim',[0 unithresh*1.6])
            % save cluster area vs Loc/cluster image as .png file
            saveas(gcf,[savefile(1:end-8) '_thresholds.png'],'png')
            
        
      
        case 'manual'
            
            
            [xy,~,~,repidx,neighborList] = VoronoiAreas(xy,true);
            [mask.BW,~] = Locs2Mask( xy(:,1:2));
            mask.BW = imresize(mask.BW,pix2nm/finalpix);
            % write a new molecule list assigning localizations to channels 0-9
            % according to their Voronoi area
            manualThresh = sort(manualThresh,'descend');
            Allthresholds = manualThresh(:);
            savefile = fullfile(spath, [fname '_manualVorSegData.mat']);
            save(savefile,...
                'clusterID','manualThresh','minNLoc','pix2nm','finalpix','showIM',... inputs
                'xy', 'Allthresholds', 'repidx', ... outputs
                'savefile')
            
            
        case 'cdf'
            [xy,~,~,repidx,neighborList] = VoronoiAreas(xy,true);
            [mask.BW,~] = Locs2Mask( xy(:,1:2));
            mask.BW = imresize(mask.BW,pix2nm/finalpix);
            Vareas = xy(:,3);
            [pvals, areaorder] = ecdf(Vareas);
            AreaThresh = zeros(length(cdfThresh),1);
            for ii = 1:length(cdfThresh)
                pvalInd = pvals <= cdfThresh(ii)/100;
                AreaThresh(ii) = max(areaorder(pvalInd));
            end
            % write a new molecule list assigning localizations to channels 0-9
            % according to their Voronoi area
            AreaThresh = sort(AreaThresh,'descend');
            Allthresholds = AreaThresh(:);
            savefile = fullfile(spath, [fname '_cdfVorSegData.mat']);
            save(savefile,...
                'clusterID','cdfThresh','AreaThresh','minNLoc','pix2nm','finalpix','showIM',... inputs
                'xy', 'Allthresholds', 'repidx', ... outputs
                'savefile')            
        case 'pval'
            [xy,~,~,repidx,neighborList] = VoronoiAreas(xy,true);
            [mask.BW,~] = Locs2Mask( xy(:,1:2));
            mask.BW = imresize(mask.BW,pix2nm/finalpix);
            Vareas = xy(:,3);
            % get the boundary on the points
            K = boundary(xy(:,1),xy(:,2));
            Valids = true(length(xy(:,1)),1);
            Valids(K) = 0;
            VVareas = Vareas(Valids);
            Scale = geomean(Vareas(Valids)); % scale to adjust threshold
            % hyper parameters for the distribution of random data
            a = 1.07950;
            b = 3.03226;
            c = 3.31122;
            % Threshold by pvalue
            AreaThresh = zeros(length(pvThresh),1);
            for ii = 1:length(pvThresh)
                AreaThresh(ii) = Scale*(gammaincinv(pvThresh(ii)/100,c/a)/b)^(-1/a);            
            end
            % write a new molecule list assigning localizations to channels 0-9
            % according to their Voronoi area
            AreaThresh = sort(AreaThresh,'descend');
            Allthresholds = AreaThresh(:);
            savefile = fullfile(spath, [fname '_pvalVorSegData.mat']);
            save(savefile,...
                'clusterID','pvThresh','AreaThresh','minNLoc','pix2nm','finalpix','showIM',... inputs
                'xy', 'Allthresholds', 'repidx', ... outputs
                'savefile')           
    end % main cluster switch
    writeIterativeVoronoiSegmentedLL(LL,xy,repidx,Allthresholds,spath)
    % send localizations for clustering
    cluster = cell(size(Allthresholds,1),1);
    VClustFig = gobjects(size(Allthresholds,1),1);
    for th = 1:size(Allthresholds,1)
        if showIM
            [cluster{th}, fig_h] = ...
                clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc,...
                [pix2nm,finalpix],showIM,mask.BW,LL.filename);
            % save cluster stats for application of threshold 'th'
            switch clusterID
                case 'automatic'
                    saveas(gcf,[savefile(1:end-8) '_VoronoiThresh#' num2str(th) '.png'],'png')
                case 'manual'
                    saveas(gcf,[savefile(1:end-8) '_manVoronoiThresh#' num2str(th) '.png'],'png')
                case 'cdf'
                    saveas(gcf,[savefile(1:end-8) '_cdfVoronoiThresh#' num2str(th) '.png'],'png')
                case 'pval'
                    saveas(gcf,[savefile(1:end-8) '_pvalVoronoiThresh#' num2str(th) '.png'],'png')
            end
        else
            try
                cluster{th} = ...
                    clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc);
            catch ME
                msg=('All coordinates are on the same line');
                causeException = MException('DELAUNAY:triangulationEmpty',msg);
                ME = addCause(ME,causeException);
                cluster{th} = [];
            end
        end
        % plot the cluster statistics
        switch clusterID
            case 'automatic'
                figAnnotation = [files.data{f,1}(1:end-4) ', thresh = ' num2str(Allthresholds(th,1),'%.3g') '    '];
            case 'manual'
                figAnnotation = [files.data{f,1}(1:end-4) ', thresh = ' num2str(Allthresholds(th),'%.3g') '    '];
            case 'cdf'
                figAnnotation = [files.data{f,1}(1:end-4) ', CDF Thresh = ' num2str(Allthresholds(th),'%.3g') '    '];
            case 'pval'
                figAnnotation = [files.data{f,1}(1:end-4) ', pvalue Thresh = ' num2str(Allthresholds(th),'%.3g') '    '];
        end
        if ~isempty(cluster{th}) && ~isscalar(cluster{th}.nLocs) && ~isscalar(cluster{th}.NND(:,1))
            VClustFig(th) = plotVoronoiClusterStats2(...
                cluster{th}.nLocs, ...
                cluster{th}.areas*(pix2nm^2), ...
                cluster{th}.NND(:,1)*pix2nm, ...
                figAnnotation);
            % save cluster stats for application of threshold 'th'
            switch clusterID
                case 'automatic'
                    saveas(VClustFig(th),[savefile(1:end-8) '_ClustThresh#' num2str(th) '.png'],'png')
                case 'manual'
                    saveas(VClustFig(th),[savefile(1:end-8) '_ManualVorThresh#' num2str(th) '.png'],'png')
                case 'cdf'
                    saveas(VClustFig(th),[savefile(1:end-8) '_cdfVorThresh#' num2str(th) '.png'],'png')
                case 'pval'
                    saveas(VClustFig(th),[savefile(1:end-8) '_pvalVorThresh#' num2str(th) '.png'],'png')
            end
            % save the clusters in a localization list for visualization
            writeVoronoiClusteredLL(LL,cluster{th},Allthresholds(th,1), repidx)
        else
            badIndex = [badIndex, th];
        end
    end
    % save cluster info
    save(savefile,'cluster','-append')
    
    if size(Allthresholds,1) > 1
        goodCluster = cluster;
        goodCluster(badIndex) = [];
        goodThresholds = Allthresholds;
        goodThresholds(badIndex) = [];
        if isempty(goodCluster)
            warning('No valid clusters to plot by Area');
            continue
        end
        plotVoronoiArea_nLocs(cluster,Allthresholds,pix2nm,LL.filename( find(LL.filename==filesep,1,'last')+1:end-4))
        % save cluster area vs Loc/cluster image as .png file
        saveas(gcf,[savefile(1:end-8) '_AreaVLoc.png'],'png')
    end
    % add to the xythirddata cell array
    xythirddata{f,2} = xy(:,3);
    close all;  
   % end
end
% file loop
if size(xy,1) >=3
header = xythirddata(:,1);
dataList = xythirddata(:,2);
[max_size, max_index] = max(cellfun('size', dataList, 1));
upLoadList = cell(length(dataList));
for ii = 1 : length(dataList)
    upLoadList{ii} = -1*ones(max_size,1);
    upLoadList{ii}(1:length(dataList{ii})) = dataList{ii};
end
% write output of xy to a file
commaHeader = [header repmat({','},numel(header),1)]'; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
%write header to file
saveOut = fullfile(spath,saveOutFile);
fid = fopen(saveOut,'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
%write data to end of file
dlmwrite(saveOut,upLoadList','-append');
end
