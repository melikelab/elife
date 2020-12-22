% clusterBatchFilterA.m
% script to filter a batch of voronoi.mat files and output corresponding
% .bin files with localizations removed

% user threshold params, threshold by area in px^2
maxArea = 500;
minArea = 20;

% load the .mat files
d = uigetdir(pwd, 'Select a folder');
files = dir(fullfile(d, '*.mat'));
aa=size(files,1);

%Insight3 data formatting
chidx = 12;
xyidx = [1,2,3,4];
columnValues = [5,100;6,10000;7,300;9,1;11,10000;13,1;14,1;15,1;16,-1];

% lets load the first file for now
for f = 1:aa
    fileLoc = fullfile(files(f).folder,files(f).name);
    clust = load(fileLoc);
    cluster = clust.cluster;
    fields = fieldnames(cluster{1});
    numThresh = length(cluster);
    xy = clust.xy;
    dat = zeros(size(xy,1),18);
    dat(:,xyidx) = [xy(:,[1,2]),xy(:,[1,2])];
    for ii = 1:size(columnValues,1)
        dat(:,columnValues(ii,1)) = columnValues(ii,2);
    end
    
    for th = 1:numThresh
        thresh = clust.Allthresholds(th);
        thresh = sprintf('_Th%.2g',thresh);
        thresh( thresh=='.' ) = ',';
        tempCluster = cluster{th};
        numClust = length(tempCluster.nLocs);
        filter = tempCluster.areas > maxArea | tempCluster.areas < minArea;
        for fn = 1:length(fields)
            tempCluster.(fields{fn})(filter,:) = [];
        end
        % save the clusters in a localization list for visualization
        %writeVoronoiClusteredLL(LL,tempCluster,clust.Allthresholds(th), repidx)
        dat(:,chidx) = nan; % clear out the channel information
        % define pseudo localization list for cluster centers
        % channel zero default
        nclust = size(tempCluster.center,1);
        centers = zeros(nclust,size(dat,2));
        for ii = 1:size(columnValues,1)
            centers(:,columnValues(ii,1)) = columnValues(ii,2);
        end
        centers(:,xyidx) = [tempCluster.center, tempCluster.center];
        
        % initialize new data array
        LocsInClusters = cell2mat(tempCluster.Locs);
        datnew = nan( size(LocsInClusters,1),size(dat,2) );
        j = 1;
        ch = 1;
        % %%
        for c = 1:nclust
            % extract localization indicies
            locIdx = tempCluster.Locs{c,1}(:,5);
            nlocs = size(locIdx,1);
            % set the centroid information for cluster 'c'
            datnew(j,:) = centers(c,:);
            % set the localizations information for cluster 'c'
            datnew(j+1:j+nlocs,:) = dat(locIdx,:);
            datnew(j+1:j+nlocs,chidx) = ch;
            j = j+1+nlocs;
            if ch == 9
                ch = 1;
            else
                ch = ch+1;
            end
            
        end
        % ensure all the entered data has a channel info
        if ~isempty( find(isnan(datnew(:,chidx)), 1) )
            % the channel could only be nan if it came from the dat array without
            % being assigned via the cluster.Locs{c,1}(:,5) localization index
            error('Error writing cluster data to localization list: undefined channel info')
        end
        % write a new localization list
        LLnew = Insight3();
        LLnew.setData( datnew );
        LLnew.setFilename( [fileLoc '_VFiltered.bin'] );
        LLnew.forceFileOverwrite(true);
        LLnew.write();
        
    end
end