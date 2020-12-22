
% Inputs
%   xy = either a 2-column list of [x,y] positions or the output from 
%       VoronoiAreas.m
%   thresh = upper Voronoi polygon area threshold for cluster determination
%   neighborList = output of VoronoiAreas.m that indicates which
%       localizations neighbor one another
%   minNLoc = minimum number of lcoalizations allowed in each cluster
%   convFact = 1- or 2-element array with the first value being pixels to 
%       nm conversion, for converting [x,y] positions from pixels to
%       nanometers; the second element (if desired) being the pixel size of 
%       the final mask.
%   showIM = binary value true/false to show Voronoi diagram with
%       localizations whose area > thresh and identified clusters
%   mask = binary mask generated from input [x,y] positions in columns 1 &
%       2 of xy
%   figTitle = string for the title of the figure, if showIM==true
%
% Outputs
%   cluster - structure having the following fields
%       .nLocs = number of localizations in a cluster.
%       .area = summed area of all voronoi polygons comprising a cluster.
%       .center = [x,y] position of the cluster's center of mass weighted 
%           by the density of each voronoi polygon in a cluster.
%           Center = <rho*X> of particles where rho = 1/area.
%           Units are the same as input xy array.
%       .Locs = cell matrix containing the original localizations contained
%           within the cluster in columns [x,y] in columns 1,2, zero-rank
%           area in column 3, first-rank area in column 4, index values
%           column 5.
%       .ellipticity = ellipticity ratio for the cluster
%       .HullArea = Area contained in the convex hull drawn around the
%           localizations comprising a cluster; this area is larger than
%           the summed area in the .areas field
%       .NND = Nearest neighbor distance (column 1) in units of input xy
%           array, and the cluster index of the closest neighbor (column 2), 
%           mean distance to all neighbors (column 3), median distance to
%           all neighbors (column 4); neighbors found via triangulation.
%   fig_h - handle to the figure generated if shoIM == true
%   xy - input localization coordinates with zero & first order Voronoi
%       areas in columns 3 & 4, respectively
%   thresh - if input ~isempty(thresh), then the output thresh corresponds
%       to the manually selected Voronoi area threshold
%   repidx - if requested, variable returns the indicies of repeated 
%       localizations output by VoronoiAreas.m iff size(xy,2)==2,
%       otherwise an empty matrix is returned
%
% Function Call
%   [cluster, fig_h, xy, thresh, repidx] = clusterVoronoi(xy,thresh,...
%                     neighborList,minNLoc,convFact,showIM,mask,figTitle)

function [cluster, fig_h, xy, areaThresh, repidx] = clusterVoronoi(xy,areaThresh,...
                    neighborList,minNLoc,convFact,showIM,mask,figTitle)

% initialize in case it's called but not calculated
repidx = []; 

% check format of input [x,y] to see if Voronoi poly's already calc'd
if size(xy,2) < 2 && size(xy,1) < 2
    error('unknown [x,y] localization list format')
elseif size(xy,2) == 2
    if ~exist('neighborList','var') || isempty(neighborList) || nargout==5
        [xy,~,~,repidx,neighborList] = VoronoiAreas(xy,true);
    else
        xy = VoronoiAreas(xy,false);
    end
end
% extract the Voronoi areas
Vareas = xy(:,3);
% extract number of localizations
nPts = length(xy);
% check for threshold input
if ~exist('areaThresh','var') || isempty(areaThresh)
    disp('No area threshold input, select now on plot...')
    VareasPlot = Vareas(~isnan(Vareas));

areaLim = interpercentilerange(VareasPlot,[ .99 .995]);
areaLim = areaLim(1);
[counts, centers] = histcounts(VareasPlot(VareasPlot<areaLim), 'BinMethod','fd');
stp = centers(2)-centers(1);
centers = stp*0.5+centers(1:end-1);
    
    base = [min(centers)*0.05 min(counts)*0.05];
    ylimTop = ceil(1.05*max(counts));
    accept = false;
    while ~accept
        f1 = figure;
        plot(centers,counts,'b-')
        xlabel('Voronoi Polygon Area')
        ylabel('Count')
        title('Click once to set the threshold area')
        axis([0 centers(end), 0 ylimTop])
        [x_ui,~] = ginput(1);
        areaThresh = x_ui(1);    
        hold on
        box on
        area([base(1),areaThresh],[1,1]*ylimTop-base(2),base(2),...intersection(2),...
            'LineStyle','none',...
            'FaceColor',[175 238 238]/255)
        plot(centers,counts,'b-')
        commandwindow
        accept = input(['Area threshold selected to be: ' num2str(areaThresh,'%.4g') ', do you accept? (y/n) ==> '],'s');
        if ~strcmpi(accept,'y')
            accept = false;
        end
        close(f1)
    end
end
    
% check for input conversion factors to go from pixels to nm
if ~exist('convFact','var') || isempty(convFact)
    pix2nm = 1;
    p = 1;
elseif length(convFact)==2
    pix2nm = convFact(1);
    p = convFact(2);
else
    pix2nm = convFact(1);
    p = 1;
end
% check for input neighbor list
if ~exist('neighborList','var') || isempty(neighborList)
    disp('Neighbor list for Voronoi Polygons not included, calculating now...')
    [~,~,~,~,neighborList] = VoronoiAreas(xy,true);
end

% check for the minimum # localizations/cluster
if ~exist('minNLoc','var') || isempty(minNLoc)
    minNLoc = 5;
    disp(['MINIMUM NUMBER OF LOCALIZATIONS PER CLUSTER SET TO DEFAULT = ' num2str(minNLoc)])
end

%% check to see if user wants to plot the mask with localizations
if ~exist('showIM','var') || isempty(showIM)
    showIM = false;
end
if ~exist('figTitle','var') || isempty(figTitle)
    figTitle = '';
else
    if strcmp(figTitle(end-3),'.') % then input is full filename
        figTitle = figTitle(1:end-4); % remove extension
    end
    % remove folder separators if present
    idx = strfind(figTitle,filesep);
    if ~isempty(idx)
        figTitle = figTitle(idx(end)+1:end);
    end
    % now replace underscores with spaces
    figTitle(figTitle=='_') = ' ';
    figTitle = [figTitle '; Threshold = ' num2str(areaThresh,'%.3g')];
end
%% there had better be a mask if they want to plot it
if (~exist('mask','var') || isempty(mask)) && showIM
    disp(['No mask input, calculating now using conversions: pix2nm = ' ...
        num2str(pix2nm) ', final pixel size = ' num2str(p) ' nm'])
    [mask,~] = Locs2Mask( xy(:,1:2));
    mask = imresize(mask,pix2nm/p);
end

%% Begin Algorithm

% define the localizations above the threshold that are to be kept
kppt = Vareas <= areaThresh; % X(:,3)<=thresh; % nanometers

% these are the [x,y] values that will be plotted
xypt = xy(:,1:2).*pix2nm./p; % X(:,1:2)./p;
% check to make sure the mask is big enough for plotting
if showIM && ((max(xypt(:,1)) > size(mask,2)) || (max(xypt(:,2)) > size(mask,1)))
    % then the mask is too small, so re-size it
    mask = imresize(mask,pix2nm/p);
    disp('mask resized for input conversion factors')
end

fig_h = [];
if showIM
    fig_h = figure('units','pixels','Name','Voronoi Clusters');
    disp('Figure will resize when calculation completes, please be patient...')
    imshow(mask)
    fig_ax = gca;
    hold on
    % plot voronoi diagram
    hv = voronoi(xypt(:,1),xypt(:,2),'g');
    % change color of voronoi diag
    for i=1:2
        set(hv(i),'Color',[0.9290    0.6940    0.1250])%[0 0.6 0])
    end
    % blue - plot all localizations having high density
    plot(xypt(kppt,1),xypt(kppt,2),'bs','MarkerSize',3)
    [y,x] = find(mask);
    lm = [min(x) max(x) min(y) max(y)];
    if lm(1) > 3, lm(1) = lm(1)-2; end
    if lm(2) < size(mask,2)-2, lm(2) = lm(2)+2; end
    if lm(3) > 3, lm(3) = lm(3)-2; end
    if lm(4) < size(mask,1)-2, lm(4) = lm(4)+2; end
    set(fig_h,'position',[501 125 918 841])
    axis(fig_ax,lm)
    set(fig_ax,'units','normalized')
    set(fig_ax,'position',[0.0928 0.0777 0.8126 0.8870])
    % update the title if necessary
    title(fig_ax,figTitle)
end

%% get cluster
cluster = getCluster(xy,p,pix2nm,areaThresh,Vareas,xypt,kppt,neighborList,nPts,minNLoc,showIM);
%% find NND for each cluster
cluster = NND_Cluster(cluster);

% future expansion idea: use the external-most verticies of a cluster to
% determine the nearest cluster, then calculate the NND between centers

end % of function

%% Get clusters
function cluster = getCluster(xy,p,pix2nm,areaThresh,Vareas,xypt,kppt,neighborList,nPts,minNLoc,showIM)
c = 0;
cluster = [];
usedPt.nClusters = zeros(nPts,1);
usedPt.cluster = cell(nPts,1);
for pt = 1:nPts
    if kppt(pt) && ~usedPt.nClusters(pt)
        % get the neighbors surrounding the seed-point
        idxNebHi = neighborList{pt,1}(kppt(neighborList{pt,1}));
        if ~isempty(idxNebHi)
            nNbHi = 0;
            sz = size(idxNebHi,1);
            %% find all connected neighbors above thresh
            while sz ~= nNbHi
                %%
                nNbHi = length(idxNebHi);
                tmp = cell(nNbHi,1);
                for n = 1:nNbHi
                    tmp{n} = neighborList{idxNebHi(n),1};
                end
                idxNebHi2 = unique(cell2mat(tmp));
                if ~any(builtin('_ismemberhelper',idxNebHi2,idxNebHi))
                    idxNebHi = sort([idxNebHi2;idxNebHi]);
                else
                    idxNebHi = idxNebHi2;
                end
                idxNebHi = idxNebHi(kppt(idxNebHi));
                sz = size(idxNebHi,1);
            end
        else
            idxNebHi = pt;
        end
        %%
        nLocCluster = length(idxNebHi);
        if nLocCluster >= minNLoc
            %% this implies a core is decided, loop over core points and check for
            % nearest neighbors outside of the cluster
            ExtraCandidates = [];
            for jj = 1:nLocCluster
                localIndex = idxNebHi(jj);
                locNeighb = neighborList{localIndex,1};
                candNeighb = locNeighb(~builtin('_ismemberhelper',locNeighb,idxNebHi));
                keepCand = ~kppt(candNeighb); % only take radii that fail the threshold
                keepNeighb = candNeighb(keepCand);
                distX = xy(localIndex,1)-xy(keepNeighb,1);
                distY = xy(localIndex,2)-xy(keepNeighb,2);
                distR2 = distX.^2 + distY.^2;
                ExtraCandidates = [ExtraCandidates; keepNeighb distR2];
            end
            % find threshold based off a weighted distance from core points
            uniqueCandidates = unique(ExtraCandidates(:,1));
            validCandidates = false(length(uniqueCandidates),1);
            for jj = 1:length(uniqueCandidates)
                uniqueTemp = uniqueCandidates(jj);
                R2vector = ExtraCandidates(ExtraCandidates(:,1)==uniqueTemp,2);
                if ~isempty(R2vector)
                    if pi*min(R2vector) < areaThresh
                        validCandidates(jj)=true;
                    end
                end
            end
            fullIdxNeb = [idxNebHi; uniqueCandidates(validCandidates)]; % add close non-core points
            newLocCluster = length(fullIdxNeb);
            %% extract cluster data
            c = c+1;
            % mark all of the indicies idxNebHi as used in a cluster
            for i = 1:newLocCluster %nLocCluster
                usedPt.cluster{fullIdxNeb(i)} = [c;usedPt.cluster{fullIdxNeb(i)}];
                usedPt.nClusters(fullIdxNeb(i)) = length( usedPt.cluster(fullIdxNeb(i)) );
            end
            % 1) number of localizations in the cluster
            cluster.nLocs(c,1) = newLocCluster; %nLocCluster;
            % 2) total area
            nebAreas = Vareas(fullIdxNeb);
            nebAreas(nebAreas>areaThresh) = areaThresh; % we truncate area size of non-core points
            nebAreas(isnan(nebAreas)) = areaThresh; % NaN values are set to threshold (non-core point issue)
            cluster.areas(c,1) = sum(nebAreas);
            % 3) center of cluster
            % method 3: <rho*X> of particles where rho = 1/area
            rho = 1./nebAreas;
            rXbar = [sum(rho.*xypt(fullIdxNeb,1))/sum(rho),sum(rho.*xypt(fullIdxNeb,2))/sum(rho)];
            % cluster centroid found by weighting each localization by it's voronoi
            % density
            cluster.center(c,:) = rXbar*p/pix2nm;
            % 4) listing of which localizations are in the cluster
            cluster.Locs{c,1} = [xy(fullIdxNeb,:), fullIdxNeb];
            % 5) send localizations for calculation of cluster ellipticity
            cluster.ellipticity(c,1) = ...
                clusterEllipticity(cluster.Locs{c,1}(:,1:2),cluster.center(c,1:2));
            % find edges of the identified cluster
            % 2) Using the localizations comprising the cluster
            clustLoc = xypt(fullIdxNeb,:); % using localizations for perimeter
            %dtclust = delaunayTriangulation( clustLoc );
            % Convex Hull
            %[ConHull, clustVol] = convexHull(dtclust);
            %cluster.HullArea(c,1) = clustVol;

            % using boundary with spectral value of 0.5, PKR:
            [K, clustVol] = boundary(clustLoc);
            cluster.HullArea(c,1) = clustVol;
            % optionally plot Voronoi Diagram
            if showIM
                plot(rXbar(1),rXbar(2),'o','MarkerSize',8,'MarkerFaceColor','k')
                plot(clustLoc(K,1),clustLoc(K,2),'r-') % using localizations from boundary
                %plot(clustLoc(ConHull,1),clustLoc(ConHull,2),'r-') % using localizations for convHull
                % plot all points that are in the cluster
                plot(xypt(fullIdxNeb,1),xypt(fullIdxNeb,2),'c*','MarkerSize',3)
            end
        end
    end
end
if isempty(cluster)
    cluster.nLocs = 0;
    cluster.areas = 0;
    cluster.center = [];
    cluster.Locs = [];
    cluster.ellipticity = [];
    cluster.HullArea = 0;
end
end

%% Cluster NND function
function clusterOut = NND_Cluster(clusterIn)
nClust = length(clusterIn.nLocs);
if nClust>2
    DTctr = delaunayTriangulation(clusterIn.center);
    NL = findVoronoiNeighbors(DTctr,false);
end
clusterOut = clusterIn;
switch nClust
    case num2cell(0:1)
        clusterOut.NND = [-1, -1, -1, -1];
    case 2
        for c = 1:nClust
            % set the NND as the nearest centroid-to-centroid distance
            r = sqrt( ...
                (clusterOut.center(c,1)-clusterOut.center(:,1)).^2 + ...
                (clusterOut.center(c,2)-clusterOut.center(:,2)).^2 ...
                );
            [minNND,idx_r] = min( r(r>0) );
            clusterOut.NND(c,:) = [minNND, idx_r, -1, -1];
        end
    otherwise
        for c = 1:nClust
            % alternative: performa  delaunay Triangulation, then find the neighbors
            % and take an average NND over them
            r_neighbors = sqrt( ...
                (clusterOut.center(c,1)-clusterOut.center(NL{c,1},1)).^2 + ...
                (clusterOut.center(c,2)-clusterOut.center(NL{c,1},2)).^2 ...
                );
            [minNND,idx_r] = min( r_neighbors );
            % note, this minimum result gives a double-counting since each pair of
            % closest clusters will have the same NND
            % [minimum NND, index of closest cluster, mean & median distance  to neighbors]
            clusterOut.NND(c,:) = [minNND, NL{c,1}(idx_r), mean(r_neighbors), median(r_neighbors)]; 
        end
end

end
