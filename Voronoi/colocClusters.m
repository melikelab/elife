% colocClusters.m
% script to roughly determine colocalization via a hard overlap threshold

[fileR,pathR] = uigetfile('.mat','Select Reference Channel Data');
[fileG,pathG] = uigetfile('.mat','Select Colocalization Channel Data');
saveDir = uigetdir(pathG,'Select Save Out Folder');

% threshold to determine colocalization
threshold = 0.60;
ExtendedRadius = 50/116; % radius of colocalization search in pixels 

redData = fullfile(pathR,fileR);
greenData = fullfile(pathG,fileG);

redS = load(redData);
greenS = load(greenData);

redClusters = redS.cluster{1};
greenClusters = greenS.cluster{1};

% flip these variable defs. instead of the code below if you want to do green on
% red or vice versa
cluster0 = redClusters;
cluster1 = greenClusters;

% cluster colocalization/association matrix
%ColocMat = zeros(length(greenClusters.nLocs),length(redClusters.nLocs));
ColocList10 = cell(length(greenClusters.nLocs),1);

% let's do green on red, and then we can switch later
% get a list of red boundaries
Out0 = cell(length(cluster0.Locs),1);
xBounds0 = NaN(length(cluster0.Locs),2);
yBounds0 = NaN(length(cluster0.Locs),2);
% get the boundary of the reference cluster, extend it by the radius
for ii = 1:length(cluster0.Locs)
    xy0 = cluster0.Locs{ii}(:,[1 2]);
    OB = boundary(xy0);
    meanPos = mean(xy0);
    % extend boundary by search radius
    numEl = length(OB);
    dist2mp = zeros(numEl,2);
    for jj = 1:numEl
        dist2mp(jj,:) = xy0(OB(jj),:) - meanPos;
    end
    radii = sqrt(dist2mp(:,1).^2 + dist2mp(:,2).^2);
    theta = atan2(dist2mp(:,2),dist2mp(:,1));
    % get inner boundary points and filter
    extended_radii = radii + ExtendedRadius;
    extended_dist2mp = [extended_radii.*cos(theta), extended_radii.*sin(theta)];
    extendedXY = [extended_dist2mp(:,1)+meanPos(1), extended_dist2mp(:,2)+meanPos(2)];
    
    Out0{ii} = extendedXY;
    xBounds0(ii,:) = [min(xy0(:,1))-ExtendedRadius max(xy0(:,1))+ExtendedRadius];
    yBounds0(ii,:) = [min(xy0(:,2))-ExtendedRadius max(xy0(:,2))+ExtendedRadius];
end

% loop over second set of clusters
    for jj = 1:length(cluster1.Locs)
        xy1 = cluster1.Locs{jj}(:,[1 2]);
        valids = false(size(xy1,1),1);
        colocVec = [];
        % loop over all of the reference cluster boundaries
        for ii = 1:length(cluster0.Locs)
            % skip obviously separated clusters
            if ~any(xy1(:,1) > xBounds0(ii,1) & xy1(:,1) < xBounds0(ii,2))
                continue
            end
            % check y bounds
            if ~any(xy1(:,2) > yBounds0(ii,1) & xy1(:,2) < yBounds0(ii,2))
                continue
            end
            colocVec = [colocVec; ii];
            % add points from cluster1 that fall in reference clusters
            valids = valids | inpolygon(xy1(:,1),xy1(:,2),Out0{ii}(:,1),Out0{ii}(:,2));
        end
        percentage = sum(valids)/numel(valids);
        if percentage > threshold
            ColocList10{jj} = colocVec;
        end
    end

% get a list of the cluster1 clusters that overlap on cluster0 clusters
OverLapListC1 = ~cellfun(@isempty,ColocList10);
% get a list of cluster0 clusters that over lap with colocalized cluster1
% clusters
OverlapC0 = false(length(cluster1.nLocs),1);
for jj = 1:length(OverlapC0)
   if ~isempty(ColocList10{jj})
       OverlapC0(ColocList10{jj}) = true;
   end
end

OverlappedClusters1 = struct;
IsolatedClusters1 = struct;

OverLapListC0 = find(OverlapC0);

% loop over cluster1 and parse overlapped and isolated clusters
f = fieldnames(cluster1);
for jj = 1:length(f)
    OverlappedClusters1.(f{jj}) = cluster1.(f{jj})(OverLapListC1,:);
    IsolatedClusters1.(f{jj}) = cluster1.(f{jj})(~OverLapListC1,:);
end

saveOutName = [fullfile(saveDir,fileG(1:end-4)) '_Coloc.mat'];
save(saveOutName,'OverlappedClusters1','IsolatedClusters1','threshold','ColocList10',...
    'cluster0','cluster1');

% save out a bin file where overlapped clusters are in 1 channel, isolated
% in another channel and the reference cluster are another channel

%Insight3 data formatting
chidx = 12;
xyidx = [1,2,3,4];
columnValues = [5,100;6,10000;7,300;9,1;11,10000;13,1;14,1;15,1;16,-1];

% initialize bin file matrices and save out
% do cluster 0
xyTemp = [];
for ii = 1:length(cluster0.nLocs)
    xyTemp = [xyTemp; cluster0.Locs{ii}(:,[1 2])];
end
dat0 = zeros(size(xyTemp,1),18);
dat0(:,xyidx) = [xyTemp(:,[1,2]),xyTemp(:,[1,2])];
for ii = 1:size(columnValues,1)
    dat0(:,columnValues(ii,1)) = columnValues(ii,2);
end
dat0(:,chidx) = 1;
% do cluster 1 colocalized
xyTemp = [];
if ~isempty(OverlappedClusters1)
    for ii = 1:length(OverlappedClusters1.nLocs)
        xyTemp = [xyTemp; OverlappedClusters1.Locs{ii}(:,[1 2])];
    end
end
dat1Coloc = zeros(size(xyTemp,1),18);
dat1Coloc(:,xyidx) = [xyTemp(:,[1,2]),xyTemp(:,[1,2])];
for ii = 1:size(columnValues,1)
    dat1Coloc(:,columnValues(ii,1)) = columnValues(ii,2);
end
dat1Coloc(:,chidx) = 2;
% do cluster 1 isolated
xyTemp = [];
if ~isempty(IsolatedClusters1)
    for ii = 1:length(IsolatedClusters1.nLocs)
        xyTemp = [xyTemp; IsolatedClusters1.Locs{ii}(:,[1 2])];
    end
end
dat1Iso = zeros(size(xyTemp,1),18);
dat1Iso(:,xyidx) = [xyTemp(:,[1,2]),xyTemp(:,[1,2])];
for ii = 1:size(columnValues,1)
    dat1Iso(:,columnValues(ii,1)) = columnValues(ii,2);
end
dat1Iso(:,chidx) = 3;
% do total green clusters
xyTemp = [];
for ii = 1:length(cluster1.nLocs)
   xyTemp = [xyTemp; cluster1.Locs{ii}(:,[1,2])]; 
end
dat1Total = zeros(size(xyTemp,1),18);
dat1Total(:,xyidx) = [xyTemp(:,[1,2]),xyTemp(:,[1,2])];
for ii = 1:size(columnValues,1)
    dat1Total(:,columnValues(ii,1)) = columnValues(ii,2);
end
dat1Total(:,chidx) = 4;
% add all clusters to one matrix
finalDat = [dat0;dat1Coloc;dat1Iso;dat1Total];

% write a new localization list
LLnew = Insight3();
LLnew.setData( finalDat );
LLnew.setFilename( [fullfile(saveDir, fileG(1:end-4)) '_Colocalized.bin'] );
LLnew.forceFileOverwrite(true);
LLnew.write();