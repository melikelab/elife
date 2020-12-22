% colocPoints2ClustersBatch
% Function to perform clustering on a batch of files in a reference and
% colocalization channel given two directories

[fileR,pathR] = uigetfile('.mat','Select Reference Channel Data','MultiSelect','on');
[fileG,pathG] = uigetfile('.mat','Select Colocalization Channel Data','MultiSelect','on');
saveDir = uigetdir(pathG,'Select Save Out Folder');

if ~iscell(fileR)
    tempF{1} = fileR;
    fileR = tempF;
end

if ~iscell(fileG)
    tempF{1} = fileG;
    fileG = tempF;
end

% threshold parameters for determining full/iso/mixed clusters
FullAlpha = 0.4; % 1 = 100%
IsoBeta = 0.39999; % 0.05 = 5%
% Isolated Cluster <= IsoBeta Parameter < Mixed Cluster < FullAlpha Parameter <= Full Cluster

ExtendedRadius = 30/116; % radius of colocalization search in pixels

for k = 1:length (fileR)
    referenceData = fullfile(pathR,fileR{k});
    colocalizationData = fullfile(pathG,fileG{k});
    refS = load(referenceData);
    colocS = load(colocalizationData);
    
    cluster0 = refS.cluster{1};
    cluster1 = colocS.cluster{1};
    Localizations1 = colocS.xy;
    %% get a list of reference boundaries from reference cluster
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
    
    ColocValids = false(size(Localizations1,1),1);
    
    %% Determine points from colocalization data set that fall into reference cluster
    for ii = 1:length(cluster0.Locs)
        ColocValids = ColocValids | inpolygon(Localizations1(:,1),Localizations1(:,2),Out0{ii}(:,1),Out0{ii}(:,2));
    end
    
    CoLocalizations = Localizations1(ColocValids,:);
    IsoLocalizations = Localizations1(~ColocValids,:);
    
    %% determine % of colocalized points in the non-reference clusters
    ColocPercentages = zeros(length(cluster1.Locs),1);
    ColocMembers = find(ColocValids);
    
    % get three cluster lists
    f1n = fieldnames(cluster1)';
    f2n = f1n;
    len = length(f1n); %// number of fields
    for i=1:len
        f2n{2,i} = cell(size(cluster1')); %'// the contents are empty
    end
    FullCluster1 = struct( f2n{:} );
    IsoCluster1 = struct( f2n{:} );
    MixedCluster1 = struct( f2n{:} );
    
    for ii = 1:length(cluster1.Locs)
        tempIdx = [];
        tempValidIndex = [];
        tempNumLocs = cluster1.nLocs(ii);
        tempIdx = cluster1.Locs{ii}(:,5);
        tempValidIndex = builtin('_ismemberhelper',tempIdx,ColocMembers);
        ColocPercentages(ii) = sum(tempValidIndex)/size(tempValidIndex,1);
    end
    
    % change the logic here... so that x <= beta = 0, beta < x < 1 - alpha = 1,
    % x => 1- alpha = 2
    IsoCond = ColocPercentages <= IsoBeta;
    FullCond = ColocPercentages >= FullAlpha;
    MixCond = ~IsoCond & ~FullCond;
    ColCaser = IsoCond + 2*MixCond + 3*FullCond;
    for ii = 1:length(ColocPercentages)
        switch ColCaser(ii)
            case 1
                for jj = 1:numel(f1n)
                    IsoCluster1.(f1n{jj}) = [IsoCluster1.(f1n{jj}); cluster1.(f1n{jj})(ii)];
                end
            case 3
                for jj = 1:numel(f1n)
                    FullCluster1.(f1n{jj}) = [FullCluster1.(f1n{jj}); cluster1.(f1n{jj})(ii)];
                end
            case 2
                for jj = 1:numel(f1n)
                    MixedCluster1.(f1n{jj}) = [MixedCluster1.(f1n{jj}); cluster1.(f1n{jj})(ii)];
                end
        end
    end
    
    
    %% save out colocalization statistics as a .mat file
    saveOutName = [fullfile(saveDir,fileG{k}) '_ColocPoints.mat'];
    save(saveOutName,'cluster0','cluster1','ExtendedRadius',...
        'CoLocalizations','IsoLocalizations','ColocPercentages',...
        'IsoCluster1','FullCluster1','MixedCluster1');
    
    %% generate a .bin file
    % Current configuration:
    % channel1 = reference clusters
    % channel2 = colocalized points from coloc raw data
    % channel3 = isolated points from coloc raw data
    % channel4 = all points from coloc raw data
    % channel5 = totally colocalized clusters in coloc cluster structure
    % channel6 = isolated and colocalized points found in these coloc cluster structures
    % channel7 = totally isolated clusters in coloc cluster structure
    
    %Insight3 data formatting
    chidx = 12;
    xyidx = [1,2,3,4];
    columnValues = [5,100;6,10000;7,300;9,1;11,10000;13,1;14,1;15,1;16,-1]; 
    
    % initialize bin file matrices and save out
    % work on cluster structures first
    % create a cell array of our cluster structures
    clusterCell = {cluster0,FullCluster1,MixedCluster1,IsoCluster1};
    clusterChannels = [1, 5, 6, 7];
    clusterDat = cell(length(clusterCell),1);
    % loop over each cluster structure and add to the data cell array for
    % bin file writing
    for jj = 1:length(clusterCell)
        xyTemp = [];
        for ii = 1:length(clusterCell{jj}.nLocs)
            xyTemp = [xyTemp; clusterCell{jj}.Locs{ii}(:,[1 2])];
        end
        if isempty(xyTemp)
            continue
        end
        clusterDat{jj} = zeros(size(xyTemp,1),18);
        clusterDat{jj}(:,xyidx) = [xyTemp(:,[1,2]),xyTemp(:,[1,2])];
        for ii = 1:size(columnValues,1)
            clusterDat{jj}(:,columnValues(ii,1)) = columnValues(ii,2);
        end
        clusterDat{jj}(:,chidx) = clusterChannels(jj);
    end
    % work on localizations next
    % create a cell array of localization matrices
    LocCell = {CoLocalizations,IsoLocalizations,Localizations1};
    LocChannels = [2,3,4];
    LocDat = cell(length(LocCell),1);
    % loop over each localization matrix and add to the data cell array for
    % bin file writing
    for jj = 1:length(LocCell)
        if isempty(LocCell{jj})
            continue
        end
        LocDat{jj} = zeros(size(LocCell{jj},1),18);
        LocDat{jj}(:,xyidx) = [LocCell{jj}(:,[1,2]),LocCell{jj}(:,[1,2])];
        for ii = 1:size(columnValues,1)
            LocDat{jj}(:,columnValues(ii,1)) = columnValues(ii,2);
        end
        LocDat{jj}(:,chidx) = LocChannels(jj);
    end
    
    % add all clusters to one matrix
    finalDat = [];
    for jj = 1:length(clusterDat)
        finalDat = [finalDat;clusterDat{jj}];
    end
    for jj = 1:length(LocDat)
        finalDat = [finalDat;LocDat{jj}];
    end
    
    % write a new localization list
    LLnew = Insight3();
    LLnew.setData( finalDat );
    LLnew.setFilename( [fullfile(saveDir, fileG{k}) '_ColocPoints.bin'] );
    LLnew.forceFileOverwrite(true);
    LLnew.write();
end