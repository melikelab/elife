clear
% clc

%% hard code options - I set things here
pix2nm=116;
finalpix=20;
minNLoc=5;
Allthresholds=0.013; % hard coded threshold
tree = uigetfile_n_dir(); % select multiple sub-directories

%% start directory for loop
for kk = 1:length(tree)
    files=dir([tree{kk} '/*.bin']); 
    nfiles = size(files,1);
    
    %%%%%% begin %%%%%%
    for f = 1:nfiles
        % load localization list
        LL = Insight3( fullfile( files(f).folder,files(f).name) );
        % get coordinates in units of pixels
        xy = LL.getXYcorr;
        % make directory for saving data
        [fpath,fname,~] = fileparts( LL.filename );
        spath = [fpath filesep 'VoronoiCluster-0.013'];
        if ~exist(spath,'dir')
            mkdir(spath)
        end
        badIndex = [];
        [xy,~,~,repidx,neighborList] = VoronoiAreas(xy,true);
        % write a new molecule list assigning localizations to channels 0-9
        % according to their Voronoi area
        savefile = fullfile(spath, [fname '_manualVorSegData.mat']);
        save(savefile,...
            'minNLoc','pix2nm','finalpix',... inputs
            'xy', 'Allthresholds', 'repidx', ... outputs
            'savefile')
        writeIterativeVoronoiSegmentedLL(LL,xy,repidx,Allthresholds,spath)
        % send localizations for clustering
        cluster = cell(1);
        th=1;
        try
            cluster{th} = ...
                clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc);
        catch ME
            msg=('All coordinates are on the same line');
            causeException = MException('DELAUNAY:triangulationEmpty',msg);
            ME = addCause(ME,causeException);
            cluster{th} = [];
        end
        if ~isempty(cluster{th}) && ~isscalar(cluster{th}.nLocs) && ~isscalar(cluster{th}.NND(:,1))
            % save the clusters in a localization list for visualization
            writeVoronoiClusteredLL(LL,cluster{th},Allthresholds(th,1), repidx)
        else
            badIndex = [badIndex, th];
        end
        % save cluster info
        save(savefile,'cluster','-append')
        close all;
    end % file loop
end
%% end directory for loop
