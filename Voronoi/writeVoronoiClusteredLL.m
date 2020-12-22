% Function that takes an input localization list and cluster information
% from voronoi area thresholding and will write it to a new localization
% list where different clusters are assigned to channels 1-9 and centroid
% to channel 0
%
% Input
%   LL - original Insight3 localization list
%   cluster - clustering information that is the output of the program
%       clusterVoronoi.m
%   thresh (optional) - accepted inputs include:
%       1) the value (as a double) of the threshold used to segregate
%       localizations into the identified clusters that will be included 
%       into the filename of the newly saved localization list
%       2) a string of user's choice that will be appended to the new
%       localization list
%   repidx = output from VoronoiAreas.m that is empty if all localizations
%       are unique, or contains data about repeated localizations that are
%       removed from clustering here.
%
% Output
%   No variables are returned
%   A new localization list is written in the same folder as the input
%   localization list. The filename is appended with the string:
%       '_VornClust' and the input 'thresh' value or string

function writeVoronoiClusteredLL(LL,cluster,thresh, repidx)
%%

if ~isa(LL,'Insight3')
    error('First input must be an Insight3 localization list')
end

if ~isa(cluster,'struct') && ~isfield(cluster,'Locs')
    error('Second input must be a cluster variable output from clusterVoronoi.m')
end


if ~exist('thresh','var') || isempty(thresh)
    thresh = '';
elseif isa(thresh,'double')
    thresh = sprintf('_Th%.2g',thresh);
    thresh( thresh=='.' ) = ',';
elseif isa(thresh,'string')
    disp(['Input file appendix = ' thresh])
else
    error('Third "thresh" input is not a recognized format')
end

nclust = size(cluster.center,1);

dat = LL.getData();

if ~isempty( repidx )
    % some localizations are non-unique (repeated) and must be removed
    dat = dat( repidx.unqIdx, : );
end
chidx = LL.getColumnIndex('channel');
% delete all channel information to permit data check at the end
dat(:,chidx) = nan;

% define pseudo localization list for cluster centers 
% channel zero default
centers = zeros(nclust,size(dat,2));
xyidx = [LL.getColumnIndex('x'), LL.getColumnIndex('y'), ...
         LL.getColumnIndex('xc'),LL.getColumnIndex('yc')];
setCols = {'height',100;...
           'area',10000;...
           'width',300;...
           'aspect',1;...
           'intensity',10000;...
           'fitIterations',1;...
           'frame',1;...
           'trackLength',1;...
           'link',-1};
for i = 1:size(setCols,1)
    idx = LL.getColumnIndex( setCols{i,1} );
    centers(:,idx) = setCols{i,2};
end
centers(:,xyidx) = [cluster.center, cluster.center];

% the localization list will consist of 'm'=nclust clusters, 
% cluster 'c' has nlocs(c) number of localizations
% order localization list according to:
    % row_j = Cluster 'c' centroid position
    % row_j+1 ... row_j+nlocs(c) = original localization information
    
% initialize new data array
LocsInClusters = cell2mat(cluster.Locs);
% future idea:
% make list with the localizations in clusters belonging to Ch 1-9 and the
% localizations not in clusters, i.e. below the threshold, in channel 0
% in this case the cluster centroid is not saved.
datnew = nan( size(LocsInClusters,1),size(dat,2) );
j = 1;
ch = 1;
% %%
for c = 1:nclust
%     %%
%     c = 3;% j = 1; ch = 1;
    % extract localization indicies
    locIdx = cluster.Locs{c,1}(:,5);
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
LLnew.setFilename( [LL.filename(1:end-4) '_VornClust' thresh LL.filename(end-3:end)] );
LLnew.forceFileOverwrite(true);
LLnew.write();


end % of function