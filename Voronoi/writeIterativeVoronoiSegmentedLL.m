% Function writes a new localization list for Insight3 where localizations
% are assigned to a channel depending on the Voronoi area they occupy,
% where segmentation is obtained using the function
%   iterativeVoronoiSegmentation.m
% 
% Inputs
%   LL - Insight3 object 
%   xy - array of [x,y] coordinates in columns [1,2], the zero rank Voronoi 
%       areas in column 3 and, optionally, first rank Voronoi areas in
%       column 4
%   repidx - structure variable output by VoronoiAreas.m
%   thresh - vector of Voronoi area threshold values for segmenting the 
%       localization list
%   spath - (optional) path where the localization lists should be saved
%
% Function call (first four inputs required)
%   writeIterativeVoronoiSegmentedLL(LL,xy,repidx,thresh,spath)
%
% update: Nov 16, 2017, J.O. MATLAB 2016a
%   support for new repidx structure format output from VoronoiAreas.m
%   included back-compatibility for previous repidx cell format
%   performed to remove error of repeated localizations

function writeIterativeVoronoiSegmentedLL(LL,xy,repidx,thresh,spath)

[fpath,fname,ext] = fileparts( LL.filename );
if ~exist('spath','var') || exist('spath','dir')
    spath = fpath;
end
saveName = fullfile( spath, [fname '_iterVorSeg' ext] );
%%
% load localization list & extract channel information
dat = LL.getData();
ch = LL.getColumnIndex('channel');
% check for legacy repidx
if ~isempty(repidx) && ~isstruct( repidx ) && iscell( repidx )
    repidxIn = repidx;
    clear repidx
    repidx.info = repidxIn;
    repidx.origIdx = ( 1:LL.numMolecules )';
    for j = 1:size(repidx.info,1)
        repidx.origIdx( repidx.info{j,1} ) = repidx.info{j,1}(1);
    end
end
% if there were repetitions detected previously, update xy to include them
if ~isempty( repidx )
    xy = xy( repidx.origIdx,:);
end
% delete all channel information
dat(:,ch) = nan;
%% find localizations on the edge with infinite Voroni area
idxSet = isnan(xy(:,3)); % channel 0
dat(idxSet,ch) = 0;
% %%
for t = 1:size(thresh,1)
    % find localizations above the threshold
    idx = xy(:,3) > thresh(t,1);
    % avoid previously chosen points
    idx(idxSet) = false;
    % reset the channel assignment 
    dat(idx,ch) = t-1;
    % update the list of utilized localizations
    idxSet(idx) = 1;
end

dat(~idxSet,ch) = t;

if ~isempty( find(isnan(dat(:,ch)), 1) )
    error('  Uh oh, some localizations didn''t get their channel reset')
    % problem arises when the duplicates are removed for the Voronoi xy
    % calculation, but remain in the original .bin list
end
%%
% save the thresholded .bin list
LL.setFilename( saveName );
LL.forceFileOverwrite(true);
LL.setData( dat );
LL.write;

end % of function