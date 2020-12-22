% Input
%   DT - delaunayTriangulation object for a set of [x,y] points
%   kppt - logical array indicating the [x,y] points that should have their
%       neighbors identified;  length(kppt) = size(DT.Points,1)
%   dispStrOut - logical value whether to write updates to command window
%
% Output
%   neighborList - a 2-column cell matrix where each row corresponds to one
%       input position in [x,y]. The first column contains a list of
%       neighboring positions for each input as determined by
%       triangulation. The second column is a count for the number of
%       neighbors.  
%      The rows are empty for those points where kppt(point) == false.

function neighborList = findVoronoiNeighbors(DT, kppt)

if ~isa(DT,'delaunayTriangulation')
    error('Input must be a delaunayTriangulation object')
end
nPoints = size(DT.Points,1);

if ~exist('kppt','var') || isempty(kppt)
    kppt = true(nPoints,1);
elseif ~islogical(kppt) && length(kppt)>1
    unq_kppt = unique(kppt);
    if length(unq_kppt) == 2 && unq_kppt(1) == 0 && unq_kppt(2) == 1
        kppt = logical(kppt);
    else
        error('Input point-selection array is of an unknown format')
    end
elseif length(kppt)==1 % then second input is intended as third
    kppt = true(nPoints,1);
end

%% initialize output
neighborList = cell(nPoints,2);

% using the sortrows function to identify neighbors is much faster than
% using the 'find' function
[ConSort,idx] = sortrows(DT.ConnectivityList(:));
nCon = length(ConSort);
szConMtx = size(DT.ConnectivityList);
prevIdx = 1;
strCR = '';
%%
for pt=1:nPoints
    
    if kppt(pt)
        % find 1st rank neighbors using delaunay triangulation
        nextIdx = prevIdx+1;
        while nextIdx<=nCon && pt == ConSort(nextIdx)
            nextIdx = nextIdx+1;
        end
        nextIdx = nextIdx-1; % remove increment since test was actually for ConSort(nextIdx-1)
        [r,~]=ind2sub(szConMtx,idx(prevIdx:1:nextIdx));
        % determine the indicies for the point in question & neighbors
        %pt_neighbor = unique(DT.ConnectivityList([r;r+szConMtx(1);r+2*szConMtx(1)]));
        pt_neighbor = sort(DT.ConnectivityList([r;r+szConMtx(1);r+2*szConMtx(1)]));
        pt_neighbor(pt_neighbor((1:end-1)')==pt_neighbor((2:end)')) = [];
        % remove the point in question and save a listing of neighbors for each point
        neighborList{pt,1} = pt_neighbor(pt_neighbor~=pt);
        % save the number of neighbors
        neighborList{pt,2} = nextIdx-prevIdx+1;
    else
        nextIdx = find(ConSort==pt,1,'last');
    end
        prevIdx = nextIdx+1; % start next search from subsequent index
end

end % of function