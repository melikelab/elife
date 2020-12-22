% Function for performing 2D Voronoi Tesellation and calculating the zero 
% and first rank area for each [x,y] input point.
%   J.Otterstrom MATLAB 2016a
%
% INPUTS
%   x (required) - either a vector of 'x' positions or a matrix of '[x,y]' 
%       positions
%   y (optional) - a vector of 'y' positions; include this input ONLY if x 
%       is a vector
%   rankAreas (optional) - a logical true (1) or false (0) input to either
%       calculate only the zero rank areas (false input) or also the first
%       rank areas (true input). Default value is false (0).
%   showText (optional) - logical true (1) or false (0) input to determine
%       whether the user is updated regarding the status of the Voronoi
%       area calculations
%
% OUTPUT
%   X - either a 3-column (rankAreas=false) or 4-column(rankAreas=true)
%       matrix with columns corresponding to:
%       column 1: input 'x' vector, or first column of x matrix
%       column 2: input 'y' vector, or second column of x matrix
%       column 3: zero rank area per point
%       column 4: first rank area per point iff rankAreas=true
%           NOTE: All localizations in the returned 'X' variable are
%           unique, hence there can be fewer localizations returned than
%           input.  See repidx variable for info about dealing with
%           repetitions.
%   DT - Delaunay Triangulation object output by MATLAB function:
%       delaunayTriangulation()
%   neighborList - a 2-column cell matrix where each row corresponds to one
%       input position in [x,y]. The first column contains a list of
%       neighboring positions for each input as determined by
%       triangulation. The second column is a count for the number of
%       neighbors.
%   VorDat - 1x2 cell array {V, C} containing the Voronoi verticies as the 
%       first element and the Voronoi regions are the second element, both 
%       being returned from the voronoiDiagram method for 
%       delaunayTriangulation objects
%   repidx - output reporting any non-unique localizations
%           If all localizations are unique, it is an empty variable {}
%           If some localizations are repeated, then it is a structure 
%               variable with two fields:
%       .info: mx4 cell array reporting the indicies of m number of 
%           repeated localizations in the input xy localization list. 
%          column 1 = indicies of repeated localizations
%           column 2 = [x,y] values of the repeated localizations (for checking)
%           column 3 = difference between localizations reported in column 2, 
%               all elements should be [0,0], otherwise there's a problem
%          column 4 = difference between the indicies of repeated localizations
%               there should not be any values of zero
%       .origIdx = vector of indicies pointing from the new (output) xy
%           list to the original localization order.  
%       .unqIdx = vector of indicies pointing from the original (input) xy
%           list to the output localization order.
%           In this case, the following calls can be used to obtain
%           the original order:
%               X_input( :, 1:2) == X_output( repidx.origIdx, 1:2 );
%                   used in writeIterativeVoronoiSegmentedLL.m to include
%                   all localizations for channel assignment based on
%                   Voronoi polyon area
%           Or the new order:
%               X_input( repidx.unqIdx, 1:2 ) = X_output( :, 1:2 );
%                   used in writeVoronoiClusteredLL.m to truncate the 
%                   localization list and exclude repetitions from 
%                   identified clusteres, which are later written into 
%                   separate .bin files
%
% to call function use commands:
%   Full function
% [X,DT,VorDat,repidx,neighborList] = VoronoiAreas(x,y,rankAreas,showText)
%   abreviated versions
% X = VoronoiAreas(x,y,rankAreas) where x and y are vectors and rankAreas 
%   is logical true (1) or false (0) 
% X = VoronoiAreas(x,rankAreas) where x is a matrix and rankArea is logical
%   1 or 0
%   equivalently: X = VoronoiAreas(x,[],rankAreas)
% X = VoronoiAreas(x,y) where x and y are vectors
% X = VoronoiAreas(x) where x is a matrix
%
% Finished: 160913
% update: 160915, included neighborList output
% update: 161010, use delaunayTriangulation function output for Voronoi
%       calculation
% update: 161017, included showText variable
% update: 161021, include check for unique localizations and return the
%       repeated indicies in repidx variable
% update: 171116, new format for 'repidx' output as a structure
%       now includes .info, .unqIdx and .origIdx fields, which areused in 
%       the functions 'writeIterativeVoronoiSegmentedLL.m' and
%       'writeVoronoiClusteredLL.m'



function [X,DT,VorDat,repidx,neighborList] = VoronoiAreas(x,y,rankAreas)

%% Process inputs
if nargin < 2 || isempty(y) || length(x)~=length(y)% (length(y)==1 && (~exist('rankAreas','var') || isempty(rankAreas)) )
    if exist('y','var') && ~isempty(y)
        if length(y)==1
            rankAreas = y;
            clear y
        else
            error('unrecognized second input')
        end
    end
    % input assumed to be a 2-column matrix of [x,y] values
    if size(x,2)>size(x,1)%~iscolumn(x) % then the data is along columns, so flip it
        rot = 1;
        X = transpose(x);
    else
        rot = 0;
        X = x;
    end
else
    if ~iscolumn(x)
        x = transpose(x); 
    end
    if ~iscolumn(y) 
        y = transpose(y); 
    end
    X = [x,y];
    rot = 0;
end

if nargout > 4 && (~exist('rankAreas','var') || isempty(rankAreas) || rankAreas~=1) % then they want the neighbor list
    rankAreas = true;
    disp('You request the neighbor list, so the 1st rank areas will also be calculated')
elseif ~exist('rankAreas','var') || isempty(rankAreas)
    rankAreas = false;
end

%% Begin algorithm
nPoints = size(X,1);

tic

% Estimate the amount of time needed for the calculation
if rankAreas
    X = [X,nan(nPoints,2)];
    estT = polyval([1.0015 -3.5355],log10(nPoints));
    % initialize output
    neighborList = cell(nPoints,2);
else
    X = [X,nan(nPoints,1)];
    estT = polyval([0.97777 -3.9502],log10(nPoints));
end
estT = 10^estT;

% need to ensure the positions are unique
[xu,ixy,ixu] = unique(X(:,1:2),'rows','stable');

if length(xu) ~= nPoints
    % find the indicies of repeated localizations
    bincts = histc(ixu,unique(ixu)); %size(find(bincts==2))
    k = find(bincts > 1);
    repidx.info = cell(length(k),4);
    repidx.origIdx = ixu;
    repidx.unqIdx = ixy;
    for j = 1:length(k)
        repidx.info{j,1} = intersect(find(X(:,1)==xu(k(j),1)),find(X(:,2)==xu(k(j),2)));
        repidx.info{j,2} = X(repidx.info{j,1},1:2); % original Locs
        repidx.info{j,3} = diff(repidx.info{j,2}); % difference between [x,y] of repeats, should be [0,0]
        repidx.info{j,4} = diff(repidx.info{j,1}); % index separation
    end
    X = X(ixy,:);
    nPoints = size(X,1);
else
    repidx = {};
end

% perform Voronoi Tessellation via Delaunay Triangulation method
if X~=0 
    %if size(X,1) >=3
DT = delaunayTriangulation(X(:,1:2));

[V,C] = voronoiDiagram(DT);
VorDat = {V, C};
% Calculate the area for each Voronoi polygon
strCR = '';

for pt=1:nPoints
    xt = V(C{pt},1);
    yt = V(C{pt},2);
    X(pt,3) = abs(sum( (xt([2:end 1]) - xt).*(yt([2:end 1]) + yt))*0.5);
end

if rankAreas
    % use DT to find neighboring polynomials and calculate summed 1st rank
    % area
    neighborList = findVoronoiNeighbors(DT, []);
    strCR = '';
    for pt=1:nPoints 
        % if the polygon is closed, calculate first-rank area
        if ~isnan(X(pt,3))
            neighbPts = [pt;neighborList{pt,1}];
            neighborAreas = X(neighbPts,3);
            neighborAreas = neighborAreas(~isnan(neighborAreas));
            X(pt,4) = sum(neighborAreas);
        end
    end
end


if rot, X = transpose(X); end

    %end
end
% subfunction
    function [tin,tmunits] = detTunits(tin)
        if tin<60
            tmunits = ' seconds';
        elseif tin<60*60
            tin=tin/60;
            tmunits = ' minutes';
        elseif tin<60*60*24
            tin=tin/(60*60);
            tmunits = ' hours';
        else
            tin=tin/(60*60*24);
            tmunits = ' days';
        end
    end

end
