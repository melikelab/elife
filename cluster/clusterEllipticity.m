% Calculates the ellipticity for a single a cluster given the localizations
% comprising it, and, optionally, the center of the cluster.
%
% Inputs
%  xyval - a 2xn row-matrix of n localizations with x values in the first
%    row and y values in the second row.  If the input is a nx2
%    column-matrix, then the transpose is used.  If there are more than 2
%    rows or columns, then only the first two are used.
%  center - (optional) this is the center [x,y] of the cluster comprised by
%    the localizations in xyval.  If no input is given, then it is
%    calculated as the mean(x) and mean(y) positions
%
% Output
%  ellipticity - calculated as sqrt(1-(b/a)^2) where b is the short radius
%    and a is the long radius of an ellipse.
%    see: http://mathworld.wolfram.com/Ellipticity.html

function ellipticity = clusterEllipticity(xyval,center)

%%
% xyval = cluster.Locs{j};
% ctr = [cluster.center(j,1),cluster.center(j,2)];

% find if the matrix orientation
if size(xyval,2) < size(xyval,1)
    % more rows than columns, so use transpose
    xyval = transpose(xyval);
end

if size(xyval,1) > 2 % take just first two rows
    xyval = xyval(1:2,:);
end

if ~exist('center','var') || isempty(center)
    center = mean(xyval,2);
elseif length(center) ~= 2
    error('Must input [x,y] coordinates of cluster center')
end

if min(size(xyval)) == 1
    % input was a single [x,y] corrdinate
    a = 1;
    b = 1;
else
    %% set center of [x;y] at [0;0]
    for r = 1:2
        xyval(r,:) = xyval(r,:)-center(r);
    end
    % get slope, angle and rotation matrix
    m = xyval(2,:)'\xyval(1,:)';
    theta = -atan(m);
    R = [cos(theta),-sin(theta); sin(theta), cos(theta)];
    % % define reference line & rotate it
    % ln = [min(xyval(1,:)),max(xyval(1,:))];
    % ln(2,:) = ln(1,:).*m;
    % ln2 = R*ln
    % rotate localizations
    xyval2 = R*xyval;
    % calculate the ellpticity with rotated points
    xsz = max(xyval2(1,:))-min(xyval2(1,:));
    ysz = max(xyval2(2,:))-min(xyval2(2,:));
    if xsz < ysz,
        b = xsz; a=ysz;
    else
        b = ysz; a=xsz;
    end
    
end

ellipticity = sqrt( 1-(b/a)^2 );


end