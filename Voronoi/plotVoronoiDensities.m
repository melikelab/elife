% plotVoronoiDensities.m
% function that receives the output of calling VoronoiAreas.m and will plot
% the Vornoi Tessellation diagram color-coding each cell according to the
% area it occupies.
% The area can be converted from pixels to nm by an optional input
%
% INPUTS
%   x - 3 or 4 column list in the format 
%       [x_coord, y_coord, area_zeroRank]
%       that is output by the VoronoiAreas.m function
%  thresholds = [minArea, maxArea]
%  VorDat -- voronoi class from MATLAB, generated if empty
%  pix2nm -- pixel to nm conversion factor, triggers scale par

% Modified by PKR, UPENN March 2019, everything unnecessary removed.

function fig_h = plotVoronoiDensities(x, thresholds, pix2nm, VorDat, barWidth, range)

if nargin < 5
    barWidth = 100; % width of scale bar, hard coded for now - PKR
end

if nargin < 6
    range = 256; % range of colors to render
end
    
if nargin < 2 || isempty(thresholds)
    thresholds = [0,inf];
end

% ensure the input data has the needed Voronoi information
if size(x,2) < 3 || nargin < 4 || ~iscell(VorDat)
    [x, ~,VorDat] = VoronoiAreas(x(:,1:2));
end

% set the pixel conversion and scale bar if input
if nargin < 3 || isempty(pix2nm)
    pix2nm = 1;
    units = 'AU^2';%[];
    areaFmt = '%4.4f';
    includeScaleBar = false;
else
    % PKR -- I'm assuming the user is an adult, put in the right units.
    units = 'nm^2';
    areaFmt = '%4.0f';
    includeScaleBar = true;    
end

%%
% find the max & min [x,y] values to set the axis limits
mnx = min(x(:,1));
mxx = max(x(:,1));
mny = min(x(:,2));
mxy = max(x(:,2));
nloc = size(x,1);

if includeScaleBar
    % calculate the bar width and check its appropriateness now, 
    %   since plotting the voronoi areas can be slow
    w = barWidth/pix2nm;
    % do not allow width to be larger than 40% of image size
    while w > 0.4*(mxx-mnx) % (1-2*(1-loc(1)))*(mxx-mnx)
        warning('off','backtrace')
        warning(['Scale bar of ' num2str(barWidth) ' ' units(1:end-2) ...
            ' is too large, introduce a smaller one in pop-up window'])
        warning('on','backtrace')
        newBar = inputdlg(['Scale bar size smaller than ' num2str(barWidth) ...
            ' [' units(1:end-2) ']'],'New size',[1 45]);
        if ~isempty(newBar) && ~strcmp(newBar,'')
            barWidth = str2double( newBar );
            w = barWidth/pix2nm;
        end
    end
end

% extract verticies
v = VorDat{1};
% extract connectivities
c = VorDat{2};

% make a color mapping based on the central percentiles
vorAreas = extractVorAreas(x);
[arsort,id] = sortrows( vorAreas );

%%

% convert connectivity cell-matrix to array of faces
nVert = cellfun(@size,c,'UniformOutput',false); 
nVert = cell2mat(nVert);
F = nan(size(c,1),max(nVert(:,2)));
Cf = zeros(size(F,1),1);
txtlim = 5000;
if nloc > txtlim*2
    strCR = '';
    txtOut = true;
else
    txtOut = false;
end
strFormat = 'Preparing Voronoi polygons: %3.1f%%%% \n';
for i = 1:size(F,1)
    if txtOut && rem(i,txtlim)==0
        strOut = sprintf(strFormat,100*i/nloc);
        fprintf([strCR strOut])
        strCR = repmat('\b',1,length(strOut)-1);
    end
    F(i,1:nVert(i,2)) = c{i};
end

% threshold colors
Cf = vorAreas;
Cf(Cf<log10(thresholds(1))) = log10(thresholds(1));
valids = Cf<log10(thresholds(2));
% get a matrix of jet colors to match the 
jetCol = jet(range);
jetCol = flipud(jetCol);
fade = 8;
for ii = 1:fade % make the last 4 entries dimmer versions to fade to black
    jetCol(end+1,:) = jetCol(end,:)/2;
end
% match valids to range values
vAreas = Cf(valids);
vRange = max(vAreas)-min(vAreas);
vStep = vRange/(range+fade-1);
svAreas = vAreas-min(vAreas);
svRankColors = fix(svAreas/vStep)+1;
ColorOut = zeros(size(F,1),3);
ColorOut(valids,:) = jetCol(svRankColors,:);

%Cf(Cf>log10(thresholds(2))) = log10(thresholds(2));

%%% Begin plotting

% make figure and set size according to screen size
scrn = get(0,'Screensize'); 
sz = scrn(4)*0.7;
fig_h.fig = figure('position',[scrn(3)*0.05 scrn(4)*0.65*0.3 [1.15 1]*sz]);
fig_h.ax = subplot(1,1,1,'Parent',fig_h.fig);

% plot polygons as patches
hp = patch(fig_h.ax,'Faces',F,'Vertices',v,...
    'FaceVertexCData',ColorOut,'FaceColor','flat');
hp.EdgeColor = 'none';

if txtOut 
    strOut = sprintf(strFormat,100*i/nloc);
    fprintf([strCR strOut])
	disp('Rendering image, please be patient')
end

% adjust axis to make an image
axis(fig_h.ax,'ij')
axis(fig_h.ax,'image')
axis(fig_h.ax,[mnx mxx mny mxy])
set(fig_h.ax,'Color','k'...    
    ,'Yticklabel',[]...
    ,'Xticklabel',[]) 

% include a colorbar
fig_h.cb = colorbar(fig_h.ax);
%t=get(fig_h.cb,'Limits');
colorScale = linspace(0,1,5);
ActualScale = logspace(log10(thresholds(1)),log10(thresholds(2)),5);
set(fig_h.cb,'Ticks',colorScale);
set(fig_h.cb,'TickLabels',ActualScale);
% setting colormap
colormap(jetCol);

% Set the ticks of the colorbar to cover the 0-1 range
% unfortunately, it does not show tick marks at 0 or 1 exactly
%set(fig_h.cb,'Ticks',linspace(0.02,0.99,10))
title(fig_h.cb,'Area/Polygon')
% make 10 increments along the color bar
%ytic = 1-linspace(pctls(1)/100,pctls(2)/100,10);
%ytic = flipud( ytic' )';
% set(fig_h.cb,'Ticks',ytic) % issue: colorbar always goes 0 to 1
% find the corresponding values given the spacing between the percentiles
%ticVal = nan(1,length(ytic));
%for i = 1:length(ytic)/2
%    ticVal([length(ytic)-i+1, i]) = pix2nm^2*10.^...
%        interpercentilerange( vorAreas( ~isnan(vorAreas) ), 1-ytic([length(ytic)-i+1, i]) );
%end
% ticVal
%ticStr = get(fig_h.cb,'TickLabels');
%for i = 1:length(ticVal)
%    ticStr{i} = [num2str(ticVal(i),areaFmt) ' ' units];
%end
% ticStr
%set(fig_h.cb,'TickLabels',ticStr)

if includeScaleBar
    % place bottom-right corner at loc(1)% width & loc(2)% height
    loc = [10 7]/100;
    tr = [(mxx-mnx)*loc(1) loc(2)*(mxy-mny)];
    % thickness = 1.6% y-height
    h = tr(2)*0.016/loc(2);
    bl = [mxx-(tr(1)+w) mxy-(tr(2)+h)];
    fig_h.scale = rectangle(fig_h.ax,'Position',[bl w h],...
        'EdgeColor','none',...
        'FaceColor','w');
    % report the scale bar length in the title
    title(fig_h.ax,['Scalebar = ' num2str( barWidth ) ' ' units(1:end-2)])
    
end % adding scalebar

%% Subfunction
    % function to extract the voronoi areas and transform them to
    % log10-values
    function vorAreas = extractVorAreas(x)
        vorAreas =  x(:,3);
        vorAreas = log10( vorAreas );
    end % of subfunction

    function pctls = range2percentile( rng, x )
        vorAr = extractVorAreas(x);
        % get CDF of Voronoi areas
        [f,ar] = ecdf( vorAr );
        % interpolate CDF values to obtain the corresponding percentile
        % values
        pctls = interp1(ar(2:end),100*f(2:end),rng,'pchip');
    end

end % of function