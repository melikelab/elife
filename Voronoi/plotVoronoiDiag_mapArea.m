% function that receives the output of calling VoronoiAreas.m and will plot
% the Vornoi Tessellation diagram color-coding each cell according to the
% area it occupies.
% The area can be converted from pixels to nm by an optional input
%
% INPUTS
%   x - 3 or 4 column list in the format 
%       [x_coord, y_coord, area_zeroRank, area_firstRank]
%       that is output by the VoronoiAreas.m function
%   VorDat - 2-element cell array.  The first cell contains the vornoi
%       verticies, 'V', and the second element contains the connectivity
%       matrix for plotting the voronoi polygons.
%       These values will be calculated if an empty vector is passed in.
%   pix2nm - (optional) conversion from pixel units to nanometers.
%       Values < 1 are considered to be in units of mirometers, and so the
%       conversion is updated as: pix2nm = pix2nm*1000
%       to make the units nm again
%   barWidth - (optional) adds a scale bar to the image of Voronoi polygons
%       if a value of pix2nm is entered, then the scale bar will be in 
%       units of nanometers, otherwise the scale bar will be in units of
%       pixels
% Name-value pair inputs
%   'percentile', pctls - (optional) requires either 'percentile' or 'pctl'
%       string combined with a 1x2 row-vector containing the percentile 
%       range of the data to be used for color-coding the Voronoi polygons. 
%       Default is [1 95] such that the 1st percentile of areas (smallest) 
%       are color-coded yellow, while the 95th percentile of areas 
%       (largest) are color-coded blue.  
%       NOTE:
%       - This input is permitted only if the area range is not also 
%           stated.
%       - The program requests that the percentile range be larger than 5% 
%           of the total data, which can be bypassed.
%   'range', range - (optional) requires either 'range' or 'rng' string
%       combined with a 1x2 row-vector containing the lower and upper
%       bounds for Voronoi polygon are color-coding, small being yellow and
%       large being blue.
%       NOTE:
%       - The range values input must correspond to the pix2nm conversion, 
%           if it is input, since it is rescaled according to: 
%           range/(pix2nm^2).
%       - This input is permitted only if the percentile range is not also 
%           stated.
%       - The program requests that the percentile range corresponding to 
%           the specified area range consist of at least 5% of the the 
%           total data, which can be bypassed.
%   'Override', value - (optional) requires either 'Over' or 'over' or
%       'Override' or 'override' followed by a single logical 
%       true/false/1/0 value that permits bypassing the request that the 
%       percentile range for color-coding span 5% of the data.
% 
% OUTPUTS
%   fig_h - structure containing three fields:
%       .fig = figure handle
%       .ax = axis hangle
%       .cb = color bar handle
%       .scale = scale bar handle (if requested); see the rectangle.m
%           function for a list of properties
%
% 
% Function call examples:
%   fig_h = plotVoronoiDiag_mapArea(x);
%   fig_h = plotVoronoiDiag_mapArea(x,[]);
%       the 'VorDat' variable containing voronoi information is calculated
%       internal to the program, but is not returned.
%       area units will be in pixel^2
%
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat);
%       area units will be in pixel^2
% 
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, pix2nm);
%       area units will be in nm^2
% 
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, pix2nm, barWidth);
%       area units will be in nm^2, scale bar will be in units of nm
% 
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, [], barWidth); 
%       area units will be in pixel^2, scale bar will be in units of pixels
% 
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, [], [], 'percentile', pctls);
%       area units will be in pixel^2, No scale bar 
%       area color-coding according to input percentiles
% 
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, [], [], 'percentile', 50);
%       area units will be in pixel^2, No scale bar 
%       makes a binary-type image with area < median as yellow, rest as blue
%       User is querried regarding the narrow data range colored
% 
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, [], [], 'percentile',[0 1], 'Override',true);
%       area units will be in pixel^2, No scale bar 
%       plots the 1% smallest areas as yellow, no scale bar, area in pixel^2
%       User is not querried regardign the narrow data range colored
%
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, pix2nm, [],'range',[20 1000]);
%       area units will be in nm^2, No scale bar
%       Color-codes all polygons having areas between 20-1000 nm^2
%       Area in nm^2 and user can be quierried regarding narrow ranges
%
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, pix2nm, [],'range',[20 50],'Override',true);
%       area units will be in nm^2, No scale bar
%       Color-codes all polygons having areas between 20-50 nm^2
%       Area in nm^2 and user will not be quierried regarding narrow ranges
%
%   fig_h = plotVoronoiDiag_mapArea(x, VorDat, [], [],'range',[0.002 0.02]);
%       area units will be in pixel^2, No scale bar
%       Color-codes all polygons having areas between 0.002-0.02 pixel^2
%       Area in pixel^2 and user can be quierried regarding narrow ranges
%
%
%   Created: Oct 31, 2017; J.Otterstrom, Matlab 2016a
%   Updates: 
%       Nov 15, 2017, J.O., added rectangle handle to output

function fig_h = plotVoronoiDiag_mapArea(x, VorDat, pix2nm, barWidth, varargin)

% minimum percentage of the data that must be covered by the percentile
% range
minDataCover = 5;

% ensure the input data has the needed Voronoi information
if size(x,2) < 3 || nargin < 2 || ~iscell(VorDat)
    [x, ~,VorDat] = VoronoiAreas(x(:,1:2));
end


% set the pixel conversion, if input
if nargin < 3 || isempty(pix2nm)
    pix2nm = 116;
    units = 'pix^2';%[];
    areaFmt = '%4.4f';
else
    if pix2nm < 1
        disp(['converting input pixel size of ' num2str(pix2nm) ' um to ' num2str(pix2nm) ' nm'])
        pix2nm = pix2nm*1000;
    end
        units = 'nm^2';
        areaFmt = '%4.0f';
end

% check to see if there should be a scalebar
if nargin < 4 || isempty(barWidth)
    includeScaleBar = false;
else
    includeScaleBar = true;    
end

% if nargin < 5 || isempty(pctls)
%     pctls = [1 95];
% else
%     % ensure correct order & range
%     pctls = check_pctls( pctls );
% end
% 
% if nargin < 6 
%     pctlOverride = false;
% elseif ~islogical( pctlOverride )
%     pctlOverride = logical( pctlOverride(1) );
% end

% initialize all variable inputs 
inputRange = [];
pctls = [];
pctlOverride = false;

% process variable input 
if nargin >= 5
    % check that there is an even number of inputs
    if rem(nargin-4,2)~=0
        error('Variable input arguments must come in name-value pairs')
    end
    
    % process the inputs
    for i=1:2:length(varargin)
        % check, again, that inputs are properly formated name-value pairs
        if ~ischar( varargin{i} ) || ~ismatrix( varargin{i+1} )
            error('Variable input arguments must come in name-value pairs')
        end
        % extract current name & value
        currName = varargin{i};
        currValue = varargin{i+1};
        if strcmpi(currName,'percentile') || strcmpi(currName,'pctl')
            currName = 'percentile';
        end
        if strcmpi(currName,'range') || strcmpi(currName,'rng')
            currName = 'range';
        end
        if ~isempty( strfind(currName,'over') ) || ~isempty( strfind(currName,'Over') )
            currName = 'Override';
        end
        switch currName
            case 'percentile'
                pctls = checkVectIn( currValue, 'percentile' );
                % ensure correct order & range
                pctls = check_pctls( pctls );
            case 'range'
                inputRange = checkVectIn( currValue, 'range' );
                % ensure correct order & range
                inputRange = check_range( inputRange, x, pix2nm );
                               
            case 'Override'
                pctlOverride = logical( currValue );
        end
    end
end
   
% catch to see if the percentile values were input
if isempty( pctls ) && isempty( inputRange )
    % neither was input by user
    pctls = [2 95];
end
    

%% Begin data preparation

% check to see if the user has input both percentiles and inputRange
if ~isempty(inputRange) && ~isempty(pctls)
    % both are input, so user must choose one or the other
    useInput = [];
    choices = {' Percentile ',' Areas '};
    dispRange = pix2nm^2*10.^inputRange;
    while isempty(useInput)
        useInput = questdlg({['You have input both an area range and a percentile '...
            'range for color-coding.  Entered values are:'];...
            ['Percentile range: ' num2str(pctls(1)) ', ' num2str(pctls(2))];...
            ['Area value range: ' num2str(dispRange(1),areaFmt) ', ' num2str(dispRange(2),areaFmt)];...
            'Choose one to use.'} ...
            ,'Choose input to use' ...
            ,choices{1},choices{2}, choices{1});
    end
    switch useInput
        case choices{1}
            inputRange = [];
        case choices{2}
            pctls = [];
    end
end

% if the input is an inputRange, then extract the percentiles it
% corresponds to
if ~isempty(inputRange) 
    pctls = range2percentile( inputRange, x );
end

% ensure the percentile range covers more than 5% of the input data, unless
% the user has indicated they wish to bypass this check
if ~pctlOverride 
    while pctls(2) - pctls(1) < minDataCover
        lowCov = [];
        % ask if user wants to continue with this small data-range
        while isempty(lowCov)
            lowCov = questdlg({['Current percentile range from ' num2str(pctls(1)) ' to ' ...
                num2str(pctls(2)) ' constitutes less than 5% of the data. '...
                'Do you wish to continue with these values?']}...
                ,'Low data coverage'...
                ,'Yes','No','No');
        end
        
        switch lowCov
            case 'Yes'
                % continue with current values
                break
            case 'No'
                % ask user to input new values
                newpctls = [];
                while isempty(newpctls)
                    newpctls = inputdlg(...
                        {'Lower Voroni area percentile for yellow colors'; ...
                        'Upper Voronoi area percentile for blue colors'}...
                        ,'New Percentile Range'...
                        ,[1 50; 1 50]...
                        ,{num2str(pctls(1));num2str(pctls(2))} );
                end
                pctls = str2double(newpctls);
                % ensure correct order & range
                pctls = check_pctls( pctls );
        end
    end
end
npctls = length(unique(pctls));

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
% % find the minimum verticie in x
% [~,mnxi] = min( x( id(1:size(colr,1)),1 ) );
% mnx = min( v( c{ id(mnxi) }, 1) );
% 
% % find the maximum verticie in x
% [~,mxxi] = max( x( id(1:size(colr,1)),1 ) );
% mxx = max( v( c{ id(mxxi) }, 1) );
% 
% % find the minimum verticie in y
% [~,mnxi] = min( x( id(1:size(colr,1)),2 ) );
% mny = min( v( c{ id(mnxi) }, 2) );
% 
% % find the maximum verticie in y
% [~,mxxi] = max( x( id(1:size(colr,1)),2 ) );
% mxy = max( v( c{ id(mxxi) }, 2) );

% extract verticies
v = VorDat{1};
% extract connectivities
c = VorDat{2};

% make a color mapping based on the central percentiles
vorAreas = extractVorAreas(x);
[arsort,id] = sortrows( vorAreas );
pctls
%%
if npctls > 1
    % get the percentile range for color-coding
    rng = interpercentilerange( vorAreas( ~isnan(vorAreas) ),pctls/100);
    endidx = find(isnan(arsort),1,'first');
    idx = [find(arsort>rng(1),1,'first'),...
           find(arsort<rng(2),1,'last'),...
           endidx];
    if isempty( endidx ), idx = [idx, idx(2)]; end
    colrSort = flipud( parula(idx(2)-idx(1)) );
    colrSort(end,:) = [0,0,0];
    colrStart = colrSort(1,:);
    colrEnd = colrSort(end,:);
else
    % get the percentile value
    if pctls(1)+1 > 100
        rng = interpercentilerange( vorAreas( ~isnan(vorAreas) ),[pctls(1)-1 pctls(1)]/100);
        rng = rng(2);
    else
        rng = interpercentilerange( vorAreas( ~isnan(vorAreas) ),[pctls(1) pctls(1)+1]/100);
        rng = rng(1);
    end
    idx = find(arsort>rng,1,'first');
    idx = [1 1]*idx;
    idx = [idx, find(isnan(arsort),1,'first')];
    colrSort = flipud( parula(2) );
    colrSort(end,:) = [0,0,0];
    colrStart = colrSort(1,:);
    colrEnd = colrSort(end,:);
    colrSort = [];
end
% %%
% colrStart = colrSort(1,:);
% colrEnd = colrSort(end,:);
colrSort = [repmat(colrStart,idx(1)-1,1); ... 0-2.5 pctl are same as first color
        colrSort; ... main data 
        repmat(colrEnd,idx(3)-idx(2),1); ... last 2.5 pctl same as last color
        repmat([0 0 0],length(arsort)-idx(3)+1,1)]; % NaN values are black
colr = nan(length(arsort),3);
%%
for i = 1:length(arsort)
    colr(id(i),:) = colrSort(i,:);
end

% convert connectivity cell-matrix to array of faces
nVert = cellfun(@size,c,'UniformOutput',false); 
nVert = cell2mat(nVert);
F = nan(size(c,1),max(nVert(:,2)));
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


%%% Begin plotting
% make figure and set size according to screen size
scrn = get(0,'Screensize'); 
sz = scrn(4)*0.7;
fig_h.fig = figure('position',[scrn(3)*0.05 scrn(4)*0.65*0.3 [1.15 1]*sz]);
fig_h.ax = subplot(1,1,1,'Parent',fig_h.fig);


% plot polygons as patches
hp = patch(fig_h.ax,'Faces',F,'Vertices',v,...
    'FaceVertexCData',colr,'FaceColor','flat');
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
% axis(fig_h.ax,[min(v(2:end,1)) max(v(2:end,1)) min(v(2:end,2)) max(v(2:end,2))])
set(fig_h.ax,'Color','k'...    
    ,'Yticklabel',[]...
    ,'Xticklabel',[]) 

if npctls > 1
    % include a colorbar
    fig_h.cb = colorbar(fig_h.ax);
    % Set the ticks of the colorbar to cover the 0-1 range
        % unfortunately, it does not show tick marks at 0 or 1 exactly
    set(fig_h.cb,'Ticks',linspace(0.02,0.99,10))
    title(fig_h.cb,'Area/Polygon')
    % make 10 increments along the color bar
    ytic = 1-linspace(pctls(1)/100,pctls(2)/100,10);
    ytic = flipud( ytic' )';
    % set(fig_h.cb,'Ticks',ytic) % issue: colorbar always goes 0 to 1
    % find the corresponding values given the spacing between the percentiles
    ticVal = nan(1,length(ytic));
    for i = 1:length(ytic)/2
        ticVal([length(ytic)-i+1, i]) = pix2nm^2*10.^...
            interpercentilerange( vorAreas( ~isnan(vorAreas) ), 1-ytic([length(ytic)-i+1, i]) );
    end
    % ticVal
    ticStr = get(fig_h.cb,'TickLabels');
    for i = 1:length(ticVal)
        ticStr{i} = [num2str(ticVal(i),areaFmt) ' ' units];
    end
    % ticStr
    set(fig_h.cb,'TickLabels',ticStr)
end % add colorbar


if includeScaleBar
    % place bottom-right corner at loc(1)% width & loc(2)% height
    loc = [10 7]/100;
    tr = [(mxx-mnx)*loc(1) loc(2)*(mxy-mny)];
    % rectangle function requires [x y w h]
        % [x y] of bottom-left corner
        %  w = width, defined above
        %  h = height
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
    % function to check that the percentile range is appropriate
    function pctls = check_pctls( pctls )
        % ensure [smallest biggest]
        pctls = unique(pctls);%[min(pctls) max(pctls)];
        if length(pctls)==1
            pctls = [pctls, pctls];
        else
            pctls = [min(pctls) max(pctls)];
        end
        if pctls(1) < 0, pctls(1) = 0; end
        if pctls(2) > 100, pctls(2) = 100; end
    end % of check_pctls subfunction
    
    % check to see that the area range is acceptable and convert to log10
    function range = check_range( range, x, pix2nm )
        vorAr = extractVorAreas(x);
        range = log10(range/(pix2nm^2));
        % ensure [smallest biggest]
        range = unique(range);
        if length(range)==1
            range = [range, range];
        else
            range = [min(range) max(range)];
        end
        if range(1) < min(vorAr), range(1) = min(vorAr); end
        if range(2) > max(vorAr), range(2) = max(vorAr); end
    end % of check_range subfunction

    function vectOut = checkVectIn( vectIn, name )
        vectOut = vectIn;
        % check if it is a matrix
        if ~isvector(vectOut)
            warning('off','backtrace')
            warning(['Detected a ''' name ''' input in matrix format. ' ...
                ' Only the first row will be utilized.'])
            warning('on','backtrace')
            vectOut = vectOut(1,:);
        % make sure it is a column vector
        elseif iscolumn(vectOut)
            vectOut = vectOut';
        end
    end % of subfunction

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