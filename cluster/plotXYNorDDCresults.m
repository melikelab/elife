function [medianVals,h_figure,h_subplots,h_legend,legendStrings] = ...
                     plotXYNorDDCresults(files,groups,type,extractIterStr, iterStrRg,...
                                    h_figure,h_subplots,h_legend,legendStrings)
                                
% J.Otterstrom Matlab 2013b

% function to plot the CDF, PDF and boxplots of median values for cluster
% results saved in .xyn files output by FindClusters.m

% function call:
% [medianVals,h_figure,h_subplots,h_legend,legendStrings] = ...
%    plotXYNresults(files,groups,type,h_figure,h_subplots,h_legend,legendStrings)

% INPUTS
%  *required:
% files - structure variable containing fields indicated by the input
%       variable 'groups' where these fields are each nx2 cells where the
%       first cell column contains the .xyn filename and the second cell
%       column contains the folder location
%       > ideally this input is the output of the functions
%           Select2DataGroups.m
%           Select1DataGroup.m (run multiple times)
% groups - nx1 or 1xn cell containing the names of the fields in the 
%       'files' structure that indicate the filename & folder location of 
%       the .xyn files containing the cluster data to be plotted
% type - nx1 or 1xn cell containing strings for labeling each of the groups
%       in the legend
% 
%  **optional
% extractIterStr - String pattern for extracting iteration information from 
%   the filename to include in the figure legend. Used in conjunction with
%   input 'iterStrRg'. Default Value = '_list';
% iterStrRg - index value range for string extraction starting 
%   with the first character present in "extractIterStr". Default = [-3 -1];
% h_figure = handle of the figure where the plots should be placed
% h_subplots = 3x3 cell of axes handles for the 9 subplots where data
%       should be plotted
% h_legend = handle array for the items to be plotted
% legendStrings = cell array with strings for each item to be plotted, if
%       not included it is generated from the 'type' input

% OUTPUTS
% medianVals = structure variable with the median values of the .xyn files
%       contained in fields indicated by the 'groups' input
% h_figure = handle to the plotted figure
% h_subplots = 3x3 cell of axes handles corresponding to each of the 9 
%       subplots
% h_legend = handle array of the legend included in the first subplot
% legendStrings = cell array of strings included in the legend

% plot results of cluster analysis
% 3x3 subplots
% top row is for Number of localizations/cluster
% middle row is for nearest in-island neighbor distance
% bottom row is for area per cluster
% Left column is the CDF of all values obtained
% middle column is the PDF of all values
% Right column is a bar-chart of the median values

% initialize general plotting parameters
dataTypes = {'numLoc','# Localizations/cluster','Number of localizations per cluster';...
             'nnd',   'nnd (nm)',              {'Nearest ','neighbor distance'};...
             'area',  'Area (nm^2)',            'Area per cluster'};
axisFont = 8;
titleFont = 10;
lineThickness = 2;
makerType = 's';
markerSize = 10;

posnY = [0.72 0.395 0.075];
barplotwidth = 0.2;
xax = [0 50; ...
       0 100; ...
       0 3000];
upperXpdf = [1 50; 10 500; 10 2e4];
nrows = 3;
ncols = 3;

% order for subplot combining
sporder = [1 4 7 2 5 8 3 6 9];
% sporder{1} = {[1 4]; [7 10]; [13 16]};
% sporder{2} = [2 8 14 5 11 17];
% sporder{3} = {[3 6]; [9 12]; [15 18]};


% Check inputs
% check if figure is already created
if ~exist('fg','var') || ~ishandle(h_figure) 
    h_figure = figure('Position',[240 100 1250 770]);
end
if ~exist('groups','var') || isempty(groups)
    groups = {'data1';'data2'};
end
if ~exist('h_subplots','var') || isempty(h_subplots)
    h_subplots = cell(size(dataTypes,1),ncols);
end

numfiles = 0;
for g = 1:size(groups,1)
    numfiles = max([numfiles, size(files{g}.(groups{g}),1)]);
end

% this is the character string to be used for iteration info
if ~exist('extractIterStr','var') || isempty(extractIterStr)
    extractIterStr = '_list';
end

% these are the character indicies for extracting experiment-iteration info
if ~exist('iterStrRg','var') || isempty(iterStrRg)
    iterStrRg = [-3 -1];
elseif max(size(iterStrRg)) ~= 2 && min(size(iterStrRg)) ~= 1
    error('Improper values included as "iterStrRg" input')
end

if ~exist('legendStrings','var') || isempty(legendStrings)
    legendStrings = cell(size(groups,1)*numfiles,1);
else
    legendStrings = [legendStrings; cell(size(groups,1)*numfiles,1)];
end
if ~exist('h_legend','var') || isempty(h_legend)
    h_legend = nan(size(groups,1)*numfiles,1);
else
    h_legend = [h_legend; nan(size(groups,1)*numfiles,1)];
end

if length(groups) == 2
    colrs = [0 0 0.65;0.95 0 0];
else
    colrs = lines(length(groups));
end

% extract the file type from extension to describe NND type
[~,~,extension]=fileparts(files{1}.(groups{1}){1});
switch extension
    case '.xyn'
        nndType = 'in-island ';
    case '.ddc'
        nndType = 'global ';
    otherwise
        nndType = '';
end
% now reset the title of the NND plots
dataTypes{2,3} = [dataTypes{2,3}{1} nndType dataTypes{2,3}{2}];

% initialize vars to save the median values
  % all arrays must have same size for boxplot to work properly
for g = 1:size(groups,1)
    for p = 1:size(dataTypes,1)
        medianVals.(groups{g}).(dataTypes{p}) = nan(numfiles,1);
    end
end


%make the subplots for the CDF
for p = 1:size(dataTypes,1)
    h_subplots{p,1} = subplot(nrows,ncols,sporder(p),'Parent',h_figure);
%     sp{p,1} = subplot(nrows,ncols,sporder{1}{p},'Parent',fg);
    hold(h_subplots{p,1},'on')
    box(h_subplots{p,1},'on')
    axis(h_subplots{p,1},[xax(p,:) 0 1])
    xlabel(h_subplots{p,1},dataTypes{p,2})
    ylabel(h_subplots{p,1},'Cumulative Fraction')
    title(h_subplots{p,1},dataTypes{p,3},'FontSize',titleFont)
    curpos = get(h_subplots{p,1},'Position');
    set(h_subplots{p,1},'FontSize',axisFont)%,'Position',[0.08 posnY(p) 0.26 curpos(4)])
end

%make the subplots for the PDF && box plots
for g = 1:size(groups,1)
    for p = 1:size(dataTypes,1)
        % PDF
        if p == 1
            h_subplots{p,2}(g) = subplot(nrows,ncols,sporder(p+size(dataTypes,1)),...
                'Parent',h_figure,'XScale','linear','XMinorTick','on');
%             sp{p,2}(g) = subplot(nrows,ncols,sporder{2}(p+(g-1)*size(dataTypes,1)),...
%                 'Parent',fg,'XScale','linear','XMinorTick','on');
            axis(h_subplots{p,2}(g),[upperXpdf(p,:) 0 Inf])
        else
            h_subplots{p,2}(g) = subplot(nrows,ncols,sporder(p+size(dataTypes,1)),...
                'Parent',h_figure,'XScale','log','XMinorTick','on');
%             sp{p,2}(g) = subplot(nrows,ncols,sporder{2}(p+(g-1)*size(dataTypes,1)),...
%                 'Parent',fg,'XScale','log','XMinorTick','on');
            axis(h_subplots{p,2}(g),[upperXpdf(p,:) 0 Inf])
            set(h_subplots{p,2}(g),'XScale','log','XTickLabel',{'10' '100' '1000' '10000'})
        end
        set(h_subplots{p,2}(g),'FontSize',axisFont)
        hold(h_subplots{p,2}(g),'on')
        box(h_subplots{p,2}(g),'on')
        xlabel(h_subplots{p,2}(g),dataTypes{p,2})
%         set(sp{p,2}(g),'XScale','log')
        ylabel(h_subplots{p,2}(g),'Frequency')
        title(h_subplots{p,2}(g),['Aggregated ' dataTypes{p,3}],'FontSize',titleFont)
        
        % box plot
%         sp{p,3} = subplot(nrows,ncols,sporder{3}{p},'Parent',fg);
        h_subplots{p,3} = subplot(nrows,ncols,sporder(p+2*size(dataTypes,1)),'Parent',h_figure);
        hold(h_subplots{p,3},'on')
        box(h_subplots{p,3},'on')
        xlabel(h_subplots{p,3},dataTypes{p,2})
        title(h_subplots{p,3},['Median ' dataTypes{p,3}])
        
    end
end


% put data in prepared subplots
m = 0;
aggData = cell(size(groups,1),1);
% typeFields = fieldnames(type);
maxXval = zeros(3,1);
for g = 1:size(groups,1)
    
%     filegroup = files.(groups{g}); % for use with 'Select2DataGroups.m'
    filegroup = files{g}.(groups{g}); % for use with 'Select1DataGroup.m'
    % file-specific initialize plotting parameters
    linestyl = repmat({'-','--','-.',':'},1,ceil(size(filegroup,1)/4));
    
    for f = 1:size(filegroup,1)
        
        plotData = extractPlotData(fullfile(filegroup{f,2},filegroup{f,1}));
        
        % plot data in each subplot
        m = m+1;
        for p = 1:size(dataTypes,1)
            currData = plotData.(dataTypes{p,1});
            aggData{g}.(dataTypes{p,1}){f,1} = currData;
            % CDF
            [y,x] = ecdf( currData );
            H = plot(h_subplots{p,1},x,y,...
                'Color',colrs(g,:),...
                'LineStyle',linestyl{f},...
                'LineWidth',lineThickness);
            if p == 1, h_legend(m) = H; end
            % PDF, plot using stairs
            if f == size(filegroup,1)
                bindata = cell2mat( aggData{g}.(dataTypes{p,1}) );
                binstart = log10(floor(min(bindata)));
                binend = log10(ceil(max(bindata)));
                good = 0;
                while good == 0
                    bins = logspace( ...
                        binstart, ...
                        binend, ...
                        round(sqrt(length(bindata))))';
                    
                    ct = histc(bindata,bins);
                    logsep = diff( log10(bins) );
                    logsep = logsep(end);
                    if ct(end) == 0 && ct(1) == 0
                        good = 1;
                    elseif ct(end) ~= 0
                        binend = binend + logsep;
                    elseif ct(1) ~= 0
                        binstart = binstart - logsep;
                    end
                end
                space = diff(bins);
                idx = find(ct~=0);
                ct = ct(idx)./(space(idx).*sum(ct)); % normalize
                bins = bins(idx);
                stairs(h_subplots{p,2}(g),[bins(1)*0.999; bins],[0; ct],...
                    'Color',colrs(g,:),...
                    'LineStyle','-');%linestyl{f});
            end
            % a marker for the median value
            val = median( currData );
            height = y( x==val );
            if isempty(height), height = 0.5; end
            medianVals.(groups{g}).(dataTypes{p})(f) = val;
            maxXval(p) = max([maxXval(p),val]);
            plot(h_subplots{p,1}, val, height,...
                'Marker',makerType,...
                'MarkerSize',markerSize,...
                'Color',colrs(g,:),...
                'LineStyle','none',...
                'LineWidth',lineThickness+1)
            
            % include a legend entry
            numbridx = strfind(filegroup{f,1},extractIterStr);
            legendStrings{m} = [type.(groups{g}) ' - ' filegroup{f,1}(numbridx+iterStrRg(1):numbridx+iterStrRg(2))];
%             legendStrings{m} = [type.(groups{g}) ' - ' filegroup{f,1}(numbridx+1:numbridx+3)];
%             legendStrings{m} = [type{g} ' - ' filegroup{f,1}(numbridx+1:numbridx+3)]; % for use with 'Select2DataGroups.m'
        end

    end
    
end
% update x-axis 
for p = 1:size(dataTypes,1)
    set(h_subplots{p,1},'XLim',[xax(p,1) 2.5*maxXval(p)])
end

% add legend
legendStrings = legendStrings(find(~cellfun(@isempty,legendStrings)));
h_legend = h_legend( ~isnan(h_legend) );
legend(h_legend,legendStrings,'Location','SouthEast')

%% make plotting variables for boxplots
boxdata = cell(size(dataTypes,1),1);
boxOrder = size(groups,1):-1:1;
for g = boxOrder
    for p = 1:size(dataTypes,1)
        boxdata{p} = [boxdata{p}, medianVals.(groups{g}).(dataTypes{p})];
    end
end

boxDataNames = cell(1,length(groups));
for i = 1:length(groups)
    boxDataNames{i} = type.(groups{boxOrder(i)});
end
%% add box plots of all median value data

for p = 1:size(dataTypes,1)
    boxplot(h_subplots{p,3}, boxdata{p}, boxDataNames,...
            'orientation','horizontal')
%     boxplot(h_subplots{p,3}, boxdata{p}, {type{length(type):-1:1}},...
%             'orientation','horizontal') % for use with 'Select2DataGroups.m'
    % include the extracted median values above the bar
    for g = size(groups,1):-1:1
        plot(h_subplots{p,3},...
            medianVals.(groups{g}).(dataTypes{p}),...
            0.3+(size(groups,1)-g+1)*ones(size(medianVals.(groups{g}).(dataTypes{p}),1),1),...
            'Color',colrs(g,:),...
            'Marker','*',...
            'LineStyle','none');
    end
    % adjust size & position
    curpos = get(h_subplots{p,3},'Position');
    set(h_subplots{p,3},'YLim',[0 length(groups)+0.7],...
        'FontSize',axisFont,...'Position',[curpos(1) posnY(p) barplotwidth(1) curpos(4)])
        'Position',[curpos(1) curpos(2) barplotwidth(1) curpos(4)])
end

% subfunctions

    function plotData = extractPlotData(filename)
        
        %% determine the type of file selected
        [~,~,ext] = fileparts( filename );
        
        switch ext
            case '.xyn'
                % Read the XYN file
                xynInfo = XYN( filename );
                
                % extract the desired quantities from the xynInfo
                nndXY = xynInfo.data(:,10); % units = nm
                nndXY = nndXY( ~isinf(nndXY) );
                clusterArea = (xynInfo.params.original_pixel_size*xynInfo.data(:,6)).^2.*pi; % units = nm^2
                numberLocs = xynInfo.data(:,3);
                
            case '.ddc'
                ddcInfo = DDC( filename );
                nndXY = ddcInfo.data{1}(2:end,1); % units = nm
                nndXY = nndXY( ~isnan(nndXY) );
%                 nndXYZ = ddcInfo.data{1}(:,2);
%                 nndXYZ = nndXYZ( ~isnan(nndXYZ) );
                clusterArea = ddcInfo.data{1}(:,3); % units = nm^2
                clusterArea = clusterArea( ~isnan(clusterArea) );
                numberLocs = ddcInfo.data{1}(:,4);
                numberLocs = numberLocs( ~isnan(numberLocs) );
                
            otherwise
                error(['The selected "' ext '" filetype is not supported'])
        end
        
        % keep the info in a structure for ease of plotting
        plotData.nnd = unique(nndXY); % remove double counting for pairs of clusters, which give repeats of a single nnd value
        plotData.area = clusterArea;
        plotData.numLoc = numberLocs;
    end % of sub function
end