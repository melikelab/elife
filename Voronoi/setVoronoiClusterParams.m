



function ...
    [pixVals,pix2nm,finalpix,clusterID,iter,signif,minNLoc,maxloop,showIM,files] = ... 
        setVoronoiClusterParams(startdir)
%%
if ~exist('startdir','var') || isempty(startdir)
    startdir = pwd;
end

dlgTitle = 'Voronoi clustering parameters';
defaults = ...
    {'116.99999',      'nm per pixel';...                      1
     '20',       'analysis pixel size (nm)';...           2
     'manual','threshold method: auto/manual/cdf/pval';...     3
     '5',        'Num Monte Carlo iterations';...        4
     '99.5',     'MC confidence bounds (%)';...          5
     '7',        'Minimum # Loc/cluster';...             6
     '8',        'Maximum # thresholds                (in addition to uniform threshold)';...              7
     'false',    'Show Tesselation and intermediate plots? (true,1/false,0)';...8
      startdir,        'Folder for data selection'};%    9
W = 37;
nLines = ...
    [repmat([1 W],size(defaults,1)-1,1); ...
     2 2*W];

inputs = inputdlg(defaults(:,2),dlgTitle,nLines,defaults(:,1));

% nanometers per pixel
in = 1;
pix2nm = str2double(inputs{in});

% pixel size of the mask and for for cluster plotting 
in = in+1; 
finalpix = str2double(inputs{in});

pixVals = [pix2nm,finalpix];

% Cluster identification can be manual or automatic
in = in+1; 
clusterID = [];
while isempty(clusterID)
    if  sum(strcmpi(inputs{in},{'a','auto','automatic'}))
        clusterID = 'automatic'; % performs Monte Carlo simulation
    elseif  sum(strcmpi(inputs{in},{'m','man','manual'}))
        clusterID = 'manual';
    elseif sum(strcmpi(inputs{in},{'c','cdf','percentile'}))
        clusterID = 'cdf';
    elseif sum(strcmpi(inputs{in},{'p','pval','pvalue'}))
        clusterID = 'pval';
    else
        r=3;
        tmp = inputdlg(defaults(r,2),dlgTitle,nLines(r,:),defaults(r,1));
        inputs{in} = tmp{1};
    end
end

% number of Monte Carlo iterations
in = in+1; 
iter = str2double(inputs{in}); % for automatic ID only
if iter == 1
    warndlg({'Selecting 1 Monte Carlo simulation prevents multi-level thresholding.';...
        'Data will be thresholded by comparison to a uniform distribution as per S.R. Tesseler.'})
end

% percent confidence bounds for Monte Carlo output
in = in+1; 
signif = str2double(inputs{in}); % for automatic ID only

% minimum number of localizations allowed in each cluster
in = in+1; 
minNLoc = str2double(inputs{in});

% number of loops to perform when finding thresholds
in = in+1;
maxloop = str2double(inputs{in});

% true/false value for plotting the localizations, their Voronoi
% Tesselation and the identified clusters
in = in+1; 
if sum(strcmpi(inputs{in},{'true','1'}))
    showIM = true;%str2double(inputs{7});
else
    showIM = false;
end

% starting place for selecting Localization list
in = in+1; 
if exist(inputs{in},'dir')
    startpath = inputs{in};
else
    startpath = pwd;
end

% now select localization lists to be analyzed
files = Select1DataGroup('Voronoi Clustering Lists','bin',startpath);
    
    
end % of function