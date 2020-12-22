



function ...
    [pix2nm, barWidth, method, range, override, startdir] = ... 
        setVoronoiDiagParams(startdir)
%%

% method = 'percentile' or 'range'
% range
% pix2nm
% barWidth
% override

if ~exist('startdir','var') || isempty(startdir)
    startdir = pwd;
end

dlgTitle = 'Voronoi clustering parameters';
defaults = ...
    {'160',         'nm per pixel';...                      1
     'percentile',  'Voronoi diagram coloring method method. Options = "percentile" or "range"';...          2
     '5',           'Lower percentile or area-range';...     3
     '95',          'Higher percentile or area-range';...        4
     '',            'Scale bar length [nm]';...          5
     'false',       'Over-ride minimum data span for coloring (1/0,true/false)';...
      startdir,        'Folder for data selection'};%    9
W = 40;
nLines = ...
    [repmat([1 W],size(defaults,1)-1,1); ...
     3 2*W];

inputs = inputdlg(defaults(:,2),dlgTitle,nLines,defaults(:,1));

% nanometers per pixel
in = 1;
pix2nm = str2double(inputs{in});

% color-limits calculation method: percentile or fixed range
in = in+1; 
method = [];
while isempty(method)
    if  sum(strcmpi(inputs{in},{'per','pct','percentile'}))
        method = 'percentile'; % performs Monte Carlo simulation
    elseif  sum(strcmpi(inputs{in},{'r','rg','range'}))
        method = 'range';
    else
        r=3;
        tmp = inputdlg(defaults(r,2),dlgTitle,nLines(r,:),defaults(r,1));
        inputs{in} = tmp{1};
    end
end


% Lower and upper values for coloring
in = in+1; 
rng1 = str2double(inputs{in});
in = in+1; 
rng2 = str2double(inputs{in});

range = [min([rng1, rng2]), max([rng1, rng2])];


% scale bar option
in = in+1; 
barWidth = str2double(inputs{in}); 

% true/false value for plotting the localizations, their Voronoi
% Tesselation and the identified clusters
in = in+1; 
if sum(strcmpi(inputs{in},{'true','1'}))
    override = true;%str2double(inputs{7});
else
    override = false;
end

% starting place for selecting Localization list
in = in+1; 
if exist(inputs{in},'dir')
    startdir = inputs{in};
else
    startdir = pwd;
end
    
    
end % of function