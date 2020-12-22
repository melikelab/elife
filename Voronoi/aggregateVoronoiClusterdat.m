
% clear workspace 
clear, clc

pix2nm = 160; % nm

startlocn = {...
    'V:\Jason_O\Data\170712 - DNA H2B paint NSTORM\NoTSA\Cell 2'};

% define data types
type.t1 = 'No TSA';
type.t2 = '+TSA';
% type.t1 = 'T47D';
% type.t2 = 'SKBR3';
% type.t3 = 'BT549';
% type.t4 = 'MDA MB231';
% type.t5 = '160129 MeEtOH';
% type.t6 = '160203 PFA';
% type.t7 = 'Full Nuc';
% type.t8 = 'Excl Crtx';
% optional path to start looking for files


% call function allowing user to select files
groups = fieldnames(type); 
ntypes = size(groups,1); 
files = cell(ntypes,1);
if length(startlocn) ~= ntypes
    startlocn = repmat(startlocn(1),ntypes,1);
end

% have user select the .mat files to be aggregated
for t = 1:ntypes
    files{t} = Select1DataGroup(type.(groups{t}),'mat',startlocn{t});
%     files{t} = RenameField(files{t},'data',groups{t});
    nfiles = size(files{t}.data,1);
end

% go through each type and aggregate
for t = 1:ntypes
%     cluster.(groups{t}) = cell(nfiles,1);
    
    n = 0;
    for f = 1:nfiles
        clusterdat = load( fullfile(files{t}.data{f,2},files{t}.data{f,1}), 'cluster', 'Allthresholds');
        threshs = clusterdat.Allthresholds(:,1);
        if length(threshs) > 1 % multiple thresholds, ask which to use
            %%
            useall = questdlg(['Multiple clustering thresholds identified for ' ...
                type.(groups{t}) ', file number ' num2str(f) '.  Would you like to use all of them?'],...
                'use all cluster thresholds?','Yes','No','Yes');
            if isempty(useall) || strcmpi( useall,'No' )
                threshIdx = []; loop = 0; maxloop = 7;
                while isempty(threshIdx) && loop < maxloop
                    loop = loop+1;
                    threshs = round(threshs,5);
                    thstr =  num2str(threshs(1));
                    for th = 2:length(threshs), thstr = [thstr ', ' num2str(threshs(th))]; end
                    input = inputdlg({'Which of the detected thresholds would you like to use (choose only one)?'},...
                        'Multiple cluster thresholds',[1 50],{thstr});
                    if ~isempty(input)
                        input = str2double( input{1} );
                        threshIdx = find(input == threshs);
                    end
                end
            else
                threshIdx = 1:length(threshs);
            end
        else
            threshIdx = 1;
        end
        
        for i = threshIdx
            n = n+1;
            if iscell(clusterdat.cluster)
                currcluster = clusterdat.cluster{i};
            else
                currcluster = clusterdat.cluster;
            end
        
            cluster.(groups{t}){n,1} = ...
                [currcluster.nLocs, [currcluster.areas, currcluster.HullArea].*pix2nm.^2];
        end
    end
    cluster.(groups{t}) = cell2mat( cluster.(groups{t}) );
end
    
open cluster
disp(' "cluster" variable format: [number_localizations, area (nm^2), convexHull area (nm^2)]')

