% putVoronoiAreasExcel.m
% script to aggregate voronoi files for Angel and parse them into excel
% spread sheets for each sub-folder
% Author: PKR, UPENN Sepetember 17, 2019

% user needs to set the tree directory here, everything else is automatic
%tree = 'Q:\Angel\CRG\STORM\Heterokaryon Reprogramming\HK-Voronoi_thresh005';
tree = uigetdir(path,'select tree folder for parsing analysis.');

% get directories and sub directories
dirCandidates = regexp(genpath(tree),'[^;]*','match');

for ii = 1:length(dirCandidates)
    filelist = dir(dirCandidates{ii});
    filelist = filelist(~[filelist.isdir]);
    if isempty(filelist)
        continue
    end
    areas = {};
    names = {};
    % remove files that aren't .mat
    valids = zeros(length(filelist),1);
    for jj = 1:length(filelist)
        if strcmp(filelist(jj).name(end-3:end),'.mat')
            valids(jj) = 1;
        end
    end
    filelist(~valids) = [];
    if isempty(filelist)
        continue
    end
    
    for jj = 1:length(filelist)
        names{jj} = filelist(jj).name;
        filename = fullfile(filelist(jj).folder,names{jj});
        vorData = load(filename);
        if iscell(vorData.cluster)
            areas{jj} = vorData.cluster{1}.areas*vorData.pix2nm*vorData.pix2nm;
        else
            areas{jj} = vorData.cluster.areas*vorData.pix2nm*vorData.pix2nm;
        end
    end
    % create cell array to write out excel file
    maxVal = 0;
    for jj = 1:length(areas)
        maxVal = max(maxVal,length(areas{jj}));
    end
    outputArray = cell(maxVal+1,length(names));
    for jj = 1:length(names)
       outputArray{1,jj} = names{jj}; 
    end
    for jj = 1:length(names)
        for kk = 1:length(areas{jj})
            outputArray{1+kk,jj} = areas{jj}(kk);
        end
    end
    xlswrite([dirCandidates{ii} '\aggregatedAreas.xls'], outputArray);
end