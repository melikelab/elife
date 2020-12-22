% Script creating a new Insight3 .bin file that enables visualization of
% localization density topology detected through iterative cluster
% identification using SR-Tesseler program for Voronï tesselation
% A maximum of 9 iterations of cluster identification can be visualized
% with each segmentation assigned to channels 1-9 of Insight3
%
% J.Otterstrom MATLAB 2016a
% depends on use of Insight3 class developed by Joseph Borbely
% depends on prior use also of Select2DataGroups.m
%
% INPUTS
%   LL - Insight3 object
%   SRTfiles - list of files defined using Select2DataGroups.m
%
% OUTPUTS
% the outputs are two Insight3 .bin files having the same filename as the
% input but appended with '_Topology' or '_Base' for the localization lists 
% showing either the topology via channel assignment or the base image,
% respectively.
% They are saved to the same folder as the input localization list

function Vis_SRT_topology(I3file,SRTfiles)

%% if ~isa(LL,'Insight3')
%     error('First input must be an Insight3 object')
% end
if ~isstruct(I3file) || ~isfield(I3file,'data')
    error('First input must be a file list obtained using Select1DataGroup')
end
if ~isstruct(SRTfiles) || ~isfield(SRTfiles,'data1') || ~isfield(SRTfiles,'data2')
    error('Second input must be a file list obtained using Select2DataGroups')
end
% select the Insight3 Localization List
% I3file = Select1DataGroup('Insight3 Reference list','.bin');
LL = Insight3(fullfile(I3file.data{1,2},I3file.data{1,1}));
% copy the original LL datatable
Loc = LL.getData();
idxCh = LL.getColumnIndex('channel');

% select the base object and cluster objects selected by increasing
%   thresholds
% SRTfiles = Select2DataGroups('Base Object','Clusters','.xls');
% read base object
objectData = dlmread( fullfile(SRTfiles.data1{1,2},SRTfiles.data1{1,1}),...
    '\t',1,0);
disp(['Read base object file: ' SRTfiles.data1{1,1}])

% Detect the cluster objects
nClustGrp = size(SRTfiles.data2,1);
if nClustGrp > 9
    error('maximum of 9 cluster groups permitted due to I3 visualization limitations')
else
    disp(['Detected ' num2str(nClustGrp) ' cluster groupings'])
end

% check to see if a previously saved index pairing has been saved
fileList = recursiveFindFile( SRTfiles.data1{1,2}, '*_SRT_i3_trf.mat');
choose = 0;
if ~isempty(fileList)
    nFiles = size(fileList,1);
    disp('Found SR-Tesseler <=> Insight3 pairing files: ')
    for i = 1:nFiles
        [~,f,e] = fileparts(fileList{i});
        disp([ num2str(i) ') ' f e])
    end
    choose = nFiles+1;
    while choose > nFiles
        choose = input('Which file would you like to load? (0 to skip) ==');
        if ischar(choose) || choose > nFiles || choose < 0
            disp(['Input is not useful, please enter a number between 0 and ' num2str(nFiles) ])
        end
    end
end

if isempty(fileList) || choose == 0
    % send I3 localizations & SR-Tesseler object for index pairing
    disp('Beginning localization list pairing')
    SRT2i3_trf = pairSRT_I3(LL, objectData);
    % save the index pairing list
    [p,f,~] = fileparts(fullfile(SRTfiles.data1{1,2},SRTfiles.data1{1,1}));
    save(fullfile(p, [f '_SRT_i3_trf.mat']),'SRT2i3_trf')
else
    SRT2i3_trf = load( fileList{choose},'SRT2i3_trf' );
    SRT2i3_trf = SRT2i3_trf.SRT2i3_trf;
end
% remove localizations assigned to object 0 from Loc data table
% idxKeep = find(SRT2i3_trf(:,2)~=0);
% idxRem = find(SRT2i3_trf(:,2)==0);
% change their channel to 0
% Loc(SRT2i3_trf( idxKeep ,3),idxCh) = 0;
Loc(SRT2i3_trf( SRT2i3_trf(:,2)~=0 ,3),idxCh) = 0;
% change channel of the others for later removal
% Loc(SRT2i3_trf( idxRem ,3),idxCh) = -1;
Loc(SRT2i3_trf( SRT2i3_trf(:,2)==0 ,3),idxCh) = -1;
% Also for I3 localizations not present in the SR-T object
Loc( ~ismember((1:LL.numMolecules)',SRT2i3_trf(:,3)), idxCh) = -2;
% report to user how many localizations have been removed
disp(['Number of I3 localizations not present in SR-Tesseler data: ', ...
    num2str(size(find(Loc(:,idxCh)==-2),1),'%d') ])

% iteratively read the cluster groups 
% increase the channel assignment for those localizations included in the cluster
% channel:color correspondence, intended to be in order of increasing densities
% 0:white, 1:red, 7:orange, 4:yellow, 2:green, 6:cyan, 8:cyan-blue, 3:blue,
% 9:violet, 5:magenta
if nClustGrp==2
    ColorOrder = [1 8];
elseif nClustGrp==3
%     ColorOrder = [1 4 8];
    ColorOrder = [3 4 1];
elseif nClustGrp==4
    ColorOrder = [1 4 8 5];
elseif nClustGrp==5
    ColorOrder = [1 4 2 8 5];
elseif nClustGrp==6
    ColorOrder = [1 4 2 6 3 5];
elseif nClustGrp==7
    ColorOrder = [1 7 4 2 6 3 5];
elseif nClustGrp==8
    ColorOrder = [1 7 4 2 6 3 9 5];
elseif nClustGrp==9
    ColorOrder = [1 7 4 2 6 8 3 9 5];
end

for f = 1:nClustGrp
    clusterData = dlmread( fullfile(SRTfiles.data2{f,2},SRTfiles.data2{f,1}),...
        '\t',1,0);
    disp(['Read cluster object: ' SRTfiles.data2{f,1}])
    Loc( SRT2i3_trf(clusterData(:,2)~=0,3), idxCh ) = ColorOrder(f);
end

% now remove all localizations with channel < 0, since the index order doesn't
% matter any longer
Loc = Loc( Loc(:,idxCh) >= 0, :);

% Save a new LL for I3 visualization
LL.setData( Loc );
LL.forceFileOverwrite(true);
LL.write([LL.filename(1:end-4) '_Topology_' ...
    num2str(nClustGrp) 'level' LL.filename(end-3:end)]);

% make an additional LL for visualizing the base object
    % all localizations in Channel 0
origCh = Loc(:,idxCh);
Loc(:,idxCh) = 0;
LL.setData( Loc );
LL.forceFileOverwrite(true);
LL.write([LL.filename(1:end-4) '_Base' LL.filename(end-3:end)]);