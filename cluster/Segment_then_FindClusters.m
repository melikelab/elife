%function Segment_then_FindClusters(findClustersStruct, segment)
%
% Segments the Molecule List into subregions and the FindClusters function
% is called for each subregion.
%
% Inputs
% ------
% findClustersStruct : custom-written struct
% (see the documentation in FindClustersStruct.m)
%
% segment : struct
%  segment.numDivisions : integer
%    segment the Molecule List into numDivisions^2 subregions
%  segment.viewDivisions : boolean
%    display a plot showing how the molecule list is segmented
%
% Returns
% -------
% nothing
%
function Segment_then_FindClusters(findClustersStruct, segment)

numDivisions = segment.numDivisions;
viewDivisions = segment.viewDivisions;

% Chop the i3 molecule list into squares and send each separately for cluster analysis
MLlimits.x = [ (0:numDivisions-1)', (1:numDivisions)' ].*(findClustersStruct.image_width/numDivisions);
MLlimits.y = [ (0:numDivisions-1)', (1:numDivisions)' ].*(findClustersStruct.image_height/numDivisions);
% 1) read in molecule list

% check to see if a molecule list was input
if isa(findClustersStruct.i3file, 'Insight3')
    % then input is actually a molecule list; change the handle reference
    % from 'i3file' to 'i3'
    i3file = findClustersStruct.i3file;
elseif ischar(findClustersStruct.i3file) && (~isempty(strfind(findClustersStruct.i3file,'.bin')) || ~isempty(strfind(findClustersStruct.i3file,'.txt')))
    % then input is correctly formatted ML file location, read it
    if exist(findClustersStruct.i3file, 'file') ~= 2
        error('File does not exist\n\t%s', findClustersStruct.i3file);
    else
        i3file = Insight3(findClustersStruct.i3file);
    end
else
    % then the user input something for i3file with an unknown format
    error('input ''findClustersStruct.i3file'' is of unrecognized format')
end

if viewDivisions
    figure;
    i3file.show();
    hold on;
    counter = 1;
    for ix = 1:numDivisions
        for iy = 1:numDivisions
            plot([MLlimits.x(ix,1), MLlimits.x(ix,2), MLlimits.x(ix,2), MLlimits.x(ix,1)], [MLlimits.y(iy,1), MLlimits.y(iy,1), MLlimits.y(iy,2), MLlimits.y(iy,2)]);
            %t = sprintf('%d', counter);
            %text(MLlimits.x(ix,1)+5,MLlimits.y(iy,1)+15, ['\color{magenta} ' t]);
            counter = counter + 1;
        end
    end
    ask = questdlg('Is this grid size okay?','Segment Molecule List and Find Clusters', 'Yes', 'No', 'Yes');
    if strcmp('No', ask);
        return
    end
    pause(0.2); % wait for the dialog box to close properly
    drawnow();
end

if findClustersStruct.use_drift_corrected_xy == true
    xcol = i3file.getColumnIndex('xc');
    ycol = i3file.getColumnIndex('yc');
else
    xcol = i3file.getColumnIndex('x');
    ycol = i3file.getColumnIndex('y');
end
% keep a copy of the original data
origML = i3file.getData();
origFileName = i3file.filename;
% reset ML to contain localizations only within the limits
counter = 1;
for ix = 1:numDivisions
    boxDatax = origML( (origML(:,xcol) > MLlimits.x(ix,1)), :);
    boxDatax = boxDatax( (boxDatax(:,xcol) <= MLlimits.x(ix,2)), :);
    
    for iy = 1:numDivisions
        boxDataxy = boxDatax( (boxDatax(:,ycol) > MLlimits.y(iy,1)), :);
        boxDataxy = boxDataxy( (boxDataxy(:,ycol) <= MLlimits.y(iy,2)), :);
        
        if size(boxDataxy,1) > 1
            % update the i3file data table
            i3file.setData( boxDataxy );
            % update the i3file filename field, since it's used to save a new
            % .bin file after cluster analysis
            appendum = ['_GridX' num2str(round(MLlimits.x(ix,1))) 'to' num2str(round(MLlimits.x(ix,2))) ...
                'Y' num2str(round(MLlimits.y(iy,1))) 'to' num2str(round(MLlimits.y(iy,2)))];
            i3file.setFilename( [i3file.filename(1:end-4) appendum '.bin'] );
            
            % update the user on the status of the box analysis
            disp(['Working on box #' num2str(counter) ' of ' ...
                num2str(numDivisions^2) ' total boxes'])
            
            %% perform the cluster analysis on the localizations in the box
            % results column names:
            % X(pix), Y(pix), Number_of_Loc_per_Cluster, SigX(pix), SigY(pix),
            % Mean(sigX,sigY), sqrt(sigX^2+sigY^2), Z(nm), sigZ(nm)
            findClustersStruct.i3file = i3file;
            FindClusters(findClustersStruct);            
        else
            disp(['SKIPPING  box #' num2str(counter) ', TOO FEW LOCALIZATIONS INSIDE BOX'])
        end
        counter = counter + 1;
        
        % reset i3file filename
        i3file.setFilename( origFileName );
    end
end

% merge the reults from each grid
[pathstr,name,~] = fileparts(origFileName);
outPath = fullfile(pathstr, 'cluster_analysis');
fileList = recursiveFindFile(outPath,[name '_GridX*.bin']);
if isempty(fileList)
    return
end

i3Merge = Insight3(fullfile(outPath, [name sprintf('_segmented%d.bin', numDivisions)]));
xynFilename = fullfile(outPath, [name sprintf('_segmented%d.xyn', numDivisions)]);

fid1 = XYN(strrep(fileList{1}, '.bin', '.xyn'));
xynHeader = fid1.header;

xynData = [];
for i=1:length(fileList)
    i3temp = Insight3(fileList{i});
    i3Merge.appendData(i3temp.getData());
    delete(fileList{i});
    
    fname = [fileList{i}(1:end-4) '.xyn'];
    f = XYN(fname);
    xynData = [xynData; f.data];
    delete(fname);    
end
i3Merge.write();

% create the XYN file
fle = fopen(xynFilename, 'w');
medN = median(xynData(:,3));
fprintf(fle, 'medianNum\t%d\n', medN);
medArea = pi*(median(xynData(:,6))*findClustersStruct.original_pixel_size)^2;
fprintf(fle, 'medianArea[nm^2]\t%f\n', medArea);
fprintf(fle, 'medianDensity\t%f\n', medN/medArea);
fprintf(fle, '%s\n', strjoin(xynHeader,'\t'));
fclose(fle);
dlmwrite(xynFilename, xynData, '-append', 'precision', '%.10g', 'delimiter', '\t');

end
