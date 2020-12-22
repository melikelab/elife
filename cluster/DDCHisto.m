% function fout = DDCHisto(ddcFile12, maxDistance12, ddcFile13, maxDistance13)
%
% Reads two DDC files to find the optimal rotation that aligns the slope
% of the line passing through the localizations for 3 different channels 
% (typically channels 1, 2, 3 that are the nearest neighbour clusters) 
% to be parallel to the x-axis to be able to view histograms to 
% understand the distances between the clusters.
%
% For each cluster-distance ('Distances XYZ (nm)'), in ddcFile12,
% that is less than maxDistance12, find the corresponding centroid in 
% ddcFile13 that has the same centroid for the File1 parameter. If the
% cluster-distance between File1 and File2 in ddcFile13 is less than
% maxDistance13 then perform a rotation of the localizations 
% within each cluster so that the slope of a line passing through all 
% localizations in a cluster is zero and parallel to the x-axis. A 
% weighted average rotation angle is used to perform a final rotation of the 
% localizations in all 3 clusters (using the number of localizations in each 
% cluster as the weights) about a common origin and the rotated localizations 
% are saved to a csv file with the centroid position and the angle of 
% rotation included in the output filename.
%
% the csv file has no header but is saved as:
% x(Ch1), y(Ch1), x(Ch2), y(Ch2), x(Ch3), y(Ch3)
%
% Inputs
% ------
% ddcFile12 : string or DDC object
%   if a string then the full path to where the first ddc file is located
%   otherwise it is a DDC object
%   This would be the file that was created using the DistanceDualColor
%   function analysing Insight3 localizations in channel 1 and channel 2
% maxDistance12 : double
%   if a cluster distance, 'Distances XYZ (nm)', in ddcFile12 is geater than
%   this value then ignore these centroids and go to the next centroid pair in ddcFile12
% ddcFile13 : string or DDC object
%   if a string then the full path to where the second ddc file is located
%   otherwise it is a DDC object
%   This would be the file that was created using the DistanceDualColor
%   function analysing Insight3 localizations in channel 1 and channel 3
%   NOTE: the first file 'File1' must be equal to 'File1' in ddcFile12
% maxDistance13 : double
%   if a cluster distance, 'Distances XYZ (nm)', in ddcFile13 is geater than
%   this value then ignore these centroids and go to the next centroid pair in ddcFile12
%
% Returns
% -------
% a Nx1 cell of the file names that we created
%
function fnames = DDCHisto(ddcFile12, maxDistance12, ddcFile13, maxDistance13)

% specify the number of nm in a pixel
nm_per_pixel = 160.0;

% read the first DDC file
if isa(ddcFile12, 'DDC')
    ddc1 = ddcFile12;
else
    ddc1 = DDC(ddcFile12);
end

% read the second DDC file
if isa(ddcFile13, 'DDC')
    ddc2 = ddcFile13;
else
    ddc2 = DDC(ddcFile13);
end

% make sure that the File1 value is the same for each ddc file
assert(strcmp(ddc1.file1, ddc2.file1), 'The File1 value must be the same in the two ddc files'); 

% get the indices of the columns that we are interested in
xyzCol = ddc1.getColumnIndex('Distances XYZ (nm)');
x1Col = ddc1.getColumnIndex('X(File1)');
y1Col = ddc1.getColumnIndex('Y(File1)');
n1Col = ddc1.getColumnIndex('#Locs in Cluster (File1)');
x2Col = ddc1.getColumnIndex('X(File2)');
y2Col = ddc1.getColumnIndex('Y(File2)');
n2Col = ddc1.getColumnIndex('#Locs in Cluster (File2)');

% read the molecule lists
bin1 = Insight3([ddc1.file1(1:end-4) '.bin']);
assert(~isempty(bin1.data), ['The molecule list does not exist ' bin1.filename]);
bin2 = Insight3([ddc1.file2(1:end-4) '.bin']);
assert(~isempty(bin2.data), ['The molecule list does not exist ' bin2.filename]);
bin3 = Insight3([ddc2.file2(1:end-4) '.bin']);
assert(~isempty(bin3.data), ['The molecule list does not exist ' bin3.filename]);

% get the columns that the xc, yc and channel information are in
xcCol = bin1.getColumnIndex('xc');
ycCol = bin1.getColumnIndex('yc');
chCol = bin1.getColumnIndex('channel');

% get the total number of cluster-distance pairs
numIter = 0;
for i = 1:length(ddc1.data)
    numIter = numIter + size(ddc1.data{i}, 1);
end

cnt = 0;
numCreated = 0;
fprintf('Analysing...   0.0%% complete');

origLocsWithCh3 = [];
origLocsWithoutCh3 = [];

fnames = {};
% for each distance category
for i = 1:length(ddc1.data)
    
    % for each cluster distance
    for j = 1:size(ddc1.data{i}, 1)
        
        fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*cnt/numIter)]);
        cnt = cnt + 1;
        
        % ignore distances greater than maxDistance
        if ddc1.data{i}(j,xyzCol) > maxDistance12
            continue
        end
        
        % find the localizations associated with each 
        % cluster centroid in the first molecule list
        x = ddc1.data{i}(j,x1Col);
        y = ddc1.data{i}(j,y1Col);
        n1 = ddc1.data{i}(j,n1Col);
        [centroid1, locs1] = getLocalizations(x, y, n1, bin1, chCol, xcCol, ycCol);
        
        % find the localizations associated with each 
        % cluster centroid in the third molecule list
        [x, y, n3, dist] = findXYN(ddc2, centroid1, x1Col, y1Col, x2Col, y2Col, n2Col, xyzCol);
        ignoreCh3 = dist > maxDistance13;
        if ~ignoreCh3
            [centroid3, locs3] = getLocalizations(x, y, n3, bin3, chCol, xcCol, ycCol);
        else
            centroid3 = zeros(2,1);
        end
        
        % find the localizations associated with each 
        % cluster centroid in the second molecule list
        x = ddc1.data{i}(j,x2Col);
        y = ddc1.data{i}(j,y2Col);
        n2 = ddc1.data{i}(j,n2Col);
        [centroid2, locs2] = getLocalizations(x, y, n2, bin2, chCol, xcCol, ycCol);
        
        % calculate the weighted average centroid that the localizations will be
        % rotated about
        centroid = zeros(2,1);
        for ic=1:2
            centroid(ic) = (centroid1(ic)*n1 + centroid2(ic)*n2 + centroid3(ic)*n3) / (n1 + n2 + n3);
        end

        % find the rotation angle that makes the slope zero for each cluster
        % take a weighted average, where the weights are the number of
        % localizations in the cluster.
        theta1 = getTheta(centroid, locs1(:,[xcCol ycCol]));
        theta2 = getTheta(centroid, locs2(:,[xcCol ycCol]));
        if ignoreCh3
            theta3 = 0.0;
            n3 = 0;
        else
            theta3 = getTheta(centroid, locs3(:,[xcCol ycCol]));
        end
        %theta = (theta1 + theta2 + theta3)/3.0; %straight average
        theta = (theta1*n1 + theta2*n2 + theta3*n3)/(n1 + n2 + n3); % weighted ave
        
        % rotate each molecule list by the average theta value
        locs1rot = rotateLocs(theta, centroid, locs1(:,[xcCol ycCol]));
        locs2rot = rotateLocs(theta, centroid, locs2(:,[xcCol ycCol]));
        if ignoreCh3
            locs3rot = [];
        else
            locs3rot = rotateLocs(theta, centroid, locs3(:,[xcCol ycCol]));
        end
        
        % save the rotated molecule lists to a file
        fout = getOutputFilename(ddc1.filename, centroid, theta, ignoreCh3);
        fnames{cnt} = fout;
        try
            saveRotatedLocs(fout, locs1rot, locs2rot, locs3rot, nm_per_pixel);
            numCreated = numCreated + 1;
        catch
            disp('\n');
            error('Cannot save %s\nIs it currrently open in another program?', fout);
        end

        if ignoreCh3
            locs1(:,12) = 1;
            origLocsWithoutCh3 = [origLocsWithoutCh3; locs1];
            locs2(:,12) = 2;
            origLocsWithoutCh3 = [origLocsWithoutCh3; locs2];
        else
            locs1(:,12) = 1;
            origLocsWithCh3 = [origLocsWithCh3; locs1];
            locs2(:,12) = 2;
            origLocsWithCh3 = [origLocsWithCh3; locs2];
            locs3(:,12) = 3;
            origLocsWithCh3 = [origLocsWithCh3; locs3];
        end

        % plot the initial and rotated molecule lists
        %plotIt(locs1(:,[xcCol ycCol]), locs2(:,[xcCol ycCol]));
        %plotIt(locs1rot, locs2rot);
        %break
    end
    %break
end

fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete\n',100*cnt/numIter)]);
if numCreated == 1
    fprintf('Created %d csv file in %s\n', numCreated, fileparts(fout));
else
    fprintf('Created %d csv files in %s\n', numCreated, fileparts(fout));
end

saveOriginalLocs(ddc1.filename, origLocsWithoutCh3, origLocsWithCh3);

end


function [centroid, locs] = getLocalizations(x, y, n, bin, chCol, xcCol, ycCol)
% get the localizations in this cluster
    locs = [];
    % for each localization
    for k = 1:size(bin.data, 1)
        % the centroids are in channel 0 and the localizations that
        % are part of a cluster immediately follow in the molecule list
        if bin.data(k, chCol) == 0 && abs(bin.data(k, xcCol) - x) < 1e-2 && abs(bin.data(k, ycCol) - y) < 1e-2
            %fprintf('%d %f %f\n', bin.data(k, chCol), abs(bin.data(k, xcCol) - x), abs(bin.data(k, ycCol) - y))
            centroid = bin.data(k,[xcCol ycCol]);
            locs = bin.data(k+1:k+n,:);            
            break
        end
    end
    assert(~isempty(locs), 'Cannot find the centroid in the molecule list');
end

function [slope, intercept] = getSlopeIntercept(locs)
% fit the (x,y) coords of the localizations to a 
% linear function to determine the slope and intercept
    fitvars = polyfit(locs(:,1), locs(:,2), 1);
    slope = fitvars(1);
    intercept = fitvars(2);
end

function plotIt(loc1, loc2)
    [m1, b1] = getSlopeIntercept(loc1);
    [m2, b2] = getSlopeIntercept(loc2);
    min1 = min(loc1);
    max1 = max(loc1);
    min2 = min(loc2);
    max2 = max(loc2);
    x = linspace(min(min1(1), min2(1)), max(max1(1), max2(1)));
    y1 = m1*x+b1;
    y2 = m2*x+b2;
    figure 
    hold on
    plot(x, y1, 'r');
    plot(loc1(:,1), loc1(:,2), 'r+')
    plot(x, y2, 'b');
    plot(loc2(:,1), loc2(:,2), 'b+')
    legend('1', '1', '2', '2')
    axis('image')
end

function locs = rotateLocs(theta, centroid, locs)
% rotate the molecule list
    xc = centroid(1);
    yc = centroid(2);
    c = cos(theta);
    s = sin(theta);
    for m=1:size(locs,1)
        x = locs(m,1) - xc;
        y = locs(m,2) - yc;
        xp = c*x - s*y;
        yp = s*x + c*y;
        locs(m,1) = xp + xc;
        locs(m,2) = yp + yc;
    end        
end

function theta = getTheta(centroid, locs)
% determine the rotation angle that makes the slope zero
    locsCopy = locs;
    slope = Inf; theta = 0.0;
    while abs(slope) > 1e-8
        [slope, ~] = getSlopeIntercept(locsCopy);
        angle = -atan(slope);
        locsCopy = rotateLocs(angle, centroid, locsCopy);
        theta = theta + angle;
    end
end

function saveRotatedLocs(fname, loc1, loc2, loc3, nm_per_pixel)
% save the rotated molecule lists
   n1 = size(loc1, 1);
   n2 = size(loc2, 1);
   n3 = size(loc3, 1);
   n = max(n1, max(n2, n3));
   
   fp = fopen(fname, 'w');
   for i = 1:n
       if i <= n1
           text = sprintf('%f,%f,', nm_per_pixel*loc1(i,1), nm_per_pixel*loc1(i,2));
       else
           text = ',,';
       end
       if i <= n2
           text = strcat(text, sprintf('%f,%f,', nm_per_pixel*loc2(i,1), nm_per_pixel*loc2(i,2)));
       else
           text = strcat(text, ',,');
       end
       if ~isempty(loc3) && i <= n3
           text = strcat(text, sprintf('%f,%f', nm_per_pixel*loc3(i,1), nm_per_pixel*loc3(i,2)));
       else
           text = strcat(text, ',');
       end
       fprintf(fp,'%s\n', text);
   end
   fclose(fp);
end

function saveOriginalLocs(ddc1, origLocsWithoutCh3, origLocsWithCh3)
[path, name, ~] = fileparts(ddc1);

i3with = Insight3([path '/ddcHisto/' name '_withCh3.bin']);
i3with.setData(origLocsWithCh3);
i3with.forceFileOverwrite(true);
i3with.write();

i3without = Insight3([path '/ddcHisto/' name '_withoutCh3.bin']);
i3without.setData(origLocsWithoutCh3);
i3without.forceFileOverwrite(true);
i3without.write();

end

function fout = getOutputFilename(fname, centroid, theta, ignoreCh3)
    [path, name, ~] = fileparts(fname);
    
    subFolder = 'ddcHisto';
    if ignoreCh3
        subSubFolder = 'withoutCh3';
    else
        subSubFolder = 'withCh3';
    end
    
    if ~isequal(exist(strcat(path, '/', subFolder, '/'), 'dir'), 7)
        mkdir(strcat(path, '/', subFolder, '/'));        
    end
    
    if ~isequal(exist(strcat(path, '/', subFolder, '/', subSubFolder, '/'), 'dir'), 7)
        mkdir(strcat(path, '/', subFolder, '/', subSubFolder, '/'));
    end        
    
    fout = fullfile(path, subFolder, subSubFolder, [name sprintf('_x%.2f_y%.2f_phi%.5f.csv', centroid(1), centroid(2), theta)]);
end

function [x, y, n, d] = findXYN(ddc, centroid, x1Col, y1Col, x2Col, y2Col, n2Col, xyzCol)

x = Inf;
y = Inf;
n = Inf;
d = Inf;

% for each distance category
for i = 1:length(ddc.data)    
    % for each cluster distance
    for j = 1:size(ddc.data{i}, 1)        
        if abs(ddc.data{i}(j,x1Col) - centroid(1)) < 1e-2 && abs(ddc.data{i}(j,y1Col) - centroid(2)) < 1e-2
            x = ddc.data{i}(j,x2Col);
            y = ddc.data{i}(j,y2Col);
            n = ddc.data{i}(j,n2Col);
            d = ddc.data{i}(j,xyzCol);
            break
        end        
    end    
    if x ~= Inf
        break
    end
end

if x == Inf
    assert(x ~= Inf, 'Cannot find the centroid in the ddc file');
end

end
