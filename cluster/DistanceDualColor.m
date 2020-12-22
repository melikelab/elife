% function [cat, counts] = DistanceDualColor(file1, file2, categories, pixel_size_nm, minimum_localizations, sortResults)
%
% Calculate the distance between 2 cluster files and places the clusters in different
% categories (eg., SYN, PERI, EXTRA) based on the distance between the clusters. 
% The values in the 'categories' array represent the different distance ranges to 
% use to categorize the clusters.
%
% Inputs
% ------
%
% file1 : string or XYN object
%     The first cluster file
%     If string then the path to the XYN file
%
% file2 : string or XYN object
%     The second cluster file. Can be the same as file1
%     If string then the path to the XYN file
%
% categories : 1D array
%     The values, in nm, specify the different categories to separate 
%     the clusters into based on their separation distance. For example, [70 250] means 
%     that you want to have 3 different distance categories 
%     (i) distance < 70 nm
%     (ii) 70 < distance < 250 nm
%     (iii) distance > 250 nm
%
% pixel_size_nm : double
%     The size of a pixel, in nm
%
% minimum_localizations : int
%     Each cluster must have >= 'minimum_localizations' number of localizations in the
%     cluster to be considered as a cluster to calculate the distance for
%
% sortResults : boolean
%     Sort the output values by increasing distance
%
% Returns
% -------
% cat : struct of 2D arrays
%     Each struct is an Nx11 array for all the clusters in the different distance ranges
%     For example, if 'categories' = [70 250] (see above) then there will be 3
%     2D arrays within the struct representing the 3 different distance
%     ranges. The first 2D array, cat{1}, is a summary of all clusters that are 
%     within the first distance range (dist<70nm), the second 2D array, cat{2},
%     is a summary of all clusters that are within the second distance
%     range (70<dist<250) and the third 2D array, cat{3}, is a summary of all clusters 
%     that are within the third distance range (dist>250). 
%
%     The columns in the 2D arrays are:
%     Column 1 -> the distance to the nearest cluster, in nm, in the XY dimension
%     Column 2 -> the distance to the nearest cluster, in nm, in the XYZ dimension
%     Column 3 -> the cluster area for file-1, in nm^2
%     Column 4 -> the number of localizations in the file-1 cluster
%     Column 5 -> the cluster area for file-2, in nm^2
%     Column 6 -> the number of localizations in the file-2 cluster
%     Column 7 -> the X position of the cluster centroid for file-1, in pixels
%     Column 8 -> the Y position of the cluster centroid for file-1, in pixels
%     Column 9 -> the minimum of the X,Y sigma value for each file-1 cluster, in nm
%     Column 10 -> the Z value, in nm, of each cluster of file-1
%     Column 11 -> the sigma Z value, in nm, of each cluster of file-1
%     Column 12 -> the X position of the cluster centroid for file-2, in pixels
%     Column 13 -> the Y position of the cluster centroid for file-2, in pixels
%     Column 14 -> the minimum of the X,Y sigma value for each file-2 cluster, in pixels
%     Column 15 -> the Z value, in nm, of each cluster of file-2
%     Column 16 -> the sigma Z value, in nm, of each cluster of file-2
%
% counts : 1xM array
%     The first M-2 columns are the percentage of clusters that are within each distance range.
%     The last 2 columns are the total number of clusters within each color file.
%     For example, if 'categories' = [70 250] and counts returns [75.2 22.8 2.0 423 561] then this
%     means that 75.2% of the clusters had a distance < 70nm, 22.8% were
%     between 70 and 250nm, 2.0% had a distance greater than 250nm, there
%     were 423 clusters in color-1 and 561 clusters in color-2
function [cat, counts] = DistanceDualColor(file1, file2, categories, pixel_size_nm, minimum_localizations, sortResults)

if ~isa(file1, 'XYN')
    if ~ischar(file1)
        error('the first parameter, file1, must be a character array (i.e., a String) or an XYN object')
    end
    file1 = XYN(file1);
end
if ~isa(file2, 'XYN')
    if ~ischar(file2)
        error('the second parameter, file2, must be a character array (i.e., a String) or an XYN object')
    end
    file2 = XYN(file2);
end

isSameFile = strcmp(file1.filename, file2.filename);
data1 = file1.data(( file1.data(:,3) > minimum_localizations ),:);
data2 = file2.data(( file2.data(:,3) > minimum_localizations ),:);

noZValues = size(data1, 2) < 8 || size(data2, 2) < 8;

% calculate the distances
n1 = size(data1, 1);
n2 = size(data2, 1);
n = n1*n2;
out = zeros(n1,16);
dist = zeros(n2,16);
cnt = 0;
fprintf('Distance calculation is   0.0%% complete');
for i = 1:n1
    fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*cnt/n)]);

    for j = 1:n2
        
        cnt = cnt + 1;
        
        if isSameFile && i==j
            % ignore the expected distance of 0 for the same cluster
            for k=1:size(out,2)
                dist(j,k) = Inf;
            end
        else
            dist(j,1) = sqrt( (data1(i,1)-data2(j,1))^2 + (data1(i,2)-data2(j,2))^2);
            if noZValues
                dist(j,2) = dist(j,1);
            else
                dist(j,2) = sqrt( ((data1(i,1)-data2(j,1))*pixel_size_nm)^2 + ((data1(i,2)-data2(j,2))*pixel_size_nm)^2 + (data1(i,8)-data2(j,8))^2);
            end
            dist(j,3) = data1(i,6);
            dist(j,4) = data2(j,6);
            dist(j,5) = data1(i,1);
            dist(j,6) = data1(i,2);
            dist(j,7) = min(data1(i,4),data1(i,5));
            dist(j,8) = data1(i,3);
            if noZValues
                dist(j,9) = 0.0;
                dist(j,10) = 0.0;
            else
                dist(j,9) = data1(i,8);
                dist(j,10) = data1(i,9);
            end
            dist(j,11) = data2(j,1);
            dist(j,12) = data2(j,2);
            dist(j,13) = min(data2(j,4),data2(j,5));
            dist(j,14) = data2(j,3);
            if noZValues
                dist(j,15) = 0.0;
                dist(j,16) = 0.0;
            else
                dist(j,15) = data2(j,8);
                dist(j,16) = data2(j,9);
            end
        end
    end
    
    [~,I]=sort(dist(:,1));
    srt=dist(I,:);
    
    % get the closest cluster
    out(i,1) = srt(1,1)*pixel_size_nm;
    out(i,2) = srt(1,2);
    out(i,3) = pi * (srt(1,3)*pixel_size_nm)^2;
    out(i,4) = srt(1,8);
    out(i,5) = pi * (srt(1,4)*pixel_size_nm)^2;
    out(i,6) = srt(1,14);
    out(i,7) = srt(1,5);
    out(i,8) = srt(1,6);
    out(i,9) = srt(1,7)*pixel_size_nm;    
    out(i,10) = srt(1,9);
    out(i,11) = srt(1,10);
    out(i,12) = srt(1,11);
    out(i,13) = srt(1,12);
    out(i,14) = srt(1,13)*pixel_size_nm;
    out(i,15) = srt(1,15);
    out(i,16) = srt(1,16);
end
fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*cnt/n)]);

% sort the results
if sortResults
    [~,I]=sort(out(:,1));
    out=out(I,:);
end
fprintf('\n');

fprintf('Separating the clusters into the distance categories... ');

if categories(1) ~= 0
    categories = [0 categories];
end
categories = [categories Inf];

ni = size(out, 1);
nj = length(categories)-1;
counts = zeros(1,nj);
cat = cell(nj, 1);
for i = 1:ni    
    for j = 1:nj
        if (out(i,1) > categories(j)) && (out(i,1) < categories(j+1))            
            counts(j) = counts(j) + 1;
            cat{j} = [cat{j}; out(i,:)];
            break
        end
    end
end
fprintf('DONE\n');

[path, fname, ~] = fileparts(file1.filename);

outFileName = fullfile(path,strcat(fname,'.ddc'));

rangestr = '';
rangestr2 = cell(1,nj);
for j = 1:nj
    rangestr = strcat(rangestr, sprintf('%% in %d-%d nm,', categories(j), categories(j+1)));
    rangestr2{j} = sprintf('median in %d-%d nm,', categories(j), categories(j+1));
    counts(j) = 100*counts(j)/ni;
end
counts(nj+1) = size(data1, 1);
counts(nj+2) = size(data2, 1);

try
    fprintf('Creating output file... ');
    
    fle = fopen(outFileName, 'w');
    fprintf(fle, 'File1:\t%s\n', file1.filename);
    fprintf(fle, 'File2:\t%s\n', file2.filename);
    v = strjoin(strsplit(rangestr(1:end-1), ','), '\t');
    fprintf(fle, '%s\t%s\t%s\n', v, '#clusters (File1)', '#clusters (File2)');
    fclose(fle);

    dlmwrite(outFileName, counts, '-append', 'delimiter', '\t');
    
    header = {'Distances XY (nm)', 'Distances XYZ (nm)', 'Cluster Area File1 (nm²)', '#Locs in Cluster (File1)', 'Cluster Area File2 (nm²)', '#Locs in Cluster (File2)', 'X(File1)', 'Y(File1)', 'min Sigma(File1 nm)', 'Z(File1 nm)', 'sigZ(File1 nm)', 'X(File2)', 'Y(File2)', 'min Sigma (File2 nm)', 'Z(File2 nm)', 'sigZ(File2 nm)'};
    
    fle = fopen(outFileName, 'a');
    fprintf(fle, '%s\n', strjoin(header, '\t'));

    % calculate the median for each column, for each range
    for j = 1:nj
        s = '';
        n = size(cat{j},2);
        for i=1:n
            s = strcat(s, sprintf('%g,', median(cat{j,1}(:,i))));
        end
        if isempty(s)
            for i=1:size(out,2)
                s = strcat(s, '0,');
            end
        end
        v = strrep(s, ',', '\t');
        fprintf(fle, '%s\t%s\n', sprintf(v(1:end-2)), rangestr2{j}(1:end-1));
    end   

    % calculate the median for each column for all ranges
    s = '';
    n = size(out,2);
    for i=1:n
        column = [];
        for j = 1:nj
            if ~isempty(cat{j})
                column = [column; cat{j}(:,i)];
            end
        end
        s = strcat(s, sprintf('%g,', median(column)));
    end
    v = strrep(s, ',', '\t');
    fprintf(fle, '%s\t%s\n', sprintf(v(1:end-2)), 'median in entire range');

    % write another header
    fprintf(fle, '%s\n', strjoin(header, '\t'));
    fclose(fle);

    for j = 1:nj
        dlmwrite(outFileName, cat{j}, '-append', 'delimiter', '\t');
        fle = fopen(outFileName, 'a');
        fprintf(fle, '\n');
        fclose(fle);
    end

    fprintf('DONE\n');
    fprintf('Saved to %s\n', outFileName);
    
catch
    fprintf('\n');
    error('Cannot write to %s\nIs this file open in another program?\n', outFileName);
end

end



