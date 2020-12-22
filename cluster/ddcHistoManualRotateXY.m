% function [dataRotated, gaussianCenters] = ddcHistoManualRotateXY(fname, radians, nBins, showPlots, saveRotatedLocs, saveFitDistances)
%
% Rotate the (x,y) values in 'fname' by the angle 'radians' and 
% display and/or save the results.
%
% Inputs
% ------
% fname : string
%   a filename created from the DDCHisto function
%
% radians : double
%   the angle, in radians, to rotate the localizations by
%
% nBins : integer
%   the number of bins to create for plotting the histograms
%   if the value is <1, NaN, Inf then nBins = round(sqrt(# of data points))
%
% showPlots : boolean
%   show the figure window
%
% saveRotatedLocs : boolean
%   save the rotated localizations to a file. The filename will be the same
%   as 'fname' except the value of phi will be updated.
%
% saveFitDistances : boolean
%   save the fit results to a file. The filename will be the same
%   as 'fname' except the value of phi will be updated and the string
%   '_fitDistances' will be appended.
%
% Returns
% -------
% dataRotated: The rotated results, a N x 1 cell of rotated values.
%
% fitDistances: histogram centers deteremined from a fit to the rotated localizations
%   if there are three channels then the structure of the output is
%   x[pix], y[pix], abs(ch2-ch1)[nm], ch3-ch1[nm], ch3-ch2[nm], inside, sigma1, sigma2, sigma3, uncert_c1, uncert_c2, uncert_c3, uncert_sigma1, uncert_sigma2, uncert_sigma3
%
%   if there are two channels then the structure of the output is
%   x[pix], y[pix], abs(ch2-ch1)[nm], sigma1, sigma2, uncert_c1, uncert_c2, uncert_sigma1, uncert_sigma2
%
function [dataRotated, fitDistances] = ddcHistoManualRotateXY(fname, radians, nBins, showPlots, saveRotatedLocs, saveFitDistances)
    
    % temporarily disable warning about the legends
    warning('off','MATLAB:legend:IgnoringExtraEntries')

    assert(isa(fname, 'char'), 'fname must be a string')
    data = dlmread(fname, ',');
    
    % if the csv file was opened in Excel then there is a possiblity
    % the there are empty Excel columns next the the data... 
    % ignore the empty columns
    nCols = 0;
    for i = 1:size(data, 2)
        indices = find(data(:,i) ~= 0);
        if ~isempty(indices)
            nCols = nCols + 1;
        end
    end
    assert(nCols == 6 || nCols == 4, 'The number of columns in the file must be 4 or 6')
    nBins = round(nBins);
    
    gaussianCenters = []; % the centers (mean) of the Gaussian fits
    
    % calculate the common centroid for all localizations
    xCentroids = zeros(nCols/2,2); 
    yCentroids = zeros(nCols/2,2);
    k = 1;
    for j=1:2:nCols
        indices = find(data(:,j) ~= 0);
        xCentroids(k,1) = mean(data(indices,j));
        xCentroids(k,2) = length(indices);
        yCentroids(k,1) = mean(data(indices,j+1));
        yCentroids(k,2) = length(indices);
        k = k + 1;
    end
    xc = 0.0;
    nx = 0;
    yc = 0.0;
    ny = 0;
    for j=1:size(xCentroids,1)
        xc = xc + xCentroids(j,1) * xCentroids(j,2);
        nx = nx + xCentroids(j,2);
        yc = yc + yCentroids(j,1) * yCentroids(j,2);
        ny = ny + yCentroids(j,2);
    end
    xc = xc / nx;
    yc = yc / nx;

    dataRotated = cell(nCols/2,1);
    c = cos(radians);
    s = sin(radians);
    k = 1;
    for j=1:2:nCols
        indices = find(data(:,j) ~= 0);
        %xc = mean(data(indices,j));
        %yc = mean(data(indices,j+1));
        temp = zeros(length(indices), 1);
        for i=1:length(indices)
            idx = indices(i);
            x = data(idx,j) - xc;
            y = data(idx,j+1) - yc;
            temp(idx,1) = c*x - s*y + xc;
            temp(idx,2) = s*x + c*y + yc;
        end
        dataRotated{k} = temp;
        k = k + 1;
    end
    
    colors = ['r' 'b' 'g'];
    if showPlots % plot the results        
        figure;        
    end
    
    % original, localizations
    
    if showPlots
        subplot(2,2,1);
        hold on;
    end
    yvals = cell(nCols/2,1);
    j = 1;
    for i = 1:2:nCols
        indices = find(data(:,i) ~= 0);
        yvals{j} = data(indices,i+1);
        if showPlots
            plot(data(indices,i), data(indices,i+1), [colors(j) '+']);
        end
        j = j + 1;
    end
    if showPlots
        title('Original');
        if nCols == 6
            legend('ch1', 'ch2', 'ch3', 'Location','northoutside','Orientation','horizontal');
        else
            legend('ch1', 'ch2', 'Location','northoutside','Orientation','horizontal');
        end
    end
    
    
    % original, histograms
    if showPlots
        subplot(2,2,3);
        hold on;
    end
    plotHandle = [];
    
    for i = 1:length(yvals)        
        if nBins < 1 || isnan(nBins) || isinf(nBins)
            [n, xout] = hist(yvals{i} , linspace(min(yvals{i}), max(yvals{i}), round(sqrt(length(yvals{i})))));
            [params, xfit, yfit, uncerts] = fitBestNumber(yvals{i}, round(sqrt(length(yvals{i}))), xout);
        else
            [n, xout] = hist(yvals{i} , linspace(min(yvals{i}), max(yvals{i}), nBins));
            [params, xfit, yfit, uncerts] = fitBestNumber(yvals{i}, nBins, xout);
        end
        %[params, xfit, yfit, uncerts] = fit(xout, n);        
        if isempty(params)
            gaussianCenters = [gaussianCenters; [NaN NaN NaN NaN]];
        else
            gaussianCenters = [gaussianCenters; [params(2) params(4) uncerts(2) uncerts(4)]];
        end
        if showPlots
            bar(xout, n, colors(i));
            if ~isempty(xfit)
                h = plot(xfit, yfit, 'k');
                set(h(1),'linewidth',3);
                plotHandle = [plotHandle; h];
            end
        end
    end
    if showPlots
        string = '';
        for i = 1:length(yvals)
            if nBins < 1 || isnan(nBins) || isinf(nBins)
                string = strcat(string, num2str(round(sqrt(length(yvals{i})))), ',');
            else
                string = [num2str(nBins) ','];
                break;
            end
        end
        title(sprintf('Original, #bins (%s)', string(1:end-1)));
        xlabel('Y values');
        ylabel('Counts');
        if nCols == 6
            c1 = gaussianCenters(1,2);
            c2 = gaussianCenters(2,2);
            c3 = gaussianCenters(3,2);
            legend(plotHandle, sprintf('|2-1|: %.2f nm',abs(c2-c1)), ...
                               sprintf('3-1: %.2f nm',c3-c1), ...
                               sprintf('3-2: %.2f nm',c3-c2))
        else
            legend(plotHandle, sprintf('|2-1|: %.2f',abs(gaussianCenters(2,2)-gaussianCenters(1,2))));
        end
    end
    
    % rotated, localizations
    if showPlots
        subplot(2,2,2);
        hold on;
    end
    yvals = cell(nCols/2,1);
    for j = 1:length(dataRotated)
        indices = find(dataRotated{j}(:,1) ~= 0);
        yvals{j} = dataRotated{j}(indices,2);
        if showPlots
            plot(dataRotated{j}(indices,1), dataRotated{j}(indices,2), [colors(j) '+']);
        end
    end
    if showPlots
        title(sprintf('Rotated (%f)', radians));
        if nCols == 6
            legend('ch1', 'ch2', 'ch3', 'Location','northoutside','Orientation','horizontal');
        else
            legend('ch1', 'ch2', 'Location','northoutside','Orientation','horizontal');
        end
    end
    
    % rotated, histograms
    if showPlots
        subplot(2,2,4);
        hold on;
    end
    plotHandle = [];
    for i = 1:length(yvals)
        if nBins < 1 || isnan(nBins) || isinf(nBins)
            [n, xout] = hist(yvals{i} , linspace(min(yvals{i}), max(yvals{i}), round(sqrt(length(yvals{i})))));
            [params, xfit, yfit, uncerts] = fitBestNumber(yvals{i}, round(sqrt(length(yvals{i}))), xout);
        else
            [n, xout] = hist(yvals{i} , linspace(min(yvals{i}), max(yvals{i}), nBins));
            [params, xfit, yfit, uncerts] = fitBestNumber(yvals{i}, nBins, xout);
        end
        %[params, xfit, yfit, uncerts] = fit(xout, n);
        if isempty(params)
            gaussianCenters = [gaussianCenters; [NaN NaN NaN NaN]];
        else
            gaussianCenters = [gaussianCenters; [params(2) params(4) uncerts(2) uncerts(4)]];
        end
        if showPlots
            bar(xout, n, colors(i));
            if ~isempty(xfit)
                h = plot(xfit, yfit, 'k');
                set(h(1),'linewidth',3);
                plotHandle = [plotHandle; h];
            end
        end
    end
    
    % get the x and y values from the filename
    idx1 = regexp(fname, '_x\d+.\d+_y\d+.\d+');
    idx2 = regexp(fname, '_phi\d+.\d+');
    if isempty(idx2)
        idx2 = regexp(fname, '_phi-\d+.\d+');
    end
    xyVal = str2num(strrep(fname(idx1+2:idx2-1), '_y', ' '));
    if nCols == 6
        c1 = gaussianCenters(4,2);
        c2 = gaussianCenters(5,2);
        c3 = gaussianCenters(6,2);
        sigma1 = gaussianCenters(4,1);
        sigma2 = gaussianCenters(5,1);
        sigma3 = gaussianCenters(6,1);
        uncert_sigma1 = gaussianCenters(4,3);
        uncert_sigma2 = gaussianCenters(5,3);
        uncert_sigma3 = gaussianCenters(6,3);
        uncert_c1 = gaussianCenters(4,4);
        uncert_c2 = gaussianCenters(5,4);
        uncert_c3 = gaussianCenters(6,4);
        if c1 < c2
            inside = (c3 > c1) && (c3 < c2);
        else
            inside = (c3 > c2) && (c3 < c1);
        end
        fitDistances = [xyVal(1) xyVal(2) abs(c2-c1) c3-c1 c3-c2 inside sigma1 sigma2 sigma3 uncert_c1 uncert_c2 uncert_c3 uncert_sigma1 uncert_sigma2 uncert_sigma3];
    else
        c1 = gaussianCenters(3,2);
        c2 = gaussianCenters(4,2);
        sigma1 = gaussianCenters(3,1);
        sigma2 = gaussianCenters(4,1);
        uncert_c1 = gaussianCenters(3,4);
        uncert_c2 = gaussianCenters(4,4);
        uncert_sigma1 = gaussianCenters(3,3);
        uncert_sigma2 = gaussianCenters(4,3);
        fitDistances = [xyVal(1) xyVal(2) abs(c2-c1) sigma1 sigma2 uncert_c1 uncert_c2 uncert_sigma1 uncert_sigma2];
    end

    if showPlots
        string = '';
        for i = 1:length(yvals)            
            if nBins < 1 || isnan(nBins) || isinf(nBins)
                string = strcat(string, num2str(round(sqrt(length(yvals{i})))), ',');
            else
                string = [num2str(nBins) ','];
                break;
            end
        end
        title(sprintf('Rotated, #bins (%s)', string(1:end-1)));
        xlabel('Y values');
        ylabel('Counts');    
        if nCols == 6
            legend(plotHandle, sprintf('|2-1|: %.2f nm',abs(c2-c1)), ...
                               sprintf('3-1: %.2f nm',c3-c1), ...
                               sprintf('3-2: %.2f nm',c3-c2))
        else
            legend(plotHandle, sprintf('|2-1|: %.2f',abs(c2-c1)));
        end
    end
    
    if saveRotatedLocs
        idx = regexp(fname, '_phi\d+.\d+');
        if isempty(idx)
            idx = regexp(fname, '_phi-\d+.\d+');
        end
        fname1 = [fname(1:idx) sprintf('phi%f.csv', radians)];
        if nCols == 4
            dataRotated{3} = {};
        end
        try
            saveLocalizations(fname1, dataRotated{1}, dataRotated{2}, dataRotated{3});
        catch
            error('Cannot save %s\nIs the file open?', fname1);
        end
    end

    if saveFitDistances
        [path, name, ext] = fileparts(fname);
        idx1 = strfind(name,'_');  
        idx2 = regexp(name, '_x\d+.\d+_y\d+.\d+');
        idx3 = regexp(name, '_phi\d+.\d+');
        if isempty(idx3)
            idx3 = regexp(name, '_phi-\d+.\d+');
        end
        fname2 = fullfile(path, [name(1:idx1(1)-1) name(idx2:idx3) sprintf('phi%f_fitDistances.csv', radians)]);
        try
            saveFits(fname2, fitDistances);
        catch
            error('Cannot save %s\nIs the file open?', fname2);
        end
    end
    
    warning('on','MATLAB:legend:IgnoringExtraEntries')
    
end

function saveFits(fname, values)
    fp = fopen(fname, 'w');
    if length(values) == 3
        fprintf(fp,'%s\n', 'x[pix],y[pix],abs(ch2-ch1)[nm],sigma1,sigma2,uncert_c1,uncert_c2,uncert_sigma1,uncert_sigma2');
        fprintf(fp,'%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8), values(9));
    else
        fprintf(fp,'%s\n', 'x[pix],y[pix],abs(ch2-ch1)[nm],ch3-ch1[nm],ch3-ch2[nm],inside,sigma1,sigma2,sigma3,uncert_c1,uncert_c2,uncert_c3,uncert_sigma1,uncert_sigma2,uncert_sigma3');
        fprintf(fp,'%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', values(1), values(2), values(3), values(4), values(5), values(6), values(7), values(8), values(9), values(10), values(11), values(12), values(13), values(14), values(15));
    end    
    fclose(fp);    
    fprintf('Created: %s\n', fname);
end

function saveLocalizations(fname, loc1, loc2, loc3)
% save the rotated molecule lists
   n1 = size(loc1, 1);
   n2 = size(loc2, 1);
   n3 = size(loc3, 1);
   n = max(n1, max(n2, n3));
   
   fp = fopen(fname, 'w');
   for i = 1:n
       if i <= n1
           text = sprintf('%f,%f,', loc1(i,1), loc1(i,2));
       else
           text = ',,';
       end
       if i <= n2
           text = strcat(text, sprintf('%f,%f,', loc2(i,1), loc2(i,2)));
       else
           text = strcat(text, ',,');
       end
       if ~isempty(loc3) && i <= n3
           text = strcat(text, sprintf('%f,%f', loc3(i,1), loc3(i,2)));
       else
           text = strcat(text, ',');
       end
       fprintf(fp,'%s\n', text);
   end
   fclose(fp);
   fprintf('Created: %s\n', fname);
end

function [params, xfit, yfit, uncerts] = fit(x, y)
% fit the histogram data to a Gaussian function

% create a Gaussian object for the fit
gauss = Gaussian1D();

% set the X data for the fit function
gauss.setX(x);

% set an initial guess
sigma = std(x);
area = max(y)*sigma*2;
gauss.setGuess([0 sigma area mean(x)]);

% don't float the background
gauss.setFloatParams([0 1 1 1]);

% create the LevMar fitter
levmar = LevMar(gauss, y);

% fit
exit_code = levmar.solve();

% check if fit converged
if exit_code > 0
    params = levmar.bestParams;
    uncerts = levmar.paramUncerts;
    xfit = linspace(min(x), max(x), 100);
    gauss.setX(xfit);
    yfit = gauss.eval(params);
else
    params = [];
    xfit = [];
    yfit = [];
    uncerts = [];
end

end

function [params, xfit, yfit, uncerts] = fitBestNumber(yvals, nBins, xout)

clustHist = ClusterHisto();
try
    fitResults2 = clustHist.fit(yvals, 2, nBins, false); % fit to 2 Gaussians
    heightA = fitResults2(3,1)/(fitResults2(2, 1) * sqrt(2*pi));
    heightB = fitResults2(7,1)/(fitResults2(6, 1) * sqrt(2*pi));
    fprintf(2, 'Not enough data points to fit to 2 Gaussians\n');
catch
    heightA = NaN;
    heightB = NaN;
end
maxHeightAB = max(heightA, heightB);
fitResults1 = clustHist.fit(yvals, 1, nBins, false); % fit to 1 Gaussian
height = fitResults1(3,1)/(fitResults1(2, 1) * sqrt(2*pi));
if ~isnan(maxHeightAB)
    if height > maxHeightAB
        % use 1 Gaussian
        params = fitResults1(:,1);
        uncerts = fitResults1(:,2);
        xfit = linspace(min(xout), max(xout), 100)';
        gauss = Gaussian1D();
        gauss.setX(xfit);
        yfit = gauss.eval(params);
    else
        % use 2 Gaussians        
        if heightA > heightB
            %fprintf('A %f, %f\n', heightA, heightB)
            params = fitResults2(1:4,1);
            uncerts = fitResults2(1:4,2);
        else
            %fprintf('B %f, %f\n', heightA, heightB)
            params = fitResults2(5:8,1);
            uncerts = fitResults2(5:8,2);
        end
        gauss = {Gaussian1D(), Gaussian1D()};
        fs = FunctionSum(gauss);
        xfit = linspace(min(xout), max(xout), 100)';
        fs.setX(xfit);
        yfit = fs.eval(fitResults2(:,1));
    end
else
    % use 1 Gaussian
    params = fitResults1(:,1);
    uncerts = fitResults1(:,2);
    xfit = linspace(min(xout), max(xout), 100)';
    gauss = Gaussian1D();
    gauss.setX(xfit);
    yfit = gauss.eval(params);
end

end
