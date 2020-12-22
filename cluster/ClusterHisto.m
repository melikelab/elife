% Replacement program for the ddcHistoManualRotateXY.m function
%
classdef ClusterHisto < handle
    properties
       colors = ['r' 'b' 'g']; 
    end
    
    methods
        
        function data = getChannelData(~, filename)
        % function data = getChannelData(filename)
        %
        % Returns the localizations in each channel for the
        % specified filename.
        %
            assert(isa(filename, 'char'), 'filename must be a string')
            assert(strcmp(filename(end-3:end),'.csv'), 'Must be a csv file')
            dataIn = dlmread(filename, ',');
            
            % ignore "empty" columns
            nCols = 0;
            for i = 1:size(dataIn, 2)
                indices = find(dataIn(:,i) ~= 0);
                if ~isempty(indices)
                    nCols = nCols + 1;
                end
            end
            assert(nCols == 6 || nCols == 4, 'The number of columns in the file must be 4 or 6')
            
            % separate the data into different channels
            j = 1;
            for i=1:2:nCols
                indices = find(dataIn(:,i) ~= 0);
                data{j} = dataIn(indices, [i i+1]);
                j = j + 1;
            end
            
        end
        
        function rotatedData = rotate(~, data, radians)
        % function rotatedData = rotate(data, radians)
        %
        % Rotate the localizations by the specified angle, in radians.
        % The rotation is performed about the median (x,y) value for all
        % channels that are in data.
        %
            assert(iscell(data), 'Must pass in the data from getChannelData(filename)')
            
            % calculate the common centroid for all localizations
            numChannels = size(data, 2);
            x = []; y = [];
            for i=1:numChannels
                x = [x; data{i}(:,1)];
                y = [y; data{i}(:,2)];
            end
            xc = median(x);
            yc = median(y);
 
            % apply the rotation
            rotatedData = cell(1, numChannels);
            c = cos(radians);
            s = sin(radians);
            for j=1:numChannels
                npts = size(data{j}, 1);
                temp = zeros(npts, 2);
                for i=1:npts
                    x = data{j}(i,1) - xc;
                    y = data{j}(i,2) - yc;
                    temp(i,1) = c*x - s*y + xc;
                    temp(i,2) = s*x + c*y + yc;
                end
                rotatedData{j} = temp;
            end
            
        end
        
        function plotData(self, data, theTitle)
        % function plotData(data, theTitle)
        %
        % Displays a plot of the data with the specified title
        %
            assert(iscell(data), 'Must pass in the data from getChannelData(filename)')
            
            numChannels = size(data, 2);
            figure
            hold on            
            for i = 1:numChannels
                plot(data{i}(:,1), data{i}(:,2), [self.colors(i) '+']);
            end
            if numChannels == 3
                legend('ch1', 'ch2', 'ch3', 'Location','northoutside','Orientation','horizontal');
            else
                legend('ch1', 'ch2', 'Location','northoutside','Orientation','horizontal');
            end
            title(theTitle)
        end
        
        function showHistograms(self, data, theTitle, nBins)
        % function showHistograms(data, theTitle, nBins)
        %
        % Create and display a histogram for all the channels
        %
            assert(iscell(data), 'Must pass in the data from getChannelData(filename)')
            figure
            hold on 
            numChannels = size(data, 2);
            for i=1:numChannels
                y = data{i}(:,2);
                if nargin == 3
                    [n, c] = hist(y, linspace(min(y), max(y), round(sqrt(length(y)))));
                else
                    [n, c] = hist(y, nBins);
                end
                bar(c, n, self.colors(i));
            end
            title(theTitle)
        end
        
        function [result, chisqr, xy, xyFit] = fit(~, singleChannelData, numGaussians, nBins, showPlots, xcenters)
        % function [result, chisqr, xy, xyFit] = fit(singleChannelData, numGaussians, nBins, showPlots, xcenters)
        %
        % Creates a histogram of the singleChannelData using nBins bins and 
        % then fits the histogram data to the specified number of Gaussian 
        % functions.
        % 
        % The result is a Nx2 matrix with the number of rows = 4 * numGaussians
        % [offset uncertainty]
        % [sigma  uncertainty]
        % [area   uncertainty]
        % [mean   uncertainty]
        %  ... repeat the 4 row pattern for the next Gaussians ...
        %
            assert(~iscell(singleChannelData), 'The data must be for a SINGLE channel')
            
            nCols = size(singleChannelData, 1);
            if nCols == 2
                data = singleChannelData(:,2);
            else
                data = singleChannelData(:,1);
            end
            n = numGaussians;
            
            if nargin == 3
                nBins = round(sqrt(length(data)));
                [y, x] = hist(data, linspace(min(data), max(data), nBins));
                showPlots = true;
            elseif nargin == 4
                [y, x] = hist(data , nBins);
                showPlots = true;
            elseif nargin == 6                
                [y, x] = hist(data , xcenters);
            else
                [y, x] = hist(data , nBins);
            end
            
            xy = [x y];
            
            if showPlots
                figure
                bar(x, y, 'r');
                hold on
            end
            
            % create a Gaussian object for the fit for each gaussian
            for i=1:n
                gauss{i} = Gaussian1D();
            end

            % create a FunctionSum object
            fs = FunctionSum(gauss);

            % set the X data for the fit function
            fs.setX(x);
            
            % set an initial guess
            if n == 1
                sigma = std(x);
                area = max(y) * sigma * 2;
                center = mean(x);
                guess = [0 sigma area center]; 
            else
                minDist = floor(length(y)/n) - 1;
                minHeight = 5;
                [pks, indices] = findpeaks(y, 'MINPEAKDISTANCE', minDist, 'MINPEAKHEIGHT', minHeight);
                while length(indices) ~= n
                    minDist = minDist - 1;
                    if minDist == 0
                       fprintf(2, 'Warning in peakfinding. Cannot find %d peaks\n', n);
                       result = ones(n*4, 2)*NaN;
                       chisqr = NaN;
                       return;
                    end
                    [pks, indices] = findpeaks(y, 'MINPEAKDISTANCE', minDist, 'MINPEAKHEIGHT', minHeight);
                end
                guess = [];
                sigma = 25;
                for i=1:length(indices)
                    idx = indices(i);                    
                    area = y(idx) * sigma * 2;
                    center = x(idx);
                    guess = [guess [0 sigma area center]];
                end
            end
            %guess
            fs.setGuess(guess);

            % don't float the background of each Gaussian
            floatP = [];
            for i=1:n
                floatP = [floatP [0 1 1 1]];
            end
            fs.setFloatParams(floatP);

            % create the LevMar fitter
            levmar = LevMar(fs, y);

            % fit
            exit_code = levmar.solve();
            
            txt = sprintf('nBins = %d; centers = ',nBins);

            % check if fit converged
            if exit_code > 0
                params = levmar.bestParams;
                uncerts = levmar.paramUncerts;
                result = [params uncerts];                
                xfit = linspace(min(x), max(x), 100)';
                fs.setX(xfit);
                yfit = fs.eval(params);
                chisqr = levmar.reducedChisqr;
                if showPlots
                    plot(xfit, yfit, 'k', 'linewidth', 3);                    
                    for i=1:numGaussians
                        txt = [txt num2str(params(i*4)) ', '];
                    end
                    title(txt(1:end-2))
                end
                xyFit = [xfit yfit];
            else
                if showPlots
                    title(txt)                    
                end
                fprintf(2, 'Could not fit\n');
                result = ones(n*4, 2)*NaN;
                chisqr = NaN;
                xfit = [];
                yfit = [];
                xyFit = [];
            end
            
        end
        
    end
    
end