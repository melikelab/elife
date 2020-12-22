classdef SynapseFunctions < handle
    
    properties
        homer = '';
        bassoon = '';
        proteins = '';
        synapses = {};
        insight3 = Insight3('junk.bin');
        nm_per_pixel = 0;
        maxNumGaussians = 3;
    end
    
    properties(Hidden)
        ixc = 3; % int, the column index for drift-corrected x localization values
        iyc = 4; % int, the column index for drift-corrected y localization values
        colours = [[0 0 1]; [1 1 0]; [1 0 1]; [0 1 1]; [1 170./255. 0]; [0 169./255. 1]; [148./255. 0 1]; [1 1 1]];
    end
    
    methods
        
        function self = SynapseFunctions(homer, bassoon, proteins, nm_per_pixel)
        % function self = SynapseFunctions()
        %
        % Constructor  
        %
        % Inputs
        % ------
        % homer : string
        %   filename of .bin file from FindClusters (for Homer)
        % bassoon : string
        %   filename of .bin file from FindClusters (for Bassoon)
        % proteins : string
        %   filename of .bin file from FindClusters (for the Proteins)
        % nm_per_pixel : int
        %   the number of nanometers per pixel
        %
            assert(isa(homer,'char'), 'Must pass in a .bin filename for Homer')
            assert(isa(bassoon,'char') || isempty(bassoon), 'Must pass in a .bin filename for Bassoon')
            assert(isa(proteins,'char'), 'Must pass in a .bin filename for the Proteins')
            
            self.homer = homer;
            self.bassoon = bassoon;
            self.proteins = proteins;
            self.nm_per_pixel = nm_per_pixel;
        end
        
        function synapses = findSynapses(self, max_distance, createInsight3)
        % synapses = findSynapses(max_distance)
        %
        % Finds Synapses. A synapse is when there is a Homer and a Bassoon 
        % that are within max_distance of each other.
        %
        % Inputs
        % ------
        % max_distance : int
        %   the maximum distance, in nm, that homer and basson must be to be
        %   considered as a synapse
        % createInsight3 : logical (optional)
        %   create an Insight3 molecule list
        %
        % Returns
        % -------
        % synapses : cell array of Synapse objects
        %
            if nargin == 2
                createInsight3 = true;
            end
        
            i3a = Insight3(self.homer);
            if ~i3a.isBinFormat()
                i3a.convertFormat()
            end
            
            if isempty(self.bassoon)
                i3b = Insight3();
                i3b.setData(i3a.getData());
            else
                i3b = Insight3(self.bassoon);
                if ~i3b.isBinFormat()
                    i3b.convertFormat()
                end
            end
            
            [i3a_ch0, i3a_idx] = i3a.getChannels(0);
            i3a_idx = [i3a_idx; i3a.numMolecules+1];
            [i3b_ch0, i3b_idx] = i3b.getChannels(0);
            i3b_idx = [i3b_idx; i3b.numMolecules+1];
            
            nrows = size(i3a_ch0,1);
            synapses = cell(nrows, 1);
            cnt = 1;
            associatedBassoon = [];
            for ia = 1:nrows
                
                % find the closest neighbour
                minDist = Inf;
                for j = 1:size(i3b_ch0,1)
                    dist = self.nm_per_pixel * sqrt( (i3a_ch0(ia,self.ixc) - i3b_ch0(j,self.ixc))^2 ...
                                                   + (i3a_ch0(ia,self.iyc) - i3b_ch0(j,self.iyc))^2);
                    if dist < minDist
                        minDist = dist;
                        ib = j;
                    end
                end
                
                if minDist < max_distance
                    synapses{cnt} = Synapse();
                    
                    h = i3a.data( i3a_idx(ia)+1:i3a_idx(ia+1)-1 ,:);
                    h(:,i3a.getColumnIndex('channel')) = 1;
                    synapses{cnt}.homer = h;
                    synapses{cnt}.homerRotated = h;
                    sigmax = std(h(:,self.ixc)) * self.nm_per_pixel;
                    sigmay = std(h(:,self.iyc)) * self.nm_per_pixel;
                    synapses{cnt}.sigmaHomer = max(sigmax,sigmay);
                    synapses{cnt}.sigmaHomerDescription = 'equal to max(sigma_x, sigma_y) [in nanometers]';
                    synapses{cnt}.areaUnits = 'nm';                    
                    synapses{cnt}.areaHomer = pi * sigmax * sigmay;
                    
                    b = i3b.data( i3b_idx(ib)+1:i3b_idx(ib+1)-1 ,:);
                    b(:,i3b.getColumnIndex('channel')) = 2;
                    synapses{cnt}.bassoon = b;
                    synapses{cnt}.bassoonRotated = b;
                    sigmax = std(b(:,self.ixc)) * self.nm_per_pixel;
                    sigmay = std(b(:,self.iyc)) * self.nm_per_pixel;
                    synapses{cnt}.sigmaBassoon = max(sigmax,sigmay);
                    synapses{cnt}.sigmaBassoonDescription = 'equal to max(sigma_x, sigma_y) [in nanometers]';
                    synapses{cnt}.areaBassoon = pi * sigmax * sigmay;
                    
                    xcentroid = median([h(:,self.ixc); b(:,self.ixc)]);
                    ycentroid = median([h(:,self.iyc); b(:,self.iyc)]);
                    synapses{cnt}.origin = [xcentroid ycentroid];
                    synapses{cnt}.originDescription = '[median(homer_xc,bassoon_xc) median(homer_yc,bassoon_yc)]';
                    
                    self.insight3.appendData(h);
                    self.insight3.appendData(b);
                    
                    cnt = cnt + 1;
                    
                    % check if this Bassoon is already associated with a Homer
                    if ~isempty( find(associatedBassoon == ib) )
                        %fprintf(2, 'WARNING! Bassoon at (x,y)=(%.1f,%.1f) is already associated with a Homer\n', i3b_ch0(ib,self.ixc), i3b_ch0(ib,self.iyc));
                        error('Bassoon at (x,y)=(%.1f,%.1f) is already associated with a Homer', i3b_ch0(ib,self.ixc), i3b_ch0(ib,self.iyc));
                    end
                    associatedBassoon = [associatedBassoon; ib];
                    
                end
                
            end
            % delete the empty rows from the synapses cell array
            synapses( all(cellfun(@isempty,synapses),2), : ) = [];
            
            % create the homer_bassoon molecule list
            self.insight3.forceFileOverwrite(true);
            fidx = regexp(self.homer, '_ch\d+_roi\d+_.+_min\d+');
            if isempty(fidx)
                fidx = regexp(self.homer, '_roi\d+_th\d+_fac');
            end
            if createInsight3
                self.insight3.write([self.homer(1:fidx) 'homer_bassoon.bin']);
            end
            fprintf('Created %d synapses\n', length(synapses));
            
            self.synapses = synapses;
            
        end
        
        function synapses = addProteins(self, max_distance, createInsight3)
        % synapses = addProteins(max_distance)
        %
        % Adds the proteins to the appropriate synapse.
        %
        % Inputs
        % ------
        % max_distance : int
        %   the maximum distance, in nm, that a protein can be from the
        %   Synapse "origin" value and be tagged as belonging to that
        %   synapse
        % createInsight3 : logical (optional)
        %   create an Insight3 molecule list
        %
        % Returns
        % -------
        % synapses : cell array of Synapse objects with the proteins added
        %

            if nargin == 2
                createInsight3 = true;
            end

            if isempty(self.synapses)
                error('There are no synapses');
            end
            
            i3 = Insight3(self.proteins);
            
            if ~i3.isBinFormat()
                i3.convertFormat()
            end
            
            [i3_ch0, i3_idx] = i3.getChannels(0);
            i3_idx = [i3_idx; i3.numMolecules+1];
            
            xc = i3.getColumnIndex('xc');
            yc = i3.getColumnIndex('yc');
            
            % don't put these protein in ch1 or ch2 in the homer_bassoon molecule list
            channels = [3 4 5 6 7 8 9 0];
            
            synapseChannels = ones(length(self.synapses), 1);
            
            channelIndex = i3.getColumnIndex('channel');
            
            % for each protein cluster
            for i = 1:size(i3_ch0,1)
                
                % find the closest synapse
                minDist = Inf;
                for j = 1:length(self.synapses)
                    xy = self.synapses{j}.origin;
                    dist = self.nm_per_pixel * sqrt( (i3_ch0(i,xc) - xy(1))^2 ...
                                                   + (i3_ch0(i,yc) - xy(2))^2);
                    if dist < minDist
                        minDist = dist;
                        idx = j;
                    end
                end
                
                % add the protein to the closest synapse
                if minDist < max_distance
                    nProteins = length(self.synapses{idx}.proteins);
                    
                    protein = i3.data( i3_idx(i)+1:i3_idx(i+1)-1 ,:);
                    self.synapses{idx}.proteins{nProteins+1,1} = protein;
                    self.synapses{idx}.proteinsRotated{nProteins+1,1} = protein;                    
                    sigmax = std(protein(:,self.ixc)) * self.nm_per_pixel;
                    sigmay = std(protein(:,self.iyc)) * self.nm_per_pixel;
                    self.synapses{idx}.sigmaXProteins{nProteins+1,1} = sigmax;
                    self.synapses{idx}.sigmaYProteins{nProteins+1,1} = sigmay;                    
                    self.synapses{idx}.areaProteins{nProteins+1,1} = pi * sigmax * sigmay;
                    
                    protein(:,channelIndex) = channels(synapseChannels(idx));
                    self.insight3.appendData(protein);
                    
                    synapseChannels(idx) = synapseChannels(idx) + 1;
                    if synapseChannels(idx) > length(channels)
                        synapseChannels(idx) = 1;
                    end
                    
                end
                
            end
            
            if createInsight3
                self.insight3.write(strrep(self.insight3.filename, '.bin', '_proteins.bin'));
            end
            
            synapses = self.synapses;
            
        end
        
        function synapses = autoAngle(self, index)
        % function synapses = autoAngle(index)
        %
        % Find the angle that makes the slope of the line that is fit to
        % the (xc,yc) values of the localizations for Homer to be 0
        %
        % Uses the "origin" value in the Synapse object as a reference to 
        % perform the rotation about.
        %
        % If you call this function with no input arguments then the angle
        % will be calculated for all synapses.
        %
        % Inputs
        % ------
        % index (optinal argument) : int
        %   The index of the synapse
        %
            if isempty(self.synapses)
                error('There are no synapses');
            end
            
            if nargin > 1
                assert(index > 0 && index <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', index, length(self.synapses));
                istart = index;
                iend = index;
            else
                istart = 1;
                iend = length(self.synapses);
            end
            
            for i = istart:iend
                data = self.synapses{i}.homer(:,3:4);
                angle = self.makeSlopeZero(data, self.synapses{i}.origin, i);
                self.setAngle(i, angle);
            end
            
            synapses = self.synapses;
            
        end
        
        function synapses = setAngle(self, index, angle)
        % function angle = setAngle(index, angle)
        %
        % Rotates the localizations in the specified synapse by the 
        % requested angle.
        %
        % Uses the "origin" value in the Synapse object as a reference to 
        % perform the rotation about.
        % 
        % Inputs
        % ------
        % index : int
        %   the synapse index in the cell array
        % angle : double
        %   the abgle to rotate the localizations in the synapse
        %
            if isempty(self.synapses)
                error('There are no synapses');
            else
                assert(index > 0 && index <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', index, length(self.synapses));
            end

            self.synapses{index}.angle = angle;
            
            origin = self.synapses{index}.origin;
            
            rotated = self.rotate(self.synapses{index}.homer(:,1:2), angle, origin);
            self.synapses{index}.homerRotated(:,1:2) = rotated;

            rotated = self.rotate(self.synapses{index}.homer(:,3:4), angle, origin);
            self.synapses{index}.homerRotated(:,3:4) = rotated;

            rotated = self.rotate(self.synapses{index}.bassoon(:,1:2), angle, origin);
            self.synapses{index}.bassoonRotated(:,1:2) = rotated;

            rotated = self.rotate(self.synapses{index}.bassoon(:,3:4), angle, origin);
            self.synapses{index}.bassoonRotated(:,3:4) = rotated;
            
            for i = 1:length(self.synapses{index}.proteins)
                rotated = self.rotate(self.synapses{index}.proteins{i}(:,1:2), angle, origin);
                self.synapses{index}.proteinsRotated{i}(:,1:2) = rotated;

                rotated = self.rotate(self.synapses{index}.proteins{i}(:,3:4), angle, origin);
                self.synapses{index}.proteinsRotated{i}(:,3:4) = rotated;
            end
            
            synapses = self.synapses;
            
        end
                
        function synapses = fitHomer(self, varargin)
        % function synapses = fitHomer(varargin)
        %
        % Bin the rotated Homer data and then fit to Gaussian functions.
        %
        % Inputs
        % ------
        % index (optional argument) : int
        %   The synapse index in the cell array.
        %   The default is to fit the Homer data in every synapse.
        % numGaussians (optional argument) : int
        %   The number of Gaussian functions to use to fit the data.
        %   The default is to fit for numGaussians = 1:self.maxNumGaussians 
        %   and to use the fit with the smallest chi-square value
        % numBins (optional argument) : int
        %   The number of bins to use for the histogram
        %   The default is to use round(sqrt(length(data)))
        %
            self.fitLoop('homer', varargin{:});
            synapses = self.synapses;
        end
        
        function synapses = fitBassoon(self, varargin)
        % function synapses = fitBassoon(varargin)
        %
        % Bin the rotated Bassoon data and then fit to Gaussian functions.
        %
        % Inputs
        % ------
        % index (optional argument) : int
        %   The synapse index in the cell array.
        %   The default is to fit the Bassoon data in every synapse.
        % numGaussians (optional argument) : int
        %   The number of Gaussian functions to use to fit the data.
        %   The default is to fit for numGaussians = 1:self.maxNumGaussians 
        %   and to use the fit with the smallest chi-square value
        % numBins (optional argument) : int
        %   The number of bins to use for the histogram
        %   The default is to use round(sqrt(length(data)))
        %
            self.fitLoop('bassoon', varargin{:});
            synapses = self.synapses;
        end

        function synapses = fitProtein(self, iSynapse, iProtein, varargin)
        % function synapses = fitProtein(iSynapse, iProtein)
        %
        % Fit the rotated Protein in this synapse data to a Gaussian
        %
        % If no inputs arguments are passed in the fits all Protein data
        % in all synpases to determine the Gaussian fit.
        %
        % If iProtein is not specified then fits all protein data in the
        % specified synapse to determine the Gaussian fit.
        %
        % Inputs
        % ------
        % iSynapse (optional argument) : int
        %   the synapse index in the cell array
        %
        % iProtein (optional argument) : int
        %   the protein index in the cell array
        %
        % varargin Param-Value options
        %   numGaussians (optional argument) : int
        %     The number of Gaussian functions to use to fit the data.
        %     The default is to fit for numGaussians = 1:self.maxNumGaussians 
        %     and to use the fit with the smallest chi-square value
        %   numBins (optional argument) : int
        %     The number of bins to use for the histogram
        %     The default is to use round(sqrt(length(data)))
        %   
            if isempty(self.synapses)
                error('There are no synapses');
            end
             
            if nargin == 1
                index = -1;
                proteinIndex = -1;
            elseif nargin == 2
                assert(iSynapse > 0 && iSynapse <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', iSynapse, length(self.synapses));
                index = iSynapse;
                proteinIndex = -1;
            else
                assert(iSynapse > 0 && iSynapse <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', iSynapse, length(self.synapses));
                index = iSynapse;
                assert(iProtein > 0 && iProtein <= length(self.synapses{iSynapse}.proteins), 'Invalid protein index of %d. Valid range is: 1 <= index <= %d', iProtein, length(self.synapses{iSynapse}.proteins));
                proteinIndex = iProtein;
            end
            
            % insert the protein index at the beginning of varargin
            varargin(2:end+1) = varargin(1:end);
            varargin(1) = {proteinIndex};
            varargin(2:end+1) = varargin(1:end);
            varargin(1) = {'proteinIndex'};

            % insert the synapse 'index' at the beginning of varargin
            varargin(2:end+1) = varargin(1:end);
            varargin(1) = {index};
            
            self.fitLoop('protein', varargin{:});
            synapses = self.synapses;
        end

        function synapses = calculateDistances(self, index)
        % function distances = calculateDistances(index)
        %
        % Calculates the distances between Homer, Bassoon and the
        % Protein(s) for the specified synapse.
        %
        % Inputs
        % ------
        % index : int
        %   the cell index of the synapse
        %
            if isempty(self.synapses)
                error('There are no synapses');
            else
                assert(index > 0 && index <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', index, length(self.synapses));
            end            

            synapse = self.synapses{index};
            
            if isempty(synapse.homerRotatedResult)
                error('Cannot calculate distances for synapse %d. Homer has not been fit. Run fitHomer()\n', index);
            end
            [centerH,~] = self.findBestPeak(synapse.homerRotatedResult);
            
            if isempty(synapse.bassoonRotatedResult)
                error('Cannot calculate distances for synapse %d. Bassoon has not been fit. Run fitBassoon()\n', index);
            end
            [centerB,~] = self.findBestPeak(synapse.bassoonRotatedResult);
            
            if ~isempty(synapse.proteins) && isempty(synapse.proteinsRotatedResult)
                error('Cannot calculate protein-homer distances for synapse %d. Run fitProtein()\n', index);
            end
            
            self.synapses{index}.distances = {};
            self.synapses{index}.distances.units = 'nm';
            self.synapses{index}.distances.bassoon_homer = self.nm_per_pixel * (centerB - centerH);
            if isempty(synapse.proteins)
                self.synapses{index}.distances.proteins_homer = NaN;
                self.synapses{index}.distances.isProteinBetween = false;
            else
                for i = 1:length(synapse.proteinsRotatedResult)
                    [centerP,~] = self.findBestPeak(synapse.proteinsRotatedResult{i});                
                    self.synapses{index}.distances.proteins_homer(i) = self.nm_per_pixel * (centerP - centerH);
                    if centerB > centerH
                        self.synapses{index}.distances.isProteinBetween(i) = centerP < centerB && centerP > centerH; 
                    else
                        self.synapses{index}.distances.isProteinBetween(i) = centerP > centerB && centerP < centerH;
                    end
                end
            end
            
            synapses = self.synapses;
            
        end

        function setSynapses(self, synapses)
        % function setSynapses(synapses)
        %
        % Set the synapse cell array
        %
            assert(~isempty(synapses), 'Empty synapses');
            assert(isa(synapses,'cell'), 'Invalid synapses. Must be a cell array of Synapse objects')
            for i = 1:length(synapses)
                assert(isa(synapses{i},'Synapse'), 'Invalid synapses. Must be a cell array of Synapse objects')
            end            
            self.synapses = synapses;
            fprintf('Updated. Created %d synapses\n', length(synapses));
        end

        function plotIDs(self)
        % function plotIDs()
        %
        % Plots all synapses with the ID number written on top of each
        % synapse
        %
            if isempty(self.insight3.data)
                error('There are no synapses');
            end
            self.insight3.show();
            for i = 1:length(self.synapses)
                xy = self.synapses{i}.origin;
                text(xy(1),xy(2),sprintf('%d', i),'Color','red','FontSize',11,'FontWeight','bold');
            end
        end
        
        function plot(self, index)
        % function plot(index)
        % 
        % Displays the localizations in the specified synpase index and, 
        % if available, the histogram
        %
            if isempty(self.synapses)
                error('There are no synapses');
            else
                assert(index > 0 && index <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', index, length(self.synapses));
            end
            
            figure('Name',sprintf('Synapse %d', index), 'NumberTitle', 'off');
            synapse = self.synapses{index};
            self.plotLocs(1, synapse, synapse.homer, synapse.bassoon, synapse.proteins)
            self.plotLocs(2, synapse, synapse.homerRotated, synapse.bassoonRotated, synapse.proteinsRotated)
            
            if isempty(synapse.homerRotatedX) && isempty(synapse.bassoonRotatedX) && isempty(synapse.proteinsRotatedX)
                return
            end
            
            subplot(1,3,3);
            hold on
            
            legendLabels = {};
            if ~isempty(synapse.homerRotatedX)
                legendLabels{1} = 'homer';
                bar(synapse.homerRotatedX, synapse.homerRotatedY, 'g');
            end
            if ~isempty(synapse.bassoonRotatedX)
                legendLabels{length(legendLabels)+1} = 'bassoon';
                bar(synapse.bassoonRotatedX, synapse.bassoonRotatedY, 'r');
            end
            if ~isempty(synapse.proteinsRotatedX)
                colors = self.colours;
                if length(synapse.proteinsRotatedX) > length(colors)
                    colors = colormap;
                end
                for i = 1:length(synapse.proteinsRotatedX)
                    legendLabels{length(legendLabels)+1} = sprintf('protein%d', i);
                    bar(synapse.proteinsRotatedX{i}, synapse.proteinsRotatedY{i}, 'EdgeColor', colors(i,:), 'FaceColor', colors(i,:));
                end
            end
            legend(legendLabels, 'Location', 'northoutside');

            if ~isempty(synapse.homerRotatedXFit)
                plot(synapse.homerRotatedXFit, synapse.homerRotatedYFit, 'k', 'linewidth', 3);
            end
            if ~isempty(synapse.bassoonRotatedXFit)
                plot(synapse.bassoonRotatedXFit, synapse.bassoonRotatedYFit, 'k', 'linewidth', 3);
            end
            if ~isempty(synapse.proteinsRotatedXFit)
                for i = 1:length(synapse.proteinsRotatedXFit)
                    plot(synapse.proteinsRotatedXFit{i}, synapse.proteinsRotatedYFit{i}, 'k', 'linewidth', 3);
                end
            end

            title(sprintf('Synapse %d', index));
        end

        function print(self)
        % function print()
        %
        % Prints a summary of the number of proteins in each synapse to the
        % Command Window.
        %
            fprintf('\n  Synapses with proteins\n=========================\n');
            j = 1;
            temp = {};
            for i = 1:length(self.synapses)
                if ~isempty(self.synapses{i}.proteins)
                    n = length(self.synapses{i}.proteins);
                    if n > 1
                        temp{j,1} = sprintf('%2d proteins in Synapse %2d\n', n, i);
                    else
                        temp{j,1} = sprintf('%2d protein  in Synapse %2d\n', n, i);
                    end
                    j = j + 1;
                end
            end
            temp = sort(temp);
            fprintf('%s', temp{:});
            fprintf('\nSynapses without proteins\n=========================\n');
            for i = 1:length(self.synapses)
                if isempty(self.synapses{i}.proteins)
                    fprintf('No proteins in Synapse %2d\n', i);
                end
            end
            fprintf('\n');
        end
        
    end
    
    methods(Hidden)
        
        function theta = makeSlopeZero(~, data, origin, index)
        % function theta = makeSlopeZero(data, origin, index)
        %
        % Rotate the points in data so that a linear fit to the 
        % points has a slope = 0. The rotation is performed about 
        % the "origin" coordinate.
        %
        % Inputs
        % ------
        % data : 2D array
        %   The data to rotate
        % origin : 1D array
        %   The origin to use for the rotation
        % index : int
        %   The synapse index
        %
        % Returns
        % -------
        % theta : double
        %   The angle that was used to rotate the data to make the slope 0
        %
            
            % make sure we don't get stuck in an infinite loop
            % typically, < 10 iterations are needed to make the slope 0
            maxCount = 1000;
            
            rotated = data;

            % determine the rotation angle that makes the slope zero
            slope = Inf; theta = 0.0; cnt = 0;
            while abs(slope) > 1e-8 && cnt < maxCount
                
                % fit the (x,y) points to a linear function to determine the slope
                fitvars = polyfit(rotated(:,1), rotated(:,2), 1);
                slope = fitvars(1);
                angle = -atan(slope);
                
                % rotate the points
                c = cos(angle);
                s = sin(angle);
                for m = 1:size(rotated, 1)
                    x = rotated(m,1) - origin(1);
                    y = rotated(m,2) - origin(2);
                    xp = c*x - s*y;
                    yp = s*x + c*y;
                    rotated(m,1) = xp + origin(1);
                    rotated(m,2) = yp + origin(2);
                end
                
                % increment the rotation angle
                theta = theta + angle;
                
                % increment the counter
                cnt = cnt + 1;
            end
            
            if cnt >= maxCount
                fprintf(2, 'Maximum rotation iterations reached for synapse %d\n', index);
            end                
            %assert(cnt < maxCount, 'Maximum rotation iterations reached for synapse %d', index)
        end
        
        function rotated = rotate(~, data, angle, origin)
        % function rotated = rotate(data, angle, origin)
        %
        % Rotate data by the specified angle. The rotation is performed 
        % about the "origin" coordinate.
        %
        % Inputs
        % ------
        % data : 2D array
        %   The data to rotate
        % angle : double
        %   The angle to use to rotate the data
        % origin : 1D array
        %   The x,y coordinates of the origin to use to rotate the data
        %
        % Returns
        % -------
        % rotated : 2D array
        %   The rotated data
        %
            c = cos(angle);
            s = sin(angle);
            rotated = data;
            for m = 1:size(rotated, 1)
                x = rotated(m,1) - origin(1);
                y = rotated(m,2) - origin(2);
                xp = c*x - s*y;
                yp = s*x + c*y;
                rotated(m,1) = xp + origin(1);
                rotated(m,2) = yp + origin(2);
            end
        end
        
        function [result, x, y, xfit, yfit, chisq] = fit(~, data, numGaussians, nBins, flag, indexSynapse, indexProtein, nGaussians)
        % function [result, x, y, xfit, yfit, chisq] = fit(values, numGaussians, nBins)
        %
        % Creates a histogram of the data using nBins bins and then fits 
        % the histogram data to the specified number of Gaussian functions.
        % 
        % The result is a Nx2 matrix with the number of rows = 4 * numGaussians
        % [offset uncertainty]
        % [sigma  uncertainty]
        % [area   uncertainty]
        % [mean   uncertainty]
        %  ... repeat the 4 row pattern for the next Gaussians ...
        %
            if strcmp(flag, 'protein')
                fprintf('Fitting Protein %d in Synapse %d to %d Gaussians\n', indexProtein, indexSynapse, nGaussians)
            else
                fprintf('Fitting %s in Synapse %d to %d Gaussians\n', flag, indexSynapse, nGaussians)
            end
        
            [y, x] = hist(data , nBins);
            
            result = []; xfit = []; yfit = []; chisq = NaN;
            
            % create a Gaussian object for the fit for each gaussian
            for i=1:numGaussians
                gauss{i} = Gaussian1D();
            end

            % create a FunctionSum object
            fs = FunctionSum(gauss);

            % set the X data for the fit function
            fs.setX(x);
            
            % set an initial guess
            if numGaussians == 1
                sigma = std(x);
                area = max(y) * sigma * 2 * 1.5;
                center = mean(x);
                guess = [0 sigma area center]; 
            else
                minDist = floor(length(y)/numGaussians) - 1;
                [~, indices] = findpeaks(y, 'MINPEAKDISTANCE', minDist);
                while length(indices) ~= numGaussians
                    minDist = minDist - 1;
                    if minDist == 0
                       fprintf(2, 'Cannot find %d peaks\n', numGaussians);
                       return;
                    end
                    [~, indices] = findpeaks(y, 'MINPEAKDISTANCE', minDist);
                end
                guess = [];
                sigma = 0.5;
                for i=1:length(indices)
                    idx = indices(i);                    
                    area = y(idx) * sigma * 2 * 1.5;
                    center = x(idx);
                    guess = [guess [0 sigma area center]];
                end
            end
            fs.setGuess(guess);

            % don't float the background of each Gaussian
            floatP = [];
            for i=1:numGaussians
                floatP = [floatP [0 1 1 1]];
            end
            fs.setFloatParams(floatP);

            % create the LevMar fitter
            levmar = LevMar(fs, y);

            % fit, and check if fit converged
            if levmar.solve() > 0
                if any(isnan(levmar.paramUncerts))
                    return
                end
                xfit = linspace(min(x), max(x), 100)';
                fs.setX(xfit);
                yfit = fs.eval(levmar.bestParams);
                chisq = levmar.reducedChisqr;
                result = [levmar.bestParams levmar.paramUncerts];
            end
            
        end
        
        function fitLoop(self, flag, varargin)
        % function fitLoop(flag, varargin)    
        %
        % flag : string
        %   either 'homer', 'bassoon' or 'protein'
        % varargin = varargin{:} from the calling function
        %
            if isempty(self.synapses)
                error('There are no synapses');
            end
            
            p = inputParser;
            
            defaultIndex = -1; % do all synapses
            checkIndex = @(x) x >= -1 && x <= length(self.synapses);
            
            defaultNumGaussians = -1; % fit for numGaussians = 1:self.maxNumGaussians and use the best fit
            checkNumGaussians = @(x) x > 0 && x <= self.maxNumGaussians;
            
            defaultNumBins = -1; % determine based on the number of data points
            
            defaultProteinIndex = -1;
            
            addRequired(p,'flag',@ischar);
            addOptional(p,'index',defaultIndex,checkIndex);
            addParamValue(p,'numGaussians',defaultNumGaussians,checkNumGaussians)
            addParamValue(p,'numBins',defaultNumBins,@isnumeric)
            addParamValue(p,'proteinIndex',defaultProteinIndex,@isnumeric)
            
            parse(p,flag,varargin{:})
            
            if p.Results.index == -1
                istart = 1;
                iend = length(self.synapses);
            else
                istart = p.Results.index;
                iend = p.Results.index;
                assert(iend > 0 && iend <= length(self.synapses), 'Invalid synapse index of %d. Valid range is: 1 <= index <= %d', iend, length(self.synapses));
            end
            
            if p.Results.numGaussians == -1
                iGauss1 = 1;
                iGauss2 = self.maxNumGaussians;
            else
                iGauss1 = p.Results.numGaussians;
                iGauss2 = p.Results.numGaussians;
            end

            for i = istart:iend
                
                if strcmp(p.Results.flag, 'protein')
                    if p.Results.proteinIndex == -1
                        istartP = 1;
                        istopP = length(self.synapses{i}.proteinsRotated);
                    else
                        istartP = p.Results.proteinIndex;
                        istopP = p.Results.proteinIndex;
                    end
                else
                    % then the protein loop does not make sense since we are
                    % only gitting homer or bassoon. do 1 'dummy' loop
                    istartP = 1;
                    istopP = 1;
                end
                
                for ip = istartP:istopP
                    
                    if strcmp(p.Results.flag, 'homer')
                        data = self.synapses{i}.homerRotated(:,4);
                    elseif strcmp(p.Results.flag, 'bassoon')
                        data = self.synapses{i}.bassoonRotated(:,4);
                    else
                        data = self.synapses{i}.proteinsRotated{ip}(:,4);
                    end

                    if p.Results.numBins == -1
                        nBins = round(sqrt(length(data)));
                    else
                        nBins = p.Results.numBins;
                    end

                    bestResult = []; bestX = []; bestY = []; bestXfit = []; bestYfit = [];
                    minChisq = Inf;
                    for iGauss = iGauss1:iGauss2
                        try
                            [result, x, y, xfit, yfit, chisq] = self.fit(data, iGauss, nBins, p.Results.flag, i, ip, iGauss);
                        catch
                            chisq = NaN;
                        end
                        if ~isnan(chisq) && chisq < minChisq
                            bestResult = result;
                            bestX = x';
                            bestY = y';
                            bestXfit = xfit;
                            bestYfit = yfit;
                            minChisq = chisq;
                        end
                    end

                    if strcmp(p.Results.flag, 'homer')
                        self.synapses{i}.homerRotatedResult = bestResult;
                        [~,widthH] = self.findBestPeak(bestResult);
                        self.synapses{i}.widthHomer = widthH;
                        if isempty(bestX)
                            self.synapses{i}.homerRotatedX = x';
                            self.synapses{i}.homerRotatedY = y';
                        else
                            self.synapses{i}.homerRotatedX = bestX;
                            self.synapses{i}.homerRotatedY = bestY;
                        end
                        self.synapses{i}.homerRotatedXFit = bestXfit;
                        self.synapses{i}.homerRotatedYFit = bestYfit;
                    elseif strcmp(p.Results.flag, 'bassoon')
                        self.synapses{i}.bassoonRotatedResult = bestResult;
                        [~,widthB] = self.findBestPeak(bestResult);
                        self.synapses{i}.widthBassoon = widthB;
                        if isempty(bestX)
                            self.synapses{i}.bassoonRotatedX = x';
                            self.synapses{i}.bassoonRotatedY = y';
                        else                            
                            self.synapses{i}.bassoonRotatedX = bestX;
                            self.synapses{i}.bassoonRotatedY = bestY;
                        end
                        self.synapses{i}.bassoonRotatedXFit = bestXfit;
                        self.synapses{i}.bassoonRotatedYFit = bestYfit;
                    else
                        self.synapses{i}.proteinsRotatedResult{ip,1} = bestResult;
                        if isempty(bestX)
                            self.synapses{i}.proteinsRotatedX{ip,1} = x';
                            self.synapses{i}.proteinsRotatedY{ip,1} = y';
                        else
                            self.synapses{i}.proteinsRotatedX{ip,1} = bestX;
                            self.synapses{i}.proteinsRotatedY{ip,1} = bestY;
                        end
                        self.synapses{i}.proteinsRotatedXFit{ip,1} = bestXfit;
                        self.synapses{i}.proteinsRotatedYFit{ip,1} = bestYfit;                        
                    end

                    if isinf(minChisq)
                        if strcmp(p.Results.flag, 'protein')
                            fprintf(2, 'Cannot fit Protein %d in Synapse %d\n', ip, i);                            
                        else
                            fprintf(2, 'Cannot fit %s in Synapse %d\n', p.Results.flag, i);
                        end
                    end
                end   
            end
            
        end
       
        function [center, width] = findBestPeak(self, data)
        % function [center, width] = findBestPeak(data)
        %
        % Find the Gaussian with the largest peak height and return the
        % corresponding Gaussian center value, in pixels, and width value, in nm.
        %
            if isempty(data)
                center = NaN;
                width = NaN;
            else
                areas = data(3:4:end,1);
                sigmas = data(2:4:end,1);
                heights = areas ./ (sqrt(2*pi) * sigmas );
                [~, imax] = max(heights);
                center = data(imax*4,1);
                width = 2 * sqrt(2*log(2)) * data(imax*2,1) * self.nm_per_pixel;
            end
        end
        
        function plotLocs(self, plotIndex, synapse, h, b, p)
            
            colors = self.colours;
            if length(p) > length(colors)
                colors = colormap;
            end
            if isempty(synapse.homerRotatedX) && isempty(synapse.bassoonRotatedX) && isempty(synapse.proteinsRotatedX)
                subplot(1,2,plotIndex);
            else
                subplot(1,3,plotIndex);
            end
            
            hold on;
            
            if plotIndex == 2
                legendLabels = {'homerRotated', 'bassoonRotated'};
            else
                legendLabels = {'homer', 'bassoon'};
            end
            scatter(h(:,self.ixc), h(:,self.iyc), 'g+');%''green', 'filled');
            scatter(b(:,self.ixc), b(:,self.iyc), 'r+');%'red', 'filled');
            for i = 1:length(p)
                if plotIndex == 2
                    legendLabels{2+i} = sprintf('protein%dRotated', i);
                else
                    legendLabels{2+i} = sprintf('protein%d', i);
                end
                scatter(p{i}(:,self.ixc), p{i}(:,self.iyc), [], colors(i,:), '+');%'filled');
            end
            legend(legendLabels, 'Location', 'southoutside');
            set(gca,'xaxislocation','top','yaxislocation','left','ydir','reverse','Color',[1 1 1]);
            axis('image');
            %axis([0 self.imageWidth 0 self.imageHeight]);
            %set(gca,'xtick',[],'ytick',[])
        end
        
    end
    
end