% The default struct to use for the FindClusters function
%
% i3file : string or an Insight3 object
%    if 'string' then the file path of the Insight3 molecule list. 
%    Can be a '*.bin' or a '*.txt' file
%
% useCorr : boolean
%     Use the drift corrected X,Y localization values to find clusters
%
% useChannel : integer, or an array of ints
%     Specify which channels (categories) in the molecule list you want to use               
%     If useChannel = -1 then use all channels
%     If useChannel = 1 then find clusters for molecules only in channel 1
%     If useChannel = [1 3] then find clusters only using molecules in channels 1 and 3
%     If useChannel = 0:8 then use all channels except for channel 9
%
% width : integer
%     The frame width, in pixels, of the original image file that was used to 
%     generate the molecule list
%
% height : integer
%     The frame height, in pixels, of the original image file that was used to 
%     generate the molecule list
%
% original_pixel_size : double
%     The pixel size, in nm, of the original image
%
% analysis_pixel_size : double
%     The pixel size, in nm, to use to perform the cluster analysis
%     A 2D histogram of the number of localizations per pixel is created
%
% roi : integer, should be an odd-valued integer (but it doesn't have to be)
%     Each pixel contains a certain number of molecules (a 2D histogram with each 
%     pixel size equal to analysis_pixel_size in nm). 
%     Sum the number of molecules in a 'roi' x 'roi' area surrounding each pixel
%     and put the sum in the pixel under consideration
%     Example roi=5,  0  0  2  1  0  the central pixel gets the sum ->  x  x  x  x  x 
%                     1  1  2  2  2                                     x  x  x  x  x  
%                     1  2  4  2  0                                     x  x 31  x  x 
%                     0  2  0  0  1                                     x  x  x  x  x                              
%                     2  3  3  0  0                                     x  x  x  x  x
%
% threshold : integer
%     After applying the 'roi' x 'roi' sum, each pixel must have a value greater
%     than 'threshold' in order to be further considered in the cluster analysis
%     i.e., if 'analysis_pixel_size'=10nm and 'roi'=5 then there must be more than 
%     'threshold' molecules in a 50nm x 50nm ROI surrounding each pixel
%
% factor : double
%     The factor to further reduce the analysis_pixel_size once you have a region 
%     to start finding clusters in 
%     i.e., if 'factor'=5 and 'analysis_pixel_size'=10nm then the program
%     looks for clusters with each pixel size being 2nm
%
% precision : double
%     The precision of a localization, in nm. This value is used as the sigma value 
%     in a 2D Gaussian point-spread function
%
% minCluster : integer
%     After the program finds a cluster a final check is performed to ensure that 
%     there must be >= 'minCluster' molecules within the cluster. If there
%     are < 'minCluster' molecules within the cluster then this is not a cluster.
%
% ignoreNumPeakThreshold : integer
%     When doing the peak finding (finding localizations on an island)
%     if the number of pixels to scan gets too large then the program will
%     take a long time to finsih. 
%     Ignore islands that have > ignoreNumPeakThreshold of pixels.
%        
% ignoreNumIslandThreshold : integer
%     If the number of localizations in a particular island gets larger than
%     this value then the program will take a long time to compute a Sum of
%     Gaussians matrix which is used to find the cluster centroid position.
%     Ignore this island and continue to analyse the next island if the
%     number of localizations exceeds ignoreNumIslandThreshold. 
%
% drawUsingZRange : boolean
%     if true then the localizations in each centroid are set to be in
%     an Insight3 channel that depends on the z value of the centroid. 
%     For example, 9 Insight3 channels are used for 9 z regions. The z 
%     regions are defined as minZ:dz:maxZ, where dz = (maxZ-minZ)/9 and 
%     depending on the z value of the centroid the localizations for 
%     this centroid will go into a particular Insight3 channel.
%
%     zRange index | Insight3 
%     minZ:dz:maxZ | channel 
%                1   5 Magenta              
%                2   9 Violet
%                3   3 Blue 
%                4   8 Cyan-Blue
%                5   6 Cyan
%                6   2 Green
%                7   4 Yellow
%                8   7 Orange
%                9   1 Red
%
%     if false then all localizations in a cluster get assigned to an
%     Insight3 channel and the channel is randomly selected for each
%     cluster to be between 1 and 9.
% 
% use_iterative_segmentation : boolean
%     If an island is too large then use an iterative segmentation
%     algorithm to reduce the island area.
%
% max_segmentation_area : integer
%     If use_iterative_segmentation is "true" then the iterative segmentation
%     algorithm will continue to reduce the size of large islands
%     until all islands have an area that is smaller than this value.
%
% show_density_map : boolean
%     Whether or not to display a figure of the density map
%
% show_mask : boolean
%     Whether to show the binary mask
%
classdef FindClustersStruct
    properties
        i3file = '';
        image_width = 256;
        image_height = 256;
        use_drift_corrected_xy = true;
        use_channels = -1;
        original_pixel_size = 160;
        analysis_pixel_size = 10;
        sum_roi_size = 5;
        sum_threshold = 5;
        factor = 5;
        localization_precision = 9.0;
        minimum_molecules_per_cluster = 5;
        ignoreNumPeakThreshold = 20e3;
        ignoreNumIslandThreshold = 15e3;
        drawUsingZRange = false;
        use_iterative_segmentation = false;
        max_segmentation_area = 20000;
        show_density_map = false;
        show_mask = false;
    end
    methods
        function write(self, filename, delimiter)
        % function write(filename, delimiter)
        %
        % Creates/Overwrites a file to write the key:value pairs to
        %
        % filename : string
        %   the name of the output file
        %
        % delimiter : string
        %   the delimiter to use to separate the key:value pair
        %   e.g., '\t' for tab delimited
            p = properties(self);
            fid = fopen(filename, 'w');
            for i = 1:length(p)
                key = p{i};
                value = self.(key);
                if isa(value, 'Insight3')
                    value = value.filename;
                elseif ~isa(value, 'char')
                    value = num2str(value);
                end
                fprintf(fid, '%s%s%s\n', key, sprintf(delimiter), value);
            end
            fclose(fid);
        end
    end    
end