% inputparam_script.m
% Call clusters struct and define all of the parameters
Zzyx = FindClustersStruct;
Zzyx.i3file = 'C:\tau\STORM_003_driftcorrected_list.bin';
Zzyx.image_width = 256;
Zzyx.image_height = 256;
Zzyx.use_drift_corrected_xy = true;
Zzyx.use_channels = -1;
Zzyx.original_pixel_size = 160;
Zzyx.analysis_pixel_size = 10;
Zzyx.sum_roi_size = 5;
Zzyx.sum_threshold = 5;
Zzyx.factor = 5;
Zzyx.localization_precision = 9.0;
Zzyx.minimum_molecules_per_cluster = 5;
Zzyx.ignoreNumPeakThreshold = 20e3;
Zzyx.ignoreNumIslandThreshold = 15e3;
Zzyx.drawUsingZRange = false;
Zzyx.use_iterative_segmentation = false;
Zzyx.max_segmentation_area = 20000;
Zzyx.show_density_map = false;
Zzyx.show_mask = false;

% call clustering algorithm
q = FindClusters(Zzyx)