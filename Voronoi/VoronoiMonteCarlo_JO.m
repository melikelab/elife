% Inputs:
%   xy - array containing [x,y] positions of localizations in columns 1 & 2.
%       If more columns are present, it is assumed that column 3 is the
%       voronoi area corresponding to each localization 
%       It is assumed that the units of [x,y] are in pixels, but this can
%       be adjusted using the pixVals input
%   iter - integer value for the number of monte carlo iterations to
%       perform
%   signif - number ranging from 1-99 for calculation of the confidence
%       bounds of the Monte Carlo simulation
%   pixVals - two-element array of order: [pix2nm, finalpix] 
%       pix2nm is the conversion of [x,y] coordinates to nanometers
%       finalpix is the size of a binary Super Res mask generated for the
%           MC simulation; smaller pixel sizes better approximate the
%           actual area covered by the input localizations, but also
%           increase computational time
%       defalut values are [1 1] which leaves all calculations in units of
%       pixels for the input [x,y] coordinates
%   areaLim - optional input argument to limit the voronoi areas used for
%       calculations

function [ Histograms, intersection, xy, neighborList, mask, DT,VorDat, repidx, Varea_rnd] = ...
                            VoronoiMonteCarlo_JO(xy,iter,signif,pixVals)
%%
npts = size(xy,1);
% generate binary mask from localizations

if ~exist('pixVals','var') || isempty(pixVals)
    finalpix = [];
    pix2nm = [];
elseif ~iscolumn(pixVals) && ~isrow(pixVals)
    error('unrecognized pixel values input; should be [pix2nm, finalpix]')
else
    finalpix = pixVals(2);
    pix2nm = pixVals(1);
end
[mask.BW,mask.Area] = Locs2Mask( xy(:,1:2), finalpix, pix2nm );
% %%
% Extract voronoi areas from input points
if size(xy,2) == 2
    [xy,DT,VorDat,repidx,neighborList] = VoronoiAreas(xy,true);
else
    neighborList = [];
    repidx = [];
end
Vareas = xy(~isnan(xy(:,3)),3);
areaLim = interpercentilerange(Vareas,[ .95 .995]);
areaLim = areaLim(1);
% % set the histogram binning limits 
Vareas = Vareas(Vareas<areaLim);

% apply the uniform distrib thresh to get bin edges, then expand them to 
% cover the data range
unithresh = 2*((finalpix/pix2nm)^2)*mask.Area/size(xy,1);
[~, centers] = histcounts(Vareas(Vareas<unithresh), 'BinMethod','fd');
stp = centers(2)-centers(1);
centers = centers(1):stp:max(Vareas);
[counts, centers] = histcounts(Vareas, centers);

% set the number of random localizations to be defined
szFOV = sqrt(numel(mask.BW)); % this way it's easy to handle non-square regions
Atot = szFOV^2;
areaRatio = Atot/mask.Area; % as finalpix decreases, areaRatio increases
N_rand_pts = ceil(areaRatio*npts);
% %%
% dilate the mask to reduce edge effects in voronoi areas of randomly
% distributed localizations
se = strel('square',4);
maskDil = imdilate(mask.BW,se);
   
% find x,y pixels in the mask where localizations were found
[yD,xD] = find(maskDil); % in the dilated mask
% find dilated masked area extrema
minxD = min(xD); maxxD = max(xD); minyD = min(yD); maxyD = max(yD);

%% % if exist ('iter', 'var') && iter > 0
h = waitbar(0, 'Voronoi Monte Carlo simulation');
counts_rnd = nan(iter,length(counts));
Varea_rnd = cell(iter,1);
for j = 1:iter
    %% generate random localizations throughout the masked FOV
    
    ptsRnd = [rand(N_rand_pts, 2) * szFOV, nan(N_rand_pts,1)];
    % first, trim random locs to a box encompasing the masked area
    ptsRnd = ptsRnd( ptsRnd(:,1)>=minxD,: );
    ptsRnd = ptsRnd( ptsRnd(:,1)<=maxxD,: );
    ptsRnd = ptsRnd( ptsRnd(:,2)>=minyD,: );
    ptsRnd = ptsRnd( ptsRnd(:,2)<=maxyD,: );
    %%
    % ptsRand has x,y values in units of finalPix
    [maskSq,~,pixLocs] = Locs2Mask(ptsRnd(:,1:2),[],[],false);
    
    %% select random locs inside the dilated mask
    kp = pixLocs( maskDil );
    kp = cell2mat(kp);
    kp = sort(kp);
    
    %% keep track of original indicies
    idxOrig = ismember(1:size(ptsRnd,1),kp);
    idxOrig = transpose(idxOrig);
    
    % calculate the Voronoi areas of the 'kept', randomized localizations
    ptsSent = VoronoiAreas((finalpix/pix2nm)*ptsRnd(kp,1:2),false,false);
    % columns 1,2 in units of pixels, column 3 in units of pix^2
    % adjust indicies
    ptsRnd(idxOrig,3) = ptsSent(:,3);
    % columns 1,2 in units of finalPix-pix, column 3 in units of pix^2
    %%
    % select random locs inside the original mask
    kp = pixLocs( mask.BW );
    kp = cell2mat(kp);
    kp = sort(kp);
    
    % choose the random localizations to keep
    ptsRnd = ptsRnd(kp,:);

    %% determine histogram for random localizations
    Varea_rnd{j} = ptsRnd(~isnan(ptsRnd(:,3)),3);
    [counts_rnd(j,:), ~] = histcounts(Varea_rnd{j},centers);
    % update MonteCarlo wait bar
    waitbar(j/iter, h, 'Voronoi Monte Carlo simulation');
    %         %%
end

centers = stp*0.5+centers(1:end-1);

MeanCounts = mean(counts_rnd,1);
StdCounts = std(counts_rnd,1,1);
% z for given significance level (signif %)
z = norminv(1-(100-signif)*0.01/2);
upperConf = MeanCounts + z * StdCounts;
lowerConf = MeanCounts - z * StdCounts;

close(h);
% %% 
% find intersection of data histogram and random points histogram
[xi, yi] = polyxpoly(centers, counts, centers, MeanCounts);
intersection = [xi, yi];
% expect at most two intersections, if there are more I presume
% it's from a noisy signal and so the data must first be smoothed
% to find intersections that may or may not exist
if size(intersection,1) > 2
    counts_sm = slideFilter('avg',counts,4,2);
    MeanCounts_sm = slideFilter('avg',MeanCounts,4,2);
    [xi, yi] = polyxpoly(centers, counts_sm, centers, MeanCounts_sm);
    intersection = [xi, yi];
end

    
%     %%
Histograms = [centers', counts', MeanCounts', (lowerConf)', (upperConf)'];

end % of function
