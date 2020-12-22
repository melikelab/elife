%function [results, i3Out] = FindClusters(findClustersStruct)
%
% Find clusters within an Insight3 molecules list. 
% Originally written by Carlo. 
%   - modified to perform an iterative segmented algorithm
% Modified by Joe to:
%   - select the Insight3 channels to use
%   - create a new Insight3 molecule list with the cluster information
%   - improve the speed of finding locatizations in each island
%   - improve the speed of the "sum of gaussian" calculation
%
% Inputs
% ------
% findClustersStruct : a custom-written struct
% (see the documentation in FindClustersStruct.m)
%
% Returns
% -------
% results : an N x 12 matrix of the cluster information, N corresponds to the number of clusters found
%     Column 1 : the X position value of each cluster, in pixels
%     Column 2 : the Y position value of each cluster, in pixels
%     Column 3 : the number of molecules in each cluster
%     Column 4 : the standard deviation of the x position value of the molecules that are part of each cluster
%     Column 5 : the standard deviation of the y position value of the molecules that are part of each cluster
%     Column 6 : the average of column 4 and 5
%     Column 7 : column 4 and 5 added in quadrature, ie. sqrt( (Column4)^2 + (Column5)^2 )
%     Column 8 : the Z values
%     Column 9 : the standard deviation for the Z values
%     Column 10 : the nearest neighbour distance, in nm, within a particular island. 
%                 A value of Inf means that only 1 cluster was in that island.
%     Column 11 : the number of clusters that were found in this island
%     Column 12 : the island index number that this cluster was found in
%
% i3Out : an Insight3 object
%    Contains the molecule list for each cluster. The centroids of each cluster are in channel 0 and the
%    localizations associated with each cluster are put in channels 1 through 9
%
function [results, i3Out] = FindClusters(findClustersStruct)                        

%%% Edits:
% 2015-August-07 Joe: implemented Carlo's iteratively-segment algorithm
%
% 2015-August-04 Joe: improve speed of peak finding and sum of gaussians
%
% 2015-July-22 Joe: added the nearest neighbour distance and information about
% the clusters in each island (columns 10, 11, 12 of results)
%
% 2015-July-08 updated the parsing of the 'i3file' input parameter so the
% program will determine if the input is a .bin file or a molecule list
% that's already been opened
% - J.Otterstrom
%
% 2015-June-28 included a second threshold value to avoid large summations
% - J.Otterstrom
%

if ~isa(findClustersStruct,'FindClustersStruct')
    error('The input variable to FindClusters is not of FindClustersStruct data type');
end

i3file = findClustersStruct.i3file;
width = findClustersStruct.image_width;
height = findClustersStruct.image_height;
useCorr = findClustersStruct.use_drift_corrected_xy;
useChannel = findClustersStruct.use_channels;
original_pixel_size = findClustersStruct.original_pixel_size;
analysis_pixel_size = findClustersStruct.analysis_pixel_size;
roi = findClustersStruct.sum_roi_size;
threshold = findClustersStruct.sum_threshold;
factor = findClustersStruct.factor;
precision = findClustersStruct.localization_precision;
minCluster = findClustersStruct.minimum_molecules_per_cluster;
ignoreNumPeakThreshold = findClustersStruct.ignoreNumPeakThreshold;
ignoreNumIslandThreshold = findClustersStruct.ignoreNumIslandThreshold;
use_iterative_segmentation = findClustersStruct.use_iterative_segmentation;
max_area = findClustersStruct.max_segmentation_area;
show_density_map = findClustersStruct.show_density_map;
show_mask = findClustersStruct.show_mask;

results = [];
i3Out = Insight3('bin');

%if exist(i3file, 'file') ~= 2
%    fprintf('File does not exist\n\t%s\n', i3file);
%    return
%end

subFolder = 'cluster_analysis'; % sub folder to save the output files to

G=@(x,y,x0,y0,rx) (1./((2*pi)*(rx*rx))).*exp(-(((x-x0)).^2)/(2*rx^2)).*...
    exp(-(((y-y0)).^2  )/(2*rx^2)); %% FUNCTION HANDLE

BW = true(width, height);

NND=[];
SIG=[];
SIGQUAD=[];
SIGX=[];
SIGY=[];
SIGZ=[];
X=[];
Y=[];
Z=[];
NN=[];
LOCS={};
NCI=[]; % number of clusters in an island
IslandIdx=[]; % the island index

updateStatus = false;

%SI1=[];
%SI2=[];
%NI1=[];
%NI2=[];

zoom_factor = original_pixel_size/analysis_pixel_size;

% check to see if a molecule list was input
if isa(i3file,'Insight3')
    % then input is actually a molecule list; change the handle reference
    % from 'i3file' to 'i3'
    i3 = i3file;
    disp('Using input Insight3 Molecule list')
elseif ischar(i3file) && (~isempty(strfind(i3file,'.bin')) || ~isempty(strfind(i3file,'.txt'))) 
    % then input is correctly formatted ML file location, read it
    if exist(i3file, 'file') ~= 2
        fprintf('File does not exist\n\t%s\n', i3file);
        return
    else
        i3 = Insight3(i3file);
    end
else
    % then the user input something for i3file with an unknown format
    error('input ''i3file'' is of unrecognized format')
end

% Joe 2016.05.02 added to make sure the molecule list is a bin format
if ~i3.isBinFormat()
    i3.convertFormat();
end

% removed 150708 ~JO % read the molecule list
% i3 = Insight3(i3file);
if isempty(i3.data)
    fprintf('Molecule list is empty\n\t%s\n', i3.getFilename());
    return
end

% get the user-specified channels from the molecule list
if useChannel > -1
    [mols, indices] = i3.getChannels(useChannel);
else
    mols = i3.getData();
    indices = transpose(1:i3.numMolecules);
end

if isempty(mols)
    fprintf('Molecule list for the specified channel is empty\n\t%s\n', i3.getFilename());
    return
end

% create the x,y locaizations in the magnified canvas
Xn = zeros(size(indices,1),3);
Xn(:,3) = indices(:); % the row indices of the molecule list, used for keeping track of which localizations are in each cluster
if useCorr
    if i3.isBinFormat()
        Xn(:,1:2) = mols(:,3:4)*zoom_factor;
    else
        Xn(:,1:2) = mols(:,4:5)*zoom_factor;
    end
else
    if i3.isBinFormat()
        Xn(:,1:2) = mols(:,1:2)*zoom_factor;
    else
        Xn(:,1:2) = mols(:,2:4)*zoom_factor;
    end
end
clear mols indices

% select only the molecules that are within the magnified image (ie. x and y > 0 )
Xn = Xn(Xn(:,1) > 0 & Xn(:,1) < width*zoom_factor, :);
Xn = Xn(Xn(:,2) > 0 & Xn(:,2) < height*zoom_factor, :);

% M is a 2D histogram of the number of molecules contained in each pixel
% example,
% (120,28) 6     -> 6 localization in pixel x=120, y=28
% (158,205) 113  -> 113 localization in pixel x=158, y=205
M=sparse(round(Xn(:,1)+1),round(Xn(:,2)+1),1);

%h=fspecial('average',roi);
%%h(find(h))=1; %comment this line in order to average the number of localizations in a 'roi' x 'roi' ROI instead of just summing the number of localizations in the ROI
%Mf=imfilter(full(M),h);

% apply the summing filter

fprintf('Applying the summing filter... ');
Mf=imfilter(full(M),ones(roi));
fprintf('DONE\n')
if show_density_map    
    fac_norm = (roi * analysis_pixel_size)^2;
    figure
    %imagesc(Mf'/fac_norm)
    [i,j] = find(M);
    xMin = min(i);
    xMax = max(i);
    yMin = min(j);
    yMax = max(j);
    imagesc(Mf(xMin:xMax,yMin:yMax)'/fac_norm)
    title(['Normalization factor = ' num2str(fac_norm) ' nm^2']);
    axis off
    axis('image');
    colorbar
    drawnow
end

%dlmwrite('C:/Users/jborbely/Desktop/ml.txt', transpose(Mf), 'delimiter', '\t');
%return
%hold on
%plot(Xn(:,2)+1,Xn(:,1)+1,'yo')
%hh(1)=gca;

fprintf('Applying the threshold... ')
Mf2=zeros(size(Mf));
%Mf2(find(Mf>minCluster))=1;  % modify if too dense
Mf2(find(Mf>threshold))=1;  % modify if too dense
fprintf('DONE\n')

if show_mask
    figure;
    imagesc(Mf2');
    %hh(2)=gca;
    %linkaxes(hh,'xy')
    axis off;
    axis('image');
    drawnow
end

if isempty(Mf2)
    fprintf('NO ISLANDS EXIST. THE THRESHOLD VALUE OF %d IS TOO BIG FOR %s\n', threshold, i3.getFilename())
    return
end

% segment
if use_iterative_segmentation
    [BWL, m] = segment_nn(Mf, Mf2, BW, threshold, max_area);
else
    fprintf('Finding islands... ')
    BWL=bwlabel(Mf2.*imresize(BW,[size(Mf2,1) size(Mf2,2)]),4);
    % figure(3)
    % imagesc(BWL)
    m=regionprops(BWL,Mf2,'PixelValues','Image','PixelList','BoundingBox');
    fprintf('DONE\n')
end


fprintf('Analysing %d islands:   0.0%% complete', length(m));
for j=1:length(m)
    fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*j/length(m))]);
    %disp('\n');
    
    BB=m(j).BoundingBox;
    %xlim([BB(1)-1 BB(1)+BB(3)+1])
    %ylim([BB(2)-1 BB(2)+BB(4)+1])

    %figure(100)
    %imagesc(m(j).Image)
    
    Mi=m(j).Image;
    xi=Xn(find(Xn(:,2)+1>=(BB(1)) & Xn(:,2)+1<=(BB(1))+BB(3) & Xn(:,1)+1>=(BB(2)) & Xn(:,1)+1<=(BB(2))+BB(4)),:);
    xi_b(:,1)=xi(:,1)-BB(2)+.5;
    xi_b(:,2)=xi(:,2)-BB(1)+.5;
    xi_b(:,3)=xi(:,3); % keep track of the row number in the molecule list
    xi_b(:,4) = round(xi_b(:,1)+1); % Joe added 2015-Aug-03
    xi_b(:,5) = round(xi_b(:,2)+1); % Joe added 2015-Aug-03
    
    %figure(100)
    %hold on
    %plot(xi_b(:,2)+1,xi_b(:,1)+1,'.k')
    
    %xi_c=NaN*ones(1,5);
    %xi_c=NaN*ones(1,3);
    [iok,jok]=find(Mi==1);

    % 
    % peak finding
    % if the length of iok gets to big then the program took too long
    if (length(iok) > ignoreNumPeakThreshold)
        fprintf(2,'\nWARNING! the peak finder will take too long. Island contains %d pixels, ignoreNumPeakThreshold=%d\n', length(iok), ignoreNumPeakThreshold);
        fprintf('Analysing %d islands:%5.1f%%%% complete', length(m), 100*j/length(m));
        clear Mi xi xi_b iok jok
        continue;
    end
    if (length(iok) > 10000)
        fprintf(2,'\nWARNING! Checking %d pixels for localizations, the localization finder may take a while...\n', length(iok));
        updateStatus = true;
        %figure(j)
        %imagesc(m(j).Image)
        %pause(0.5); % wait for the plot to show
    end
    xi_c = findLocsInIsland(xi_b, iok, jok, updateStatus);
    if updateStatus
        fprintf('\nAnalysing %d islands: %5.1f%% complete', length(m), 100*j/length(m));
        updateStatus = false;
    end
    %for i=1:length(iok)
    %    %fprintf('%d, %d\n', iok(i), jok(i));
    %    %ihigh=find( round(xi_b(:,1)+1)==iok(i) & round(xi_b(:,2)+1)==jok(i));        
    %    ihigh=find( xi_b(:,4)==iok(i) &  xi_b(:,5)==jok(i) );
    %    if ~isempty(ihigh)
    %        %fprintf('ihigh=%d\n', ihigh);
    %        xi_c=[xi_c; xi_b(ihigh,:)];
    %    end
    %    clear ihigh
    %end
    
    % ignore the NaN elements in row 1
    %xi_c=xi_c(2:end,:);
    
    % xi_c seems to be the list of molecules in "island j" but the
    % list is sorted by the order in which we search the Bounding
    % Box (BB) for molecules. Eg. the search starts in the top
    % left corner of the BB and checks each row in column 1. then
    % go to column 2 and check each row for localizations, then
    % column 3...
    
    %hold on
    %plot(xi_c(:,2)+1,xi_c(:,1)+1,'ob')
    %hold off
    %return
    
    Mi_c=zeros(ceil(size(Mi)*factor));
    % JSB 2015.08.03 comment out, calculate within the for loop below
    % If the size of the island gets too big then calculateing the Gaussian
    % over a large number of pixels is pointless since the contributions to
    % the Gaussian signal gets small after about 3 sigma from the centroid pixel location
    
    %[x,y]=meshgrid(1:size(Mi_c,2),1:size(Mi_c,1)); 
    
    % Similarly, if the length of xi_c gets to big then the program takes too long
    % 99.9% of size(xi_c,1) <= 12500
    if (size(xi_c,1) > ignoreNumIslandThreshold)
        fprintf(2,'\nWARNING! There are too many localizations on this island. Found %d localizations, ignoreNumIslandThreshold=%d\n', size(xi_c,1), ignoreNumIslandThreshold);
        fprintf('Analysing %d islands:%5.1f%%%% complete', length(m), 100*j/length(m));
        clear Mi xi xi_b iok jok xi_c
        continue;
    end 
    if (size(xi_c,1) > 10000)
        fprintf(2,'\nWARNING! There are %d localizations on this island, summing the contributions may take a while...\n', size(xi_c,1));
        updateStatus = true;
    end 
    %for i=1:size(xi_c,1)
    %    Mi_c=Mi_c+G(x, y, factor*(xi_c(i,2))+1.0+1.5, factor*(xi_c(i,1))+1.0+1.5, precision/(analysis_pixel_size/factor) );%2.5*2.5);
    %end
    
    % JSB 2015.08.03 added as a alternative way to sum gaussians
    sigma = precision/(analysis_pixel_size/factor);
    numSigma = 3;
    if updateStatus
        fprintf('\tCalculating the Sum of Gaussians for this island:   0.0%% complete');
    end        
    for i=1:size(xi_c,1)
        if updateStatus
            fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*i/size(xi_c,1))]);
        end
        x0 = factor*(xi_c(i,2))+1.0+1.5;
        y0 = factor*(xi_c(i,1))+1.0+1.5;        
        minx = max(1, floor(x0-5*sigma));
        maxx = min(size(Mi_c,2), ceil(x0+numSigma*sigma));
        miny = max(1, floor(y0-5*sigma));
        maxy = min(size(Mi_c,1), ceil(y0+numSigma*sigma));
        for xval=minx:maxx
            for yval=miny:maxy
                Mi_c(yval,xval) = Mi_c(yval,xval) + G(xval, yval, x0, y0, sigma);
            end
        end
    end
    if updateStatus
        fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100)]);
        %fprintf('\nAnalysing %d islands: %5.1f%% complete', length(m), 100*j/length(m));
        %updateStatus = false;
    end 
    
    % Mi_c is a matrix that represents the sum of 2D Gaussian Point Spread Functions
    
    clear x y
    
    %figure(4)
    %imagesc(Mi_c)
    %drawnow
    %hold off
    %return
    
    IM_BW=zeros(size(Mi_c));
    xm=[];
    ym=[];
    xm0=[];
    ym0=[];
    
    kkk=(linspace((1/10),(9/10),10));
    %kkk=(logspace(log10(1/10),log10(9/10),10));
    
    % the kkk loop finds the centroid(s) of the cluster(s) in "island j"
    if updateStatus
        fprintf('\n\tFinding the centroids for this island:   0.0%% complete');
    end        
    for k=1:length(kkk)
        if updateStatus
            fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*(k-1)/length(kkk))]);
        end
        im_bw=zeros(size(Mi_c));
        
        %im_bw(find(Mi_c>((1-k/10)*max(max(Mi_c))  )))=1;
        im_bw(find(Mi_c>((1-kkk(k))*max(max(Mi_c))  )))=1;
        
        im_bw=bwlabel(im_bw);        
        IM_BW=IM_BW+im_bw;
        
        if k==1
            s2=regionprops(im_bw,'Centroid');
            for lll=1:size(s2)
                ym=[ym round(s2(lll).Centroid(1))];
                xm=[xm round(s2(lll).Centroid(2))];
                ym0=[ym0 (s2(lll).Centroid(1))];
                xm0=[xm0 (s2(lll).Centroid(2))];
            end
        else
            %inin=[];
            inin=zeros(1,length(ym));
            for lll=1:length(ym)
                %inin=[inin im_bw(xm(lll),ym(lll))];
                inin(lll) = im_bw(xm(lll),ym(lll));
            end
            
            inin=unique(inin);
            for mm=1:length(inin)
                im_bw(find(im_bw==inin(mm)))=0;
            end
            
            im_bw(find(im_bw>0))=1;
            im_bw=bwlabel(im_bw);
            s2=regionprops(im_bw,'Centroid');
            for lll=1:size(s2)
                ym=[ym round(s2(lll).Centroid(1))];
                xm=[xm round(s2(lll).Centroid(2))];
                ym0=[ym0 (s2(lll).Centroid(1))];
                xm0=[xm0 (s2(lll).Centroid(2))];
            end
            clear s2
        end
    end
    
    if updateStatus
        fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100)]);
        fprintf('\nAnalysing %d islands: %5.1f%% complete', length(m), 100*j/length(m));
        updateStatus = false;
    end 

    
    % get the centroid of each cluster (X,Y) and the number of
    % molecules in the cluster (N)
    if ~isempty(ym0)%==0 && length(ym0)<50  % in subgraph
        %sig0=5*ones(size(xm0));
        %opt=optimset('algorithm','active-set','MaxFunEvals',1E3,'MaxIter',1E3,'Display','off');
        
        % Joe added 2015-Aug-03
        updateStatus = size(xm0,2) > 500;
        if updateStatus
            fprintf(2, '\nWARNING! There are a large number of clusters on this island. The cluster analysis for this island will converge when dlik < 1e-12\n');
        end
        
        % the function 'lik_sig' is defined below
        [xm1, N, molIDs]=lik_sig([xm0' ym0'], factor*xi_c(:,1)+1.0+1.5, factor*xi_c(:,2)+1.0+1.5, xi_c(:,3), updateStatus);
        %X=[X; xm1(:,1)/factor+BB(2)];
        %Y=[Y; xm1(:,2)/factor+BB(1)];
        
        if updateStatus
            fprintf('\nAnalysing %d islands: %5.1f%% complete', length(m), 100*j/length(m));
            updateStatus = false;
        end
        
        xyzNND = [];
        for ii=1:length(N)        
            if N(ii) >= minCluster
                mols = i3.data(molIDs{ii},:);
                if i3.isBinFormat()
                    if useCorr
                        xmols = mols(:,3);
                        ymols = mols(:,4);
                        zmols = mols(:,18);
                    else
                        xmols = mols(:,1);
                        ymols = mols(:,2);
                        zmols = mols(:,17);
                    end
                else
                    if useCorr
                        xmols = mols(:,4);
                        ymols = mols(:,5);
                        zmols = mols(:,18);
                    else
                        xmols = mols(:,2);
                        ymols = mols(:,3);
                        zmols = mols(:,17);
                    end
                end
                x_previous_way = ((xm1(ii,1) - 1.0 - 1.5)/factor - 0.5 + BB(2))*(1.0/zoom_factor);
                y_previous_way = ((xm1(ii,2) - 1.0 - 1.5)/factor - 0.5 + BB(1))*(1.0/zoom_factor);
                xmean = mean(xmols);
                ymean = mean(ymols);                
                if (abs(x_previous_way-xmean) > 0.05) || (abs(y_previous_way-ymean) > 0.05) 
                    fprintf('\n\nCentroid values disagree... contact joe.\nprevious (%f,%f), new (%f,%f)\n', x_previous_way, y_previous_way, xmean, ymean);
                    return
                end
                X = [X; xmean];
                Y = [Y; ymean];                 
                LOCS=[LOCS; mols];
                NN=[NN; N(ii)];
                xstd = std(xmols);
                ystd = std(ymols);
                SIGX=[SIGX; xstd];
                SIGY=[SIGY; ystd];
                SIG=[SIG; (xstd+ystd)/2.0];
                SIGQUAD=[SIGQUAD; sqrt((xstd*xstd + ystd*ystd)) ];
                zmean = mean(zmols);
                Z=[Z; zmean];
                SIGZ=[SIGZ; std(zmols)];
                xyzNND = [xyzNND; [xmean*original_pixel_size ymean*original_pixel_size zmean]];
                IslandIdx = [IslandIdx; j];
                
                
                % undo all the mathematical operations we did to the
                % localizations so that the cluster centroids are back in
                % the original_pixel_size canvas
                %X=[X; (xm1(ii,1)/factor - 1.0 - 1.5 + BB(2) - 0.5)*(1.0/zoom_factor)];
                %Y=[Y; (xm1(ii,2)/factor - 1.0 - 1.5 + BB(1) - 0.5)*(1.0/zoom_factor)];
                %for aaa=1:length(molIDs)
                %    LOCS=[LOCS; i3.data(molIDs{aaa},:)];
                %end

                %hold on
                %plot(xm1(:,2)/factor+BB(1),xm1(:,1)/factor+BB(2),'r*')
                %drawnow

                %hold on
                %plot(ym0/factor+BB(1),xm0/factor+BB(2),'ys')
                %drawnow

                %if size(xm1(ii),1)>1 % nearest neighbour distance
                %    nnd0=squareform(pdist(xm1(ii)))+diag(Inf*ones(size(xm1(ii),1),1));
                %    nnd=unique(min(nnd0,[],1));
                %    for iki=1:length(nnd)
                %        [i0i,j0j]=find(nnd(iki)==nnd0,1,'first');
                %        SI1=[SI1; sigm(i0i)];
                %        SI2=[SI2; sigm(j0j)];
                %        NI1=[NI1; N(i0i)];
                %        NI2=[NI2; N(j0j)];
                %        clear i0i j0j
                %    end
                %    NND=[NND; nnd'];   % value of ren_pix_size/factor, in units of nm/pix
                %else
                %    nnd=Inf;
                %    NND=[NND; nnd'];
                %end
            end
        end
        
        if ~isempty(xyzNND)
            ncisland = size(xyzNND,1);
            if ncisland > 1
                nnd0 = squareform(pdist(xyzNND)) + diag(Inf*ones(size(xyzNND,1),1));
                for jjj=1:size(nnd0,2)
                    NND = [NND; min(nnd0(:,jjj))];
                    NCI = [NCI; ncisland];                    
                end                
            else
                NND = [NND; Inf];
                NCI = [NCI; ncisland];
            end
        end        

        %clear xm1 N sigm sigQ nnd nnd0 nnd0 localizations sigx sigy
        clear N xmols ymols molIDs xstd ystd mols xmean ymean zmean xm1 nnd0
    %else
    %    X=[X; (xm0'/factor+BB(2))*(1.0/zoom_factor)]; %X=[X; xm0'/factor+BB(2)];
    %    Y=[Y; (ym0'/factor+BB(1))*(1.0/zoom_factor)]; %Y=[Y; ym0'/factor+BB(1)];
        
    %    NN=[NN; Inf*ones(size(xm0'))];
    %    SIG=[SIG; Inf*ones(size(xm0'))];
    %    SIGQUAD=[SIGQUAD; Inf*ones(size(xm0'))];
    %    SIGX=[SIGX; Inf*ones(size(xm0'))];
    %    SIGY=[SIGY; Inf*ones(size(xm0'))];
        
    %    for aaa=1:length(xm0)
    %        LOCS=[LOCS; -1.0*ones(1,18)];
    %    end
        
        %hold on
        %plot(ym0+BB(1),xm0+BB(2),'ys')
        %drawnow
    end
    
    %figure(4)
    %hold on
    %plot(ym0,xm0,'w+','MarkerSize',12)
    %hold off
    
    clear xi xi_b Mi  L iok jok ihigh Mi_c  ym0 xm0 im_bw inin nnd sigm localizations sigQ sigx sigy
    
end

fprintf('\nFound %d clusters\n', length(X));

if ~isempty(X)
    results = ones(size(NN,1), 12);
    results(:,1) = X;
    results(:,2) = Y;
    results(:,3) = NN;
    results(:,4) = SIGX;
    results(:,5) = SIGY;
    results(:,6) = SIG;
    results(:,7) = SIGQUAD;
    results(:,8) = Z;
    results(:,9) = SIGZ;
    results(:,10) = NND;
    results(:,11) = NCI;
    results(:,12) = IslandIdx;    
    [pathstr,name,~] = fileparts(i3.getFilename());
    
    if ~isequal(exist(strcat(pathstr, '/', subFolder, '/'), 'dir'),7)
        mkdir(strcat(pathstr, '/', subFolder, '/'));
    end
    
    channels = sprintf('%d,',useChannel);
    
    fprintf('Saving results to %s\\%s\\\n', pathstr,subFolder);
    if use_iterative_segmentation
        fname = sprintf('%s_roi%d_th%d_fac%g_ch%s_pix%d_prec%g_min%d_area%d.xyn', name, roi, threshold, factor, channels(1:end-1), analysis_pixel_size, precision, minCluster, max_area);
    else        
        fname = sprintf('%s_roi%d_th%d_fac%g_ch%s_pix%d_prec%g_min%d.xyn', name, roi, threshold, factor, channels(1:end-1), analysis_pixel_size, precision, minCluster);
    end
    txtFile = fullfile(pathstr, subFolder, fname);
    
    % fopen seems to only be able to open filenames with < 256 characters
    % append XXX so that it indicates that the filename is not complete
    if 256 - length(txtFile) < 0
        [p, n, e] = fileparts(txtFile);
        txtFile = [p '\' n(1:end-abs(256 - length(txtFile)-3)) 'XXX' e];
    end
    
    fprintf('Writing %s... ', fname);
    % write the parameters to the file
    findClustersStruct.write(txtFile, '\t');
    try
        % write the median N, area and density to the file        
        fle = fopen(txtFile, 'a');        
        medN = median(NN);
        fprintf(fle, 'medianNum\t%d\n', medN);
        medArea = pi*(median(SIG)*original_pixel_size)^2;
        fprintf(fle, 'medianArea[nm^2]\t%f\n', medArea);
        fprintf(fle, 'medianDensity\t%f\n', medN/medArea);
        fprintf(fle, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'X[pix]','Y[pix]','NumLocalizations','sigX[pix]','sigY[pix]','sigAve[pix]','sigQuad[pix]','Z[nm]','sigZ[nm]','NND[nm]','NumClustersInIsland', 'IslandIndex');
        fclose(fle);
        dlmwrite(txtFile, results, '-append', 'precision', '%.10g', 'delimiter', '\t')
        fprintf('DONE\n');
    catch
        fprintf('\n');
        error('ERROR! Cannot save the results. Is the following file open?\n\t%s',txtFile);
    end

    if length(fname) > 218
        fprintf(2,'WARNING! The length of the filename is > 218 characters. Excel cannot import text files with >218 charcters\n');
    end

%     % create a file to summarize the distance between cluster centroids
%     distances = repmat({''}, length(X)+2, length(X)+2);
%     for f1 = 1:length(X)
%         distances{1,f1+2} = sprintf('%d', f1);
%         distances{2,f1+2} = sprintf('(%f, %f)', X(f1), Y(f1));
%         distances{f1+2,1} = sprintf('%d', f1);
%         distances{f1+2,2} = sprintf('(%f, %f)', X(f1), Y(f1));
%         for f2 = 1:length(X)
%             distances{f1+2,f2+2} = sprintf('%f', sqrt( (X(f2)-X(f1))^2 + (Y(f2)-Y(f1))^2 ) );
%         end
%     end
%     fnameD = sprintf('%s_av%d_th%d_fac%d_ch%d_pix%d_prec%g_distances.txt', name, averagingFilterSize, threshold, factor, useChannel, render_pixel_size, localizationPrecision);
%     dlmcell([pathstr, '/', subFolder, '/', fnameD], distances)    
    
    % create an Insight3 bin file. Channel 0 will be the cluster centroids,
    % channels 1 to 9 will be used for the localizations in the cluster (colors will repeat)
    try
        %fname0 = sprintf('%s_roi%d_th%d_fac%g_ch%s_pix%d_prec%g_min%d.xyn', name, roi, threshold, factor, channels(1:end-1), analysis_pixel_size, precision, minCluster);
        fname0 = [fname(1:end-4) '.bin'];
        nMolecules=length(NN);
        for r=1:length(NN)
            nMolecules = nMolecules + size(LOCS{r},1);
        end
        binFile = fullfile(pathstr, subFolder, fname0);
        i3Out.setFilename(binFile);
        i3Out.createDefault(nMolecules);
        i3Out.forceFileOverwrite(true);
        drawUsingZRange = findClustersStruct.drawUsingZRange;
        if drawUsingZRange
            minZ = min(Z)-1; % just so that the minZ value is not at a range boundry
            maxZ = max(Z)+1; % just so that the maxZ value is not at a range boundry
            if maxZ - minZ > 2
                dz = (maxZ - minZ)/9; % create 9 Z ranges for channels 1,2,3,4,5,6,7,8,9
                zRange = minZ:dz:maxZ;
                nz = length(zRange)-1;
            else
                drawUsingZRange = false; % the are no Z ranges to use since maxZ==minZ
                channel = 1;
            end
        else
            channel = 1;
        end
        idx = 0;
        for r=1:length(NN)
            idx = idx + 1;
            i3Out.data(idx,1) = X(r);
            i3Out.data(idx,2) = Y(r);
            i3Out.data(idx,3) = X(r);
            i3Out.data(idx,4) = Y(r);
            i3Out.data(idx,12) = 0; % centroids are in channel 0
            i3Out.data(idx,17) = Z(r);
            i3Out.data(idx,18) = Z(r);
            % find which Z range this centroid is in
            if drawUsingZRange
                zval = Z(r);
                for i=1:nz
                    if (zval > zRange(i)) && (zval < zRange(i+1))
                        % localizations go in channels 1,2,3,4,5,6,7,8,9 depending in the centroid value
                        % convert this z range to the appropriate Insight3
                        % channel range so that zmin=Magenta, zmax=Red
                        if i == 1
                            channel = 5;
                        elseif i == 2
                            channel = 9;
                        elseif i == 3
                            channel = 3;
                        elseif i == 4
                            channel = 8;
                        elseif i == 5
                            channel = 6;
                        elseif i == 6
                            channel = 2;
                        elseif i == 7
                            channel = 4;
                        elseif i == 8
                            channel = 7;
                        elseif i == 9
                            channel = 1;
                        else
                            error('unknown z range');
                        end
                        break;
                    end
                end
            end            
            for v=1:size(LOCS{r},1)
                idx = idx + 1;
                i3Out.data(idx,:) = LOCS{r}(v,:);
                i3Out.data(idx,12) = channel; % localizations go in channels 1,2,3,4,5,6,7,8,9
            end
            if ~drawUsingZRange
                channel = channel + 1;
                if channel > 9
                    channel = 1;
                end
            end
        end
        i3Out.write();
    catch
        error('ERROR! Cannot save the molecule list. Is the following file open?\n\t%s',binFile);
    end
end

end



function [xx,N,molIDs]=lik_sig(xx,x,y,indices,updateStatus)

try
    
    dlik=Inf;
    %col='brgkcmybrgkcmybrgkcmybrgkcmybrgkcmybrgkcmy';

    if updateStatus
        fprintf('\tFinding clusters: dlik = Inf      ');
    end


    while dlik > 1E-12
        %clear S N SigQuad sx sy
        clear N

        for i=1:size(xx,1)
            likx(:,i)= ((x- xx(i,1)).^2);
            liky(:,i)= ((y- xx(i,2)).^2);
        end

        %close all
        %figure(1);
        %plot(xx(:,1),xx(:,2),'ob');
        %hold on
        %plot(x,y,'k');

        [l,I] = min(likx+ liky,[],2);

        lik1=sum(l);

        clear xx likx liky

        Ii=unique(I);
        molIDs = {};
        for i=1:length(Ii)
            j=find(I==Ii(i));
            molIDs{i} = indices(j);
            xx(i,1:2)=[mean(x(j));   mean(y(j))];
            %sx(i) = std(x(j));
            %sy(i) = std(y(j));
            %S(i)= (sx(i) +  sy(i))/2;
            %SigQuad(i)= sqrt(sx(i)^2 + sy(i)^2);
            N(i)=length(j);
            %plot(x(j),y(j),[col(i),'+']);
            %hold on
            clear j
        end

        clear I l Ii

        for i=1:size(xx,1)
            likx(:,i)= ((x- xx(i,1)).^2);
            liky(:,i)= ((y- xx(i,2)).^2);
        end

        [l,I] = min(likx+ liky,[],2);

        lik2=sum(l);

        clear l I likx liky

        dlik=abs(lik2-lik1);

        if updateStatus
            fprintf([repmat(sprintf('\b'), 1, 35) sprintf('\tFinding clusters: dlik = %7.3e', dlik)]);
        end

    end

catch
    clear l I likx liky Ii j S SigQuad sx sy
    fprintf(2, '\nOUT OF MEMORY. SKIPPING ISLAND.');
    xx = [];
    N = [];
    molIDs = [];
end
    
end



function xi_c = findLocsInIsland(xi_b, iok, jok, updateStatus)
    % original code
    %ihigh=find( round(xi_b(:,1)+1)==iok(i) & round(xi_b(:,2)+1)==jok(i));

    % sort by xint then yint
    [~,ind] = sortrows(xi_b(:,[4 5]));
    xi_b = xi_b(ind,:);

    % sort by iok then jok
    ijok = [iok jok];
    [~,ind] = sortrows(ijok(:,[1 2]));
    ijok = ijok(ind,:);

    nMol = size(xi_b, 1) + 1;
    nRoi = size(ijok, 1) + 1;
    locIndex = 1;
    roiIndex = 1;
    xi_c = NaN*ones(nMol,5); % allocate enough space for all localizations
    if updateStatus
        fprintf('\tFinding localizations on this island:   0.0%% complete');
    end
    cnt = 1;
    index = 0;
    while locIndex < nMol  && roiIndex < nRoi
        if updateStatus
            fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*cnt/(nMol+nRoi))]);
        end
        mx = xi_b(locIndex,4);
        my = xi_b(locIndex,5);
        rx = ijok(roiIndex,1);
        ry = ijok(roiIndex,2);
        if mx == rx && my == ry
            index = index + 1;
            xi_c(index,:) = xi_b(locIndex,:);
            locIndex = locIndex + 1;
        elseif mx > rx
            roiIndex = roiIndex + 1;
        elseif mx < rx
            locIndex = locIndex + 1;
        else
            if my < ry
                locIndex = locIndex + 1;
            else
                roiIndex = roiIndex + 1;
            end
        end
        cnt = cnt + 1;
    end
    xi_c = xi_c(1:index,:); % ignore all the NaN's from the matrix
    if updateStatus
        fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100)]);
    end
end

function [BWL, m] = segment_nn(Mf,Mf2,BW,th,max_area)

% ensure integer value
max_area = round(max_area);

%close all
%clearvars -except Mf2 BW Mf
%th=5;

BWL=bwlabel(Mf2.*imresize(BW,[size(Mf2,1) size(Mf2,2)]),4);  %%ALL
Mfn=Mf.*imresize(BW,[size(Mf2,1) size(Mf2,2)]);
Mf2n=Mf2.*imresize(BW,[size(Mf2,1) size(Mf2,2)]);

%Mf2BWL=Inf*Mf2n;

BWL2=BWL;

%m=regionprops(BWL,Mf2n,'PixelValues','Image','PixelList','BoundingBox');
%m=regionprops(BWL2,Mf2n,'Area','Image','BoundingBox');
m=regionprops(BWL2,Mf2n,'PixelValues','Image','PixelList','BoundingBox','Area');
mA=[m(:).Area];

% Joe added to provide status updates to the Command Window
currentMaxArea = max(mA);
len_str = length(num2str(currentMaxArea));

%ke=fspecial('average',5);
%ke(find(ke))=1;

% ke2=fspecial('average',2);
% ke2(find(ke2))=1;

fprintf('Finding islands using iterative segmentation (max area limit = %d), currently = %d', max_area, currentMaxArea);

%while any(mA>max_area)
while currentMaxArea > max_area
    fprintf([repmat(sprintf('\b'), 1, len_str) sprintf('%d', currentMaxArea)]);
    len_str = length(num2str(currentMaxArea));

    th=th+1;
    
    ibad=find(mA>max_area);
    for ij=1:length(ibad)
        j=ibad(ij);
        
%         figure(1)
%         imagesc(m(j).Image)
%         colormap(gray)
        %hold on
        %h(1)=gca;
        
        Mf_temp=Mfn(ceil(m(j).BoundingBox(2)):floor(m(j).BoundingBox(2))+m(j).BoundingBox(4),ceil(m(j).BoundingBox(1)):floor(m(j).BoundingBox(1))+m(j).BoundingBox(3));
        
        
        im1=zeros(size(m(j).Image));
        im1(find(Mf_temp>th))=1;
        im1=im1.*m(j).Image;
        
        
        
        BWL_temp=bwlabel(im1,4);  %%ALL
        clear im1
        BWL_temp2=zeros(size(BWL2));
        BWL_temp2(ceil(m(j).BoundingBox(2)):floor(m(j).BoundingBox(2))+m(j).BoundingBox(4),ceil(m(j).BoundingBox(1)):floor(m(j).BoundingBox(1))+m(j).BoundingBox(3))=BWL_temp;
        
        IM1=zeros(size(BWL2));
        IM1(ceil(m(j).BoundingBox(2)):floor(m(j).BoundingBox(2))+m(j).BoundingBox(4),ceil(m(j).BoundingBox(1)):floor(m(j).BoundingBox(1))+m(j).BoundingBox(3))=m(j).Image;
        
        
        %figure(1)
        %imagesc(IM1)
        %colormap(gray)
        %h(1)=gca;
        
        clear BWL_temp im1
        
        
        %         figure(3)
        %         imagesc(IM1)
        %         colormap(gray)
        %         %hold on
        %         h(2)=gca;
        %         figure(3)
        %         imagesc(BWL_temp2+IM1)
        %         colormap(lines(512*4))
        %         h(3)=gca;
        
        
        IMNones=zeros(size(BWL2));
        IMNones(find(BWL_temp2~=0))=1;
        
        while sum(sum(IM1-IMNones))>0
            %if sum(sum(IM1-IMNones))>30
            
            %return
            
            %out=nlfilter(BWL_temp2,[3 3],@mymode);
            out=mynlfilter(BWL_temp2,IM1,IMNones);
            
            
            
            IMN=out.*IM1;
            IMNones=zeros(size(IMN));
            IMNones(find(IMN~=0))=1;
            
            
            %figure(5)
            %imagesc(out.*IM1+IM1)
            %colormap(lines(512*4))
            %h(2)=gca;
            %drawnow
            
            %sum(sum(IM1-IMNones))
            BWL_temp2=IMN;
            clear IMN
            
            
            
        end
        
        
        
        
        testafter=BWL_temp2;%zeros(size(BWL_temp2));
        %testafter(find(BWL2==j))=1;
        %BWLtestafter=bwlabel(testafter,4);
        
        for ita=1:max(max(testafter))
            if ita==1
                BWL2(find(testafter==ita))=j;
            else
                BWL2(find(testafter==ita))=max(max(BWL2))+1;
            end
        end
        
        clear testafter BWLtestafter
        
        
    end
    
    %figure(5)
    %    imagesc(BWL2)
    %    colormap(lines(512*4))
    %    h(3)=gca;
    %    drawnow
    %    linkaxes(h,'xy')
        
        clear m mA
        %m=regionprops(BWL2,Mf2n,'Area','Image','BoundingBox');
        m=regionprops(BWL2,Mf2n,'PixelValues','Image','PixelList','BoundingBox','Area');
        mA=[m(:).Area];
        currentMaxArea = max(mA);
    
end
BWL=BWL2;
clear BWL2

fprintf([repmat(sprintf('\b'), 1, len_str) sprintf('%d', currentMaxArea) '\n']);

%save test_fragment_for.mat BWL th
end


function [out]=mynlfilter(BWL_temp2,IM1,IMNones)

out=BWL_temp2;
[a,b]=find(IM1-IMNones);


for i=1:length(a)
    
    x=BWL_temp2(max([a(i)-1  1]) :min([a(i)+1 size(IM1,1)]),max([b(i)-1  1]) :min([b(i)+1 size(IM1,2)]));
    
    
    x=reshape(x,[],1);
    c=find(x);
    if isempty(c)
        va=0;
    else
        va=mode(x(c));
    end
    
    out(a(i),b(i))=va;
end
end

