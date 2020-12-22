% The algorithm is as follows:
% Load the Insight3 molecule list and the SR-Tessler objects file.
% Sort the Insight3 localizations by the 'xc' column (required for the Binary Search algorithm).
% For each cluster
%   For each SR-Tesseler (x,y) point in this cluster
%     - perform the Binary Search algorithm to quickly find the index in 
%       the sorted Insight3 list where the x SR-Tesseler value is 
%       slightly less than the Insight3 xc value.
%     - use the index found from the Binary Search algorithm as the
%       starting index for performing a brute-force search for finding the 
%       (x,y) SR-Tesseler point in the sorted Insight3 molecule list.
%
% Programmed using MATLAB 2016a
%   Joe Borbely & Jason Otterstrom
%
% INPUTS
%   i3In: a previously loaded Insight3 object
%   srData: array of data from the exported .xls file saved using the  
%       SR-Tesseler program. Columns are: 
%       [Localization Index, Object/cluster Index, X, Y]
%   eps (optional): epsilon value to use for considering rounding errors 
%       when determining if the x values and the y values in the Insight3 
%       list and the SR-Tesseler list are equal 
%
% OUTPUT
%   SRTi3_trf: array having a length equivalent to srData and four columns
%   that are:
%   [srData_index, Object/cluster index, i3_index, separation_abs(i3-srData)]
%   Data is sorted such that srData_index is in ascending order

function SRT2i3_trf = pairSRT_I3(i3In, srData, eps)

% Check inputs
if nargin < 3 || isempty(eps)
eps = 0.001;
end

if ~isobject(i3In)
    error('The first input must be an Insight3 object')
end

if size(srData,2) ~= 4
    error('The second input is expected to have four columns')
end

%% Testing
% % load the Insight3 molecule list and the SR-Tesseler object list
% % i3In = Insight3('D:\SR-Tesseler\1mM 20min_007_list_AutoDC_0916-2490_allChs.bin');
% % srData = dlmread('D:\SR-Tesseler\ID_localizations_objects2.xls','\t',1,0);
% i3In = Insight3('F:\Testing Software\SR-Tesseler\1mM 20min_007_list_AutoDC_0916-2490_allChs.bin');
% srData = dlmread('F:\Testing Software\SR-Tesseler\ID_localizations_objects2.xls','\t',1,0);
%% Algorithm

% get the column indices that are needed from the Insight3 molecule list
% channelIndex = i3In.getColumnIndex('channel');
% xIndex = i3In.getColumnIndex('x');
% yIndex = i3In.getColumnIndex('y');
xcIndex = i3In.getColumnIndex('xc');
ycIndex = i3In.getColumnIndex('yc');

% sort the Insight3 molecule list by the xc values
[~,indices] = sortrows(i3In.data(:,xcIndex));
i3InDat_SortX = i3In.data(indices,:);

% % initialize the output Insight3 molecule list
% [p, n, e] = fileparts(i3In.filename);
% i3Out = Insight3();
% i3Out.setFilename([p '\' n '_SR-Tesseler_JOE' e]);
% i3Out.forceFileOverwrite(true);
% i3OutIndex = 1;

% for each cluster
maxCluster = max(srData(:,2));
% i3Out.createDefault(size(srData,1) + maxCluster + 1);
% channel = 1;

SRT2i3_trf = nan(size(srData,1),4);
n = 0;
% eps = 0.001;
tic
for ic = 0:maxCluster
    
    % get the localizations in this SR-Tesseler cluster
    clustIdx = find(srData(:,2) == ic);
    locsSR = srData( clustIdx, :);
    nSR = size(locsSR,1);

    fprintf('Cluster %03d of %03d (%06d localizations):   0.0%% complete', ic, maxCluster, nSR);

%     % insert the centroid to the output Insight3 molecule list
%     xmean = mean(locsSR(:,3));
%     ymean = mean(locsSR(:,4));
%     i3Out.data(i3OutIndex,xIndex) = xmean;
%     i3Out.data(i3OutIndex,yIndex) = ymean;
%     i3Out.data(i3OutIndex,xcIndex) = xmean;
%     i3Out.data(i3OutIndex,ycIndex) = ymean;
%     i3Out.data(i3OutIndex,channelIndex) = 0;
%     i3OutIndex = i3OutIndex + 1;
    
    %% for each localization in this SR-Tesseler cluster find the pairing
    % localization in the sorted Insight3 list
    for srIndex = 1:nSR
        xSR = locsSR(srIndex,3);
        ySR = locsSR(srIndex,4);

        % start with a Binary Search algorithm to determine which i3In index
        % to start at for a brute-force search
        % see https://en.wikipedia.org/wiki/Binary_search_algorithm
        L = 1;
        R = i3In.numMolecules;
        
        while true
            if L > R
                error('L > R in cluster %d', ic)
            end
            m = floor((L+R)/2);
            xI3 = i3InDat_SortX(m,xcIndex);
            if abs(xI3-xSR) < eps
                minit = m;
                % found an index where the x localization values are approximately equal.
                % step the m index backwards until the x values differ by 0.1
                % this m index will then be the value to start a brute-force search
                % need to reduce the m index since the y localization values are not sorted
                while true
                    m = m - 1;
                    if m < 1
                        m = 1;
                        break; % breaks the inner 'while true' loop
                    end
                    xI3 = i3InDat_SortX(m,xcIndex);
                    if abs(xI3-xSR) > 0.1
                        break; % breaks the inner 'while true' loop
                    end
                end
                break; % breaks the outer 'while true' loop, 
                        % m is the index to start a brute-force search
            elseif xI3 < xSR
                L = m + 1;
            elseif xI3 > xSR
                R = m - 1;
            end
        end
        mfinal = m;
%         %% I attempted to improve the algorithm, but only made it worse
%         % ~JO 160905
%         % Binary search gives a useful range for the m index to use to 
%         % perform a brute-force search for the closest corresponding
%         % localization
%         % note: a true brute-force method would always range from 1:i3In.numMolecules
%         if minit <= mfinal
%             idxRange = minit:mfinal;
%         elseif minit > mfinal
%             idxRange = mfinal:minit;
%         end
%         if idxRange(end) > i3In.numMolecules
%             fprintf('\n');
%             error('Cannot find (%f,%f) in cluster %d', xSR, ySR, ic)
%         end
%         xI3 = i3InDat_SortX(idxRange,xcIndex);
%         yI3 = i3InDat_SortX(idxRange,ycIndex);
%         [minVal,minIdx] = min(sqrt( (xI3-xSR).^2 + (yI3-ySR).^2 ));
%         m = idxRange(minIdx);
%         %%
        % got a good starting for the m index to use to start a brute-force search
        % note: a true brute-force method would always start with m=1
%         disp('Joe Method')
%         tic
        while true
            if m > i3In.numMolecules
                fprintf('\n');
                error('Cannot find (%f,%f) in cluster %d', xSR, ySR, ic)
            end
            xI3 = i3InDat_SortX(m,xcIndex);
            yI3 = i3InDat_SortX(m,ycIndex);
            if abs(xI3-xSR) < eps && abs(yI3-ySR) < eps      
                minVal = sqrt( (xI3-xSR)^2+(yI3-ySR)^2 );
                break; % search is done, m is the index where the localizations are equal (within eps)
            else
                m = m + 1;
            end
        end                
%         toc
%         m
%         %%
        % here save the transform array in format:
        % [SRT_index, SRT_cluster_ID, i3_index, SRT_i3_separation]
        n = n+1;
        SRT2i3_trf(n,:) = [clustIdx(srIndex), ic, indices(m), minVal];
%         [clustIdx(srIndex), ic, m, minVal]
        %%
        % m is the corresponding index in the sorted Insight3 molecule list
        % where this (x,y) SR-Tesseler point equals a (x,y) Insight3 point
        % insert the localization into i3Out
%         i3Out.data(i3OutIndex,:) = i3InDat_SortX(m,:);
%         i3Out.data(i3OutIndex,channelIndex) = channel;        
%         i3OutIndex = i3OutIndex + 1;
        
        fprintf([repmat(sprintf('\b'), 1, 15) sprintf('%5.1f%%%% complete',100*srIndex/nSR)]);

    end
    
    % increment the channel number for the next cluster
%     channel = channel + 1;    
%     if channel == 10
%         channel = 1;
%     end
    
    fprintf('\n');
    
end

%% Now double check to see if any localizations are doubly assigned
mult = 1;
while ~isempty(mult)
    %%
    unqidx = unique(SRT2i3_trf(:,3));
    % need two counting programs depending on the version of MATLAB
    if exist('histcounts','file')
        [N,edges]=histcounts(SRT2i3_trf(:,3),[unqidx;unqidx(end)+1]); %have to include extra bin b/c histcounts is stupid and groups things at the end
    else
        N=histc(SRT2i3_trf(:,3),unqidx);
        edges = unqidx;
    end            
    mult = edges(N>1);
    disp(['Found ' num2str(size(mult,1)) ' repeats'])
    for i = 1:size(mult,1)
        
        % get indicies where multiple instances occur
        tstidx = find( SRT2i3_trf(:,3) == mult(i) );
        % find which of the repeats is closest to the I3 localization
        [~,keepidx] = min(SRT2i3_trf( tstidx ,4));
        keepidx = tstidx(keepidx);
        % make test i3 index list that excludes the localization of interest
        testi3Idx = 1:i3In.numMolecules;
        testi3Idx = testi3Idx(~ismember(testi3Idx,mult(i)));
        % identify the indicies that must be re-assigned
        remidx = tstidx(~ismember(tstidx,keepidx));
        %%
        for j = 1:size(remidx,1)
            % extract x,y position in srData to search for other nearest Loc
            srDatIdx = SRT2i3_trf(remidx(j),1);
            xSR = srData( srDatIdx, 3 );
            ySR = srData( srDatIdx, 4 );
            
            xI3 = i3In.data(testi3Idx,xcIndex);
            yI3 = i3In.data(testi3Idx,ycIndex);
            
            [minVal,newIdx] = min( sqrt( (xI3-xSR).^2+(yI3-ySR).^2 ) );
            %         [srDatIdx, SRT2i3_trf(remidx(j),2), testi3Idx(newIdx), minVal]
            SRT2i3_trf(remidx(j),:) = [srDatIdx, SRT2i3_trf(remidx(j),2), testi3Idx(newIdx), minVal];
            
        end
    end
end
% save the Insight3 molecule list
% i3Out.write();

% sort the output transform matrix to have the SR-Tesseler localization in
% ascending order in the first column
SRT2i3_trf = sortrows(SRT2i3_trf);

toc
end