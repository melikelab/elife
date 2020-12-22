function [ct,xout] = bincalc(data,numContigBins,normalize)
%%
% ORIGINAL FILENAME = bin_calc6

% the 'data' variable is an array to be binned
% the 'numContigBins' variable is the desired number of continuous zeros to
% be skipped, typically 2 or 3

%This is the final iteration of my attempt to make histograms that have a
%continuous number of bins containing data.
%What I want is for floor(sqrt(N)) number of bins to be continuous and have
%data in them, so that outliers don't make the bins un-necessarily wide.

if min(size(data)) ~= 1
    error('the ''data'' input must be a column or row vector, not an array');
end

if ~exist('normalize','var'), normalize = 0; end

if ~exist('numContigBins','var') || isempty(numContigBins)
    numContigBins = 2; %default in case it's not specified
%    disp('Defalut number of 2 bins with zeros set as target')
end
%first attempt is simply to try making a histogram with the number
%floor(sqrt(N)) of bins
desiredBins = round(sqrt(max(size(data))));
if desiredBins >= 100
    range = 2;
else
    range = 1;
end
k=1; n=1;
trialBins(k) = desiredBins;
[ct,xout] = hist(data,trialBins(k));
contiguousBins(k) = contiguous(ct,numContigBins); %This script counts contiguous bins

finished = 0;
maxloops = 8000;
loopsrun = 1;
while finished == 0
    % 120326: added the if statement for == and the elseif loopsrun > 50
    % previously, the elseif loopsrun > 50 was the if statement, s.t.:
    %if contiguousBins(k)>=desiredBins-range && contiguousBins(k)<=desiredBins+range
    %    finished = 1;
    %elseif loopsrun>=maxloops
    if loopsrun == 31, trialBins(k) = desiredBins; end
    if contiguousBins(k) == desiredBins
        finished = 1;
    elseif loopsrun > 30 && contiguousBins(k)>=desiredBins-range && contiguousBins(k)<=desiredBins+range
        finished = 1;
    elseif loopsrun>=maxloops
        finished = 2;
        disp(['DID NOT CONVERGE AFTER ' num2str(maxloops) ' FITTING LOOPS'])
        plot(trialBins, contiguousBins)
        title('Behavior of contiguous bins')
        ylabel('contiguous bins')
        xlabel('number of bins')
        break
    else
        k=k+1;
        trialBins(k) = trialBins(k-1)+1;
        [ct,xout] = hist(data,trialBins(k));
        contiguousBins(k) = contiguous(ct,numContigBins);
        loopsrun=loopsrun+1;
    end
    if k > 150
        chg(n) = contiguousBins(k) - contiguousBins(k-150); %how much is the # of contiguous bins changing
        if n>10
            if sum(chg(n-10:n)) == 0 % no change for a long while
                mostContig = max(contiguousBins);
                k_idx = find(mostContig == contiguousBins,1); %index where max occurs
                mostBins = k_idx+desiredBins-1;
                [ct,xout] = hist(data,mostBins);
%                 disp(['DESIRED NUMBER OF CONTIGUOUS BINS NOT ATTAINABLE, MAXIMUM CONTIGUOUS BINS IS ' num2str(mostContig)])
%                 disp(['TOTAL NUMBER OF BINS USED IN HISTOGRAMING = ' num2str(mostBins)])
                %                plot(desiredBins:trialBins(k), contiguousBins)
                %                title('Behavior of contiguous bins')
                %                ylabel('contiguous bins')
                %                xlabel('number of bins')
                finished = 3;
            end
        end
        n=n+1;
    end
end

if normalize == 1
    spacing = xout(2)-xout(1);
    N = max(size(data)); % number of data points = sum(ct)
    ct = ct./(spacing*N); %normalize count array such that sum(ct)*spacing=1
end
%disp(['desired number of bins = ' num2str(desiredBins)])
%if finished == 1 || finished == 2
%    disp(['obtained number of bins = ' num2str(contiguousBins(k))])
%end
%disp(['loops run = ' num2str(loopsrun)])


    function [count] = contiguous(ct,n)
        
        count = 0;
        j=1;
        ros=[]; ros(1)=0;
        sequent = 0;
        
        for i=1:length(ct)
            sep = 0;
            if ct(i)~=0     %looks to see if a count = 0 for a bin
                count = count+1;
                sequent = 0; %this resets the 'sequent' var, records sequential events of ct=0
            else
                if count == 0
                    continue;
                else
                    ros(j)=i;  %'ros' is an array of the indicies in ct where ct=0
                    if j>1
                        sep = ros(j)-ros(j-1);%this is the separation btw bins where ct=0
                        if sep==1 %sep=1 indicates that two neighboring bins both have ct=0
                            sequent=sequent+1;
                            if sequent >= n %This will stop procedure if n sequential
                                break       %zeros are encountered in the 'ct' vector
                            else
                                j=j+1;
                                continue
                            end
                        end
                    end
                    j=j+1;
                end
            end
        end
    end

end