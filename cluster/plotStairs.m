

function [bins,counts] = plotStairs(bindata,numBins,spacing)

% INPUTS
% bindata - vector of data to be binned for plotting using the 'stairs'
%       function
% numBins - (optional) number of bins to include in binning.
%       Input a value of 'nan' to set number of bins according to the size
%       of bindata: numBins = sqrt(length(numBins))
% spacing - (optional) the type of bin spacing being either logarithmic
%       (default) or linear

% OUTPUTS
% bins - 'x' values for the stairs function
% counts - normalized 'y' values for the stairs function such that it
%       represents the probability corresponding to each bin, not the
%       number of counts therein

% check if user input a different number of bins
if ~exist('numBins','var') || isempty(numBins)
    numBins = 25;
elseif isnan(numBins)
    numBins = round(sqrt(length(bindata)));
end

% check if the user input a type of bin spacing and if it is correctly
% formatted
if ~exist('spacing','var') || isempty(spacing)
    spacing = 'log';
elseif sum(strcmpi(spacing,{'lin','linear','log','logarithmic'}))==0
    warning('Improper "spacing" input, set to default: logarithmic')
    spacing = 'log';
    % accept writing out the full word rather than inputting the abbrev.
elseif strcmpi(spacing,'linear'), spacing = 'lin';
elseif strcmpi(spacing,'logarithmic'), spacing = 'log';
end



binstart = floor(min(bindata));
binend = ceil(max(bindata));
switch spacing
    case 'log'
        if binstart == 0
            pwr = 0;
            while binstart==0
                pwr = pwr+1;
                multfct = 1*10^pwr;
                binstart = floor(multfct*min(bindata(bindata~=0)))/multfct;
            end
        end
        binstart = log10(binstart);
        binend = log10(binend);
end
good = 0;
while good == 0
    switch spacing
        case 'log'
            bins = logspace( ...
                binstart, ...
                binend, ...
                numBins)';
            separation = diff( log10(bins) );
            separation = separation(end);
        case 'lin'
            bins = linspace( ...
                binstart, ...
                binend, ...
                numBins)';
            separation = diff( bins );
            separation = separation(end);
    end
    
    counts = histc(bindata,bins);
    if counts(end) == 0 && counts(1) == 0
        good = 1;
    elseif counts(end) ~= 0
        binend = binend + separation;
    elseif counts(1) ~= 0
        binstart = binstart - separation;
    end
end

% define the output bins & counts
switch spacing
    case 'log'
        space = diff(bins);
        idx = find(counts~=0);
        counts = counts(idx)./(space(idx).*sum(counts)); % normalize
        counts = [0; counts];
        bins = bins(idx);
        bins = [bins(1)*0.999; bins];
        
    case 'lin'
        counts = counts./(separation.*sum(counts)); % normalize
end

end