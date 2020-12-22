function [Mbest,MAP] = binEstimator(dataIn)
%binWidthEstimator: function to estimate optimal number of histogram bins
%for a data set
% Ref: K.H. Knuth. 2012. Optimal data-based binning for histograms and
% histogram-based probability density models, Entropy.
%   Brief:
%       In: dataIn: Vector of values to histogram
%       Out: binEdges: Vector of bin edges for histogram
%           binEdges(1) is the left edge, binEdges(end) is the last right
%           edge

sortVals = sort(dataIn);
N = length(sortVals);
MaxM = ceil(N/4);

logP = zeros(MaxM,1);
for MM = 1:MaxM
    n = hist(dataIn,MM);
    part1 = N*log(MM) + gammaln(0.5*MM) - gammaln(N+0.5*MM);
    part2 = -MM*gammaln(0.5) + sum(gammaln(n+0.5));
    logP(MM) = part1 + part2;
end

[MAP, Mbest] = max(logP);

end

