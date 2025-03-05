% Authors: Michael Dror Keselman, Alexander Friedman, Leif Gibb Garrison, Ann Graybiel. Created 2012-2014. Copyright M.I.T.

function below_threshold = compute_incompleteness(peaks, tval)
%COMPUTE_INCOMPLETENESS Computes the incompleteness grade of the cluster.
%   below_threshold = COMPUTE_INCOMPLETENESS(peaks, tvals) returns the
%   percent of the cluster theoretically below threshold (i.e., how
%   incomplete the cluster is because of the threshold).
%
%   'peaks' is a vector of the peaks for spikes in the representative wire
%   for this cluster.
%
%   'tval' is the threshold values for the representative wire in
%   microvolts.

    [n, xout] = hist(peaks, 21);
    [~, max_n] = max(n);
    delta = xout(2) - xout(1);
    theo_bottom_half = xout(max_n + 1:end) - (xout(end) - xout(max_n)) - delta;
    theo_vals = [theo_bottom_half xout(max_n) xout(max_n + 1:end)];
    half = n(max_n+1:end);
    theo_n = [half(end:-1:1) n(max_n) half];
    below_threshold = sum(theo_n(theo_vals < tval))/sum(theo_n);

end

