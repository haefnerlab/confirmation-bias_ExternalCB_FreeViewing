function events = sampleNHPPEvents(rates)
%sampleNHPPEvents sample event times from a non-homogeneous Poisson Process
%(NHPP). Sampled times have units 'bins' based on the number of discrete
%bins in 'rates'.
assert(all(rates >= 0), 'NHPP requires non-negative rates');
nbins = length(rates);
meanrate = mean(rates);
h_events = sampleHPPEvents(meanrate, nbins);
inv_time_warp = cumsum(rates) / sum(rates) * nbins;

    function tt = time_warp(t)
        [~, x0] = closest(inv_time_warp, t);
        tt = fzero(@(tp) interp1(inv_time_warp, tp, 'spline', inf) - t, x0);
    end

events = arrayfun(@time_warp, h_events);
end