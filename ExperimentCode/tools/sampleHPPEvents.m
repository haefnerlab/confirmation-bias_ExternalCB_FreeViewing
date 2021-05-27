function events = sampleHPPEvents(rate, maxtime)
%sampleHPPEvents sample event times from a homogeneous Poisson Process
%(HPP). Event times and ' have same units as denominator of 'rate' (e.g. if
%'rate' is in events per millisecond, then returned 'events' will be ms
%timestamps).
assert(rate >= 0, 'Poisson Process requires a non-negative rate.');
events = [];
lastevent = 0;
while true
    nextevent = lastevent + exprnd(1 / rate);
    if nextevent > maxtime
        break
    else
        events = [events nextevent];
        lastevent = nextevent;
    end
end

end