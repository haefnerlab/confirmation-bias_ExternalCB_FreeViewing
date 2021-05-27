function [saccades, fixations] = simpleSaccadeDetect(tMS, xDVA, yDVA, xDPS, yDPS, options)
%SIMPLESACCADEDETECT simple saccade detection

if ~exist('options', 'var'), options = struct; end
if ~isfield(options, 'slow_speed_threshold'), options.slow_speed_threshold = .20; end
if ~isfield(options, 'min_fixation_ms'), options.min_fixation_ms = 15; end
if ~isfield(options, 'min_ballistic_ms'), options.min_ballistic_ms = 10; end
if ~isfield(options, 'max_ballistic_ms'), options.max_ballistic_ms = 1000; end
if ~isfield(options, 'drift_max_theta_speed'), options.drift_max_theta_speed = 75; end
if ~isfield(options, 'min_time_ms'), options.min_time_ms = -inf; end
if ~isfield(options, 'max_time_ms'), options.max_time_ms = +inf; end

msPerBin = mean(diff(tMS));

% Compute eye velocity with a +/-10ms window if not already given
if ~exist('xDPS', 'var') || isempty(xDPS)
    h = ceil(10 / msPerBin);
    [xDPS, yDPS] = estimateEyeVelocity(tMS, xDVA, yDVA, h);
end

% Convert velocity to polar coordinates (r = speed, th = turn)
[th, r] = cart2pol(xDPS, yDPS);

% Smooth the speed term somewhat using an 11-ms window
r = smooth(r, ceil(11 / msPerBin / 2));
r(isnan(th)) = nan;

% Throw away velocity data outside of bounds
% TODO - fix this
r(tMS < options.min_time_ms) = nan;
r(tMS > options.max_time_ms) = nan;
th(tMS < options.min_time_ms) = nan;
th(tMS > options.max_time_ms) = nan;

% Start by detecting fixation periods
isFixating = r < options.slow_speed_threshold;

% Some pre-clean-up asserting that whenever the 'theta' component is changing rapidly, it must be
% fixational drift.
dth = abs(circ_dist(th(2:end), th(1:end-1)));
isDrifting = find(dth > deg2rad(options.drift_max_theta_speed)) + 1;
isFixating(isDrifting) = true;

% Treat missing data as a fixation so that the transition from "missing" to "fixating" is not seen
% as "gaining fixation" i.e. the end of a saccade
%isFixating(isnan(r)) = true;

% Simply set saccade onset/offset based on changes in fixation
possible_onsets  = find(diff(isFixating) == -1) + 1;
possible_offsets = find(diff(isFixating) == +1) + 1;

% At this point, there may be some off-by-one errors where there are offsets without paired onsets
% or vice versa. Loop over onsets and pair each one with the next offset.
onsets = [];
offsets = [];
j = 1;
for i=1:length(possible_onsets)
    idxPairedOffset = find(possible_offsets > possible_onsets(i), 1);
    if ~isempty(idxPairedOffset)
        duration = tMS(possible_offsets(idxPairedOffset)) - tMS(possible_onsets(i));
        if duration > 0 && duration < 500
            onsets(j) = possible_onsets(i);
            offsets(j) = possible_offsets(idxPairedOffset);
            j = j+1;
        end
    end
end

%% The above method will have a lot of false alarms. Clean them up.

% 1. Saccades must be > 5 ms in duration; < 150 ms (or user-set threshold)
durations = tMS(offsets) - tMS(onsets);
tooQuick = durations < options.min_ballistic_ms;
tooSlow  = durations > options.max_ballistic_ms;
onsets = onsets(~tooSlow & ~tooQuick);
offsets = offsets(~tooSlow & ~tooQuick);

% 2. Inter-saccadic intervals must be at least 20 ms (or user-set threshold)
intervals = tMS(onsets(2:end)) - tMS(offsets(1:end-1));
tooShortISI = find(intervals < options.min_fixation_ms) + 1;
onsets(tooShortISI) = [];
offsets(tooShortISI) = [];

% Sanity check
assert(length(onsets) == length(offsets));

%% Loop over surviving onset/offset times to create set of saccade and fixaiton informatino

saccades = struct('onsetMS', {}, 'offsetMS', {}, 'durationMS', {}, 'meanVel', {}, 'peakVel', {}, 'amplitude', {});
fixations = struct('onsetMS', {}, 'offsetMS', {}, 'durationMS', {}, 'meanX', {}, 'meanY', {});

lastOffset = 0;
for iSaccade=1:length(onsets)
    fixations(iSaccade).onsetMS = tMS(lastOffset + 1);
    fixations(iSaccade).offsetMS = tMS(onsets(iSaccade) - 1);
    fixations(iSaccade).durationMS = fixations(iSaccade).offsetMS - fixations(iSaccade).onsetMS;
    fixations(iSaccade).meanX = nanmean(xDVA(lastOffset+1:onsets(iSaccade)-1));
    fixations(iSaccade).meanY = nanmean(yDVA(lastOffset+1:onsets(iSaccade)-1));

    saccades(iSaccade).onsetMS = tMS(onsets(iSaccade));
    saccades(iSaccade).offsetMS = tMS(offsets(iSaccade));
    saccades(iSaccade).durationMS = saccades(iSaccade).offsetMS - saccades(iSaccade).onsetMS;
    saccades(iSaccade).meanVel = nanmean(r(onsets(iSaccade):offsets(iSaccade)));
    saccades(iSaccade).peakVel = max(r(onsets(iSaccade):offsets(iSaccade)));
    saccades(iSaccade).amplitude = hypot(xDVA(offsets(iSaccade)) - xDVA(onsets(iSaccade)), ...
        yDVA(offsets(iSaccade)) - yDVA(onsets(iSaccade)));

    lastOffset = offsets(iSaccade);
end

nSaccades = length(saccades);
if lastOffset < length(tMS)
    fixations(nSaccades + 1).onsetMS = tMS(lastOffset + 1);
    fixations(nSaccades + 1).offsetMS = tMS(end - 1);
    fixations(nSaccades + 1).durationMS = fixations(nSaccades + 1).offsetMS - fixations(nSaccades + 1).onsetMS;
    fixations(nSaccades + 1).meanX = nanmean(xDVA(lastOffset+1:end-1));
    fixations(nSaccades + 1).meanY = nanmean(yDVA(lastOffset+1:end-1));
end
end