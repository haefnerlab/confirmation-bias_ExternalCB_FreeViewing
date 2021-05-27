function [frame_categories, xaxis, yaxis] = GaborStimulus(GaborData)
%GABORSTIMULUS(GaborData, trial) create (or recreate) stimulus frames based
%on parameters in GaborData and the seed, contrast, ratio, and noise on the
%given trial. If 'GaborData.iid(trial)' is true, each frame's category is
%drawn iid based on the 'ratio' parameter. Otherwise, exactly
%round(ratio*num_images) frames will match the 'true' category.
%
%This function makes no modifications to GaborData.

% Set RNG state to recreate stimulus for this trail.
rng(GaborData.seed(GaborData.current_trial), 'twister');
if ~isfield(GaborData, 'iid') || GaborData.iid(GaborData.current_trial)
    % Randomly set each frame to match (or mismatch) the correct choice
    % for this trail, using the current 'ratio' to decide.
    match_frames = rand(1, GaborData.number_of_images) <= GaborData.ratio(GaborData.current_trial);
else
    % Randomly permute whether each frame matches the true category, with
    % 'ratio' percent of them matching.
    n_match = round(GaborData.ratio(GaborData.current_trial) * GaborData.number_of_images+1);
    match_frames = [true(1, n_match) false(1, (GaborData.number_of_images) - n_match)];
    match_frames = Shuffle(match_frames);
end

frame_categories = zeros(size(match_frames));

% Choose frames based on whether correct answer this trial is Left or Right
if GaborData.correct_answer(GaborData.current_trial) == 1
    frame_categories(match_frames) = GaborData.vertical;
    frame_categories(~match_frames) = GaborData.horizontal;
else
    frame_categories(~match_frames) = GaborData.vertical;
    frame_categories(match_frames) = GaborData.horizontal;
end

% Set random seed again to keep match_frames independent of pixel noise.
rng(GaborData.seed(GaborData.current_trial), 'twister');

for i=1: GaborData.number_of_images
    if frame_categories(i)== -1
        xaxis(i) = GaborData.oval_scale * GaborData.aspect_ratio;
        yaxis(i) = GaborData.oval_scale * (1 - GaborData.aspect_ratio);
    else
        xaxis(i) = GaborData.oval_scale * (1 - GaborData.aspect_ratio);
        yaxis(i) = GaborData.oval_scale * GaborData.aspect_ratio;
    end
end
end