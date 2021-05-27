function GaborData = newGaborData(varargin)

if ~isempty(varargin) && isstruct(varargin{1})
    % Get first input as a 'template' struct.
    template = varargin{1};
    varargin = varargin(2:end);
else
    template = struct();
end

    function value = get_arg(name, default)
        % Helper function to get named arguments with a default
        idx = strcmpi(name, varargin);
        if any(idx)
            val_idx = find(idx)+1;
            value = varargin{val_idx};
            varargin(find(idx):val_idx) = [];
        elseif isfield(template, name)
            if any(strcmpi(name, {'stair_bounds', 'step_size', 'min_step_size'})) && ~isequal(GaborData.stair_fn, template.stair_fn)
                warning('Not copying field %s from template since template''s stair_fn is %s', name, func2str(template.stair_fn));
                value = default;
                return;
            end
            if isnumeric(template.(name))
                expected_size = length(default);
                value = template.(name)(1:expected_size);
            else
                value = template.(name);
            end
        else
            value = default;
        end
    end

%% User-settable params
GaborData.trials_per_block = get_arg('trials_per_block', 100);
GaborData.blocks = get_arg('blocks', 12); % can manually set this to 9 when running ExperimentGabor(newGaborData, "blocks", 9)
% GaborData.stair_fn = get_arg('stair_fn', @Staircase.ratio); % comment this out for contrast condition
GaborData.stair_fn = get_arg('stair_fn', @Staircase.ratio); % comment this out for ratio condition
GaborData.reversals_per_epoch = get_arg('reversals_per_epoch', 6);

total_trials = GaborData.trials_per_block * GaborData.blocks;

% Initial values of staircase-able parameters
GaborData.number_of_images = get_arg('number_of_images', 15);
%GaborData.contrast = zeros(total_trials,GaborData.number_of_images);
GaborData.contrast = zeros(1, total_trials);
GaborData.contrast_per_frame = zeros(total_trials,GaborData.number_of_images);
GaborData.contrast(1) = get_arg('contrast', 0.9);
GaborData.ratio = zeros(1, total_trials);
GaborData.ratio(1) = get_arg('ratio', 0.8);
GaborData.noise = zeros(1, total_trials);
GaborData.noise(1) = get_arg('noise', 5); % std of variance of noise added pixelwise
GaborData.step_size = zeros(1, total_trials);

% Staircase bounds and step size, with defaults set depending on stair_fn
GaborData.model_observer = get_arg('model_observer', '');
if isequal(GaborData.stair_fn, @Staircase.contrast)
    GaborData.stair_bounds = get_arg('stair_bounds', [0.5 0.8]);
    GaborData.step_size(1) = get_arg('step_size', 1.1); % multiplicative (in the "easier" direction)
    GaborData.min_step_size = get_arg('min_step_size', 1.02); % Default to two 'halvings' of the step size
    % GaborData.min_step_size = get_arg('min_step_size', 1+(GaborData.step_size(1) - 1)/4); % Default to two 'halvings' of the step size
    GaborData.test_threshold = get_arg('test_threshold', 0.55);
    GaborData.test_ratio = get_arg('test_ratio', 1.0);
elseif isequal(GaborData.stair_fn, @Staircase.ratio)
    GaborData.stair_bounds = get_arg('stair_bounds', [0.5 0.8]);
    GaborData.step_size(1) = get_arg('step_size', .05); % additive (in the "easier" direction)
    GaborData.min_step_size = get_arg('min_step_size', GaborData.step_size(1)/4); % Default to two 'halvings' of the step size
elseif isequal(GaborData.stair_fn, @Staircase.noise)
    % Note: noise is treated as a special case, where we use a discrete set
    % of values for kappa, and the staircase simply increments/decrements
    % the index.
    GaborData.noise_set = get_arg('noise_set', linspace(1,6, 21)); % Higher indices correspond to harder values of kappa.
    GaborData.stair_bounds = get_arg('stair_bounds', [1 length(GaborData.noise_set)]); % Not actually used; bounds implied by length of array. See Staircase.noise
    GaborData.step_size(1) = get_arg('step_size', 4); % additive (in the "easier" direction)
    GaborData.min_step_size = get_arg('min_step_size', 1); % Cannot step fewer than 1 indices in an array.
    GaborData.test_threshold = get_arg('test_threshold', 6);
    GaborData.test_ratio = get_arg('test_ratio', 0.9);
end

% Other misc. user-definable parameters relating to stimulus/rig.
%GaborData.flag_use_old_stimulus_code = false;  % Henceforth all stimuli are generated using 'correct' code.
% GaborData.number_of_images = get_arg('number_of_images', 48);
%GaborData.stimulus_fps = get_arg('stimulus_fps', 8);  % frame rate of stimuli
GaborData.blank_frames = get_arg('blank_frames', 5);  % number of blank screen frames per stimulus frame
%GaborData.cue_duration = get_arg('cue_duration', 0.2);  % Fixed duration, seconds to display cue after getting fixation.
GaborData.fixation_duration = get_arg('fixation_duration', 0.2); %previously 0.75 
GaborData.prefixation_duration = get_arg('prefixation_duraction',0.5);
GaborData.stim_duration = get_arg('stim_duration', 1.5);
% GaborData.annulus_radius = get_arg('annulus_radius', 500); % Size, in pixels, of hole in center of stimulus
GaborData.aspect_ratio = get_arg('aspect_ratio', 0.55);
GaborData.horizontal = get_arg('horizontal', -1);
GaborData.vertical = get_arg('vertical', 1);
GaborData.go_cue_time = get_arg('go_cue_time', 0.75);  % Time between final stimulus/mask frame and the targets appearing.
% BPG Stimulus parameters
GaborData.stim_size = get_arg('stim_size', 361);  % Width of the stimulus in pixels.
%GaborData.stim_sp_freq_cpp = get_arg('stim_sp_freq_cpp', 50);  % Mean spatial frequency of images in cycles per pixel.
%GaborData.stim_std_sp_freq_cpp = get_arg('stim_std_sp_freq_cpp', .1);  % Std deviation of spatial frequency in cycles per pixel.
GaborData.oval_scale = get_arg('oval_scale', 50); % multiplicative factor to axes ratio (in trialStimuliOval)

%hexagonal grids 
GaborData.grid_length = get_arg('grid_length', 270); % Shortest length between oval locations 
GaborData.x_grid_wdth = GaborData.grid_length*3; GaborData.y_grid_wdth = GaborData.grid_length* sqrt(3); % Length of a base of the smallest equilateral triangle on the grid. % Length of a height of the smallest equi. tri. 
GaborData.point_X = 0:GaborData.x_grid_wdth:1920; GaborData.point_Y = 0:GaborData.y_grid_wdth:1080; % Chooses points on the x axis and y axis in a sequence based on grid_wdth. 
GaborData.first_X = (1920 - (length(GaborData.point_X)-1)*GaborData.x_grid_wdth) / 2; % Calculates the first point on the x and y axes so that the grids are symmetrical.
GaborData.first_Y = (1080 - (length(GaborData.point_Y)-1)*GaborData.y_grid_wdth) / 2; 
GaborData.last_X = 1920 - GaborData.first_X; % Calculates the last point on the x and y axes so that the grids are symmetrical 
GaborData.last_Y = 1080 - GaborData.first_Y; 
GaborData.point_X = GaborData.first_X:GaborData.x_grid_wdth:GaborData.last_X; 
GaborData.point_Y = GaborData.first_Y:GaborData.y_grid_wdth:GaborData.last_Y;
GaborData.number_of_images = 2*length(GaborData.point_X)*length(GaborData.point_Y); % multiplicative factor 2 is for locating disks on both mid_in and mid_out (look at trialStimuliOval). 

% Preallocate fields that will be populated with data by running the
% experiment.
GaborData.iid = true(1, total_trials);
GaborData.streak = zeros(1, total_trials);
GaborData.reversal_counter = zeros(1, total_trials);
GaborData.ideal_answer = zeros(1, total_trials);
GaborData.reaction_time = zeros(1, total_trials);
GaborData.choice = zeros(1, total_trials);
GaborData.accuracy = zeros(1, total_trials);
GaborData.noisy_ratio = zeros(total_trials, GaborData.number_of_images);
GaborData.sig_ratio_noise = get_arg('sig_ratio_noise', 0.1);
GaborData.ideal_frame_signals = zeros(total_trials, GaborData.number_of_images);


% Note that 'seed' and 'correct_answer' must be preset due to esoteric
% properties of random number generators. GaborStimulus will read out these
% preset values.g
GaborData.seed = randi(100000000, 1, total_trials);
% GaborData.checksum = zeros(1, total_trials);  % For sanity-checks on seeds
GaborData.correct_answer = 1 * rand(1, total_trials) < .5;
GaborData.phase = zeros(GaborData.number_of_images+1, total_trials);
GaborData.phase(1, :) = 2*pi * rand(1, total_trials);
for i = 2:GaborData.number_of_images+1
    GaborData.phase(i, :) = GaborData.phase(i-1, :) + normrnd(0, 0.8, 1, total_trials);
end

GaborData.current_trial = 0;

GaborData.eye_tracker_points = {};

if isequal(GaborData.stair_fn, @Staircase.ratio)
    if GaborData.step_size(1) < 0
        warning('Changing sign of ratio step_size from %d to %d', GaborData.step_size(1), -GaborData.step_size(1));
        GaborData.step_size = -GaborData.step_size;
    end
end

if isequal(GaborData.stair_fn, @Staircase.noise)
    if GaborData.step_size(1) < 0
        warning('Changing sign of noise step_size from %d to %d', GaborData.step_size(1), -GaborData.step_size(1));
        GaborData.step_size = -GaborData.step_size;
    end
    if ~(iseffectiveinteger(GaborData.step_size(1)) && iseffectiveinteger(GaborData.min_step_size))
        error('In Noise condition, step_size and min_step_size act on indices and must be integers');
    end
end

if isequal(GaborData.stair_fn, @Staircase.contrast)
    if GaborData.step_size(1) < 0
        error('Contrast staircase is multiplicative; step size of %f doesn''t make sense', GaborData.step_size(1));
        %     elseif GaborData.step_size(1) < 1
        %         warning('Chaning contrast step_size < 1 to 1/step_size');
        %         GaborData.step_size(1) = 1 / GaborData.step_size(1);
    end
end

if GaborData.ratio(1) > 1 || GaborData.ratio(1) < 0
    error('Ratio should be between 0.5 and 1');
end

if ~isempty(GaborData.model_observer) && ~any(strcmpi(GaborData.model_observer, {'ideal', 'oracle', 'bernoulli'}))
    warning('%s is not a known model observer', GaborData.model_observer);
end

if strcmp(GaborData.model_observer, 'bernoulli')
    GaborData.sigmoid_slope = get_arg('sigmoid_slope', 20);
    default_pk = ones(1, GaborData.number_of_images) / GaborData.number_of_images;
    GaborData.model_pk = get_arg('model_pk', default_pk);
    if sum(GaborData.model_pk) ~= 1
        warning('Recommended that GaborData.model_pk sum to 1');
    end
end

if ~isempty(varargin)
    warning('Unkown arguments given to newGaborParams');
end

disp(GaborData);


    function TF = iseffectiveinteger(v)
        TF = (v == floor(v));
    end

end