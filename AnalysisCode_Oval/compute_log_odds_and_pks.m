function [logits, spatial_pk, temporal_fixation_pk, temporal_duration_pk, temporal_pk_across_trials] = compute_log_odds_and_pks(data, oval_coordinates, time_stamps,num_frames, params, bins, trial_duration)
eps = 0.01;
scaling = prctile(params(:,1).^2,50);
dist_param = prctile(exp(params(:,2)),50);
mu = prctile(-1*(params(:,3).^2),50);
sigma = prctile(params(:,4).^2 + eps,50);
temporal_fixation = prctile(params(:,5),50);% * eps,50);
bias = prctile(params(:,6),50);
% lapse = 1e-4+(1-1e-4)*sigmoid(params(7));%
lapse = prctile(params(:,7).^2,50);
temporal_duration = prctile(params(:,8),50);% * eps,50);
disp(['Parameters [scaling dist_param mu sigma temporal-fixation bias lapse temporal-duration]: ' num2str([scaling dist_param mu sigma temporal_fixation bias lapse temporal_duration])]);

trials = size(data,1);
num_im = fix(size(oval_coordinates,3)/num_frames);
for tr=1:trials
    for fx=1:num_frames
        temporal_weights(tr,((fx-1) * num_im + 1):fx * num_im) = (ones(1,num_im) * exp(temporal_fixation * fx)) * exp(temporal_duration * time_stamps(tr,fx));
        temporal_weights_across_trials(tr,fx) = exp(temporal_fixation * fx) * exp(temporal_duration * time_stamps(tr,fx));
    end
    elliptical_pixel_dist = squeeze(sqrt(oval_coordinates(tr,1,:).^2 + dist_param * oval_coordinates(tr,2,:).^2));
    spatial_weights(tr,:) = exp((-1*(elliptical_pixel_dist-mu).^2)./(sigma^2)) * scaling;
    logits(tr) = dot(data(tr,:),spatial_weights(tr,:) .* temporal_weights(tr,:)) + bias;
end
mid_x = 1920/2; 
mid_y = 1080/2;
mid_x = mid_x/1920;
mid_y = mid_y/1080;
x_bins = linspace(0,1920,bins(1)+1);
x_bins = x_bins/1920;
y_bins = linspace(0,1080,bins(2)+1);
y_bins = y_bins/1080;
for xx=1:length(x_bins)
    for yy=1:length(y_bins)
        dist = squeeze(sqrt((x_bins(xx)-mid_x)^2 + dist_param * (y_bins(yy)-mid_y)^2));
        spatial_pk(yy,xx) = (exp((-1*(dist-mu).^2)./(sigma^2)));
    end
end
temporal_fixation_pk = exp(linspace(1,num_frames,num_frames) * temporal_fixation);% * scaling;
temporal_duration_pk = exp(linspace(0,trial_duration,10) * temporal_duration);
temporal_pk_across_trials = prctile(temporal_weights_across_trials * scaling,50);
end
