function [best_hprs, log_likelihoods] = xValidatePK_with_lapse_functionalPK(data, responses, oval_coordinates, time_stamps, frames, hpr_ridge, standardize, folds)
%CUSTOMREGRESSION.XVALIDATEPK Perform folded cross-validation on logistic regression with varying
%hyperparameters. Evaluates the cartesian product of hpr_ridge x hpr_ar1 x hpr_curvature. Each term
%is evaluated over 'folds' splits.

% Begin by randomly shuffling both 'data' and 'responses' so that the x-validation splits are not
% sensitive to any ordering of the data.
[trials, ~] = size(data);
shuffle = randperm(trials);
data = data(shuffle, :);
responses = 1.0 * responses(:);
responses = responses(shuffle);

% Standardize each regressor first, across the whole dataset (this is identical to the
% CustomRegression.PsychophysicalKernel function)
switch standardize
    case 0
        % do nothing
    case 1
        % assume 0 mean (nothing to subtact) and iid (std taken over all data)
        data = data / std(data(:));
    case 2
        data = zscore(data);
    otherwise
        error('Expected argument ''standardize'' to be one of [0, 1, 2]');
end

% Get start, end indices of each split
splitStart = round(linspace(1, trials+1, folds+1));
splitEnd = splitStart(2:end) - 1;

log_likelihoods = zeros(length(hpr_ridge), folds);
sz = size(log_likelihoods);

parfor ii=1:numel(log_likelihoods)
    [iRidge, iFold] = ind2sub(sz, ii);
    
    dataIdx = [1:splitStart(iFold)-1 splitEnd(iFold)+1:trials];
    [params, ~, ~, ~] = ...
    CustomRegression.PsychophysicalKernelwithlapse_functionalPK(data(dataIdx, :), responses(dataIdx),...
        oval_coordinates(dataIdx,:,:), time_stamps(dataIdx,:), frames,  ...
        hpr_ridge(iRidge), standardize);
    log_likelihoods(ii) = bernoulli_log_likelihood(data(splitStart(iFold):splitEnd(iFold), :), ...
        responses(splitStart(iFold):splitEnd(iFold)), oval_coordinates(splitStart(iFold):splitEnd(iFold),:,:), time_stamps(splitStart(iFold):splitEnd(iFold),:),frames, params);
end

avg_ll = mean(log_likelihoods, 2);
[~, imax] = max(avg_ll(:));
[iRidge] = ind2sub(sz(1), imax);
% Err on the side of less regularization by choosing smoothing that is one order of magnitude less
% than the best.
% iRidge = max(iRidge-1, 1);
best_hprs = [hpr_ridge(iRidge)];

end

function LL = bernoulli_log_likelihood(data, responses, oval_coordinates, time_stamps, num_frames, params)
eps = 0.01;
scaling = params(1)^2;
dist_param = exp(params(2));
mu = -(params(3)^2);
sigma = params(4)^2 + eps;
temporal_fixation = (params(5));% * eps;
bias = params(6);
% lapse = 1e-4+(1-1e-4)*sigmoid(params(7));%
lapse = params(7)^2;
temporal_duration = (params(8));% * eps;
trials = size(data,1);
num_im = fix(size(oval_coordinates,3)/num_frames);

for tr=1:trials
    for fx=1:num_frames
        temporal_weights(tr,((fx-1) * num_im + 1):fx * num_im) = (ones(1,num_im) * exp(temporal_fixation * fx)) * exp(temporal_duration * time_stamps(tr,fx));
    end
    elliptical_pixel_dist = squeeze(sqrt(oval_coordinates(tr,1,:).^2 + dist_param * oval_coordinates(tr,2,:).^2));
    spatial_weights(tr,:) = exp((-1*(elliptical_pixel_dist-mu).^2)./(sigma^2)) * scaling;
    logits(tr) = dot(data(tr,:),spatial_weights(tr,:) .* temporal_weights(tr,:)) + bias;
end
log_bernoulli = -1*(-1*(log(0.5*lapse+(1-lapse)*sigmoid(logits(:))).*responses)-1*(log(1-0.5*lapse-(1-lapse)*sigmoid(logits(:))).*(1-responses)));
LL = sum(log_bernoulli);
end