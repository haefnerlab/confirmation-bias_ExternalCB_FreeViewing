function [weights, postVal, errors, map_ridge] = ...
    PsychophysicalKernelwithlapse_functionalPK(data, responses, oval_coordinates, time_stamps, num_frames,  hpr_ridge, standardize)
%PSYCHOPHYSICALKERNEL Regress PK
%
% [ weights, postVal, errors ] = PSYCHOPHYSICALKERNEL(data, responses)
% where data is [trials x regressors] and responses is [trials x 1] of
% booleans, computes regression weights 'weights', the value of the
% log posterior, and the estimated errors on each weight. Errors may be NaN
% if the Hessian is degenerate at the MAP solution.
%
% PSYCHOPHYSICALKERNEL(data, responses, hpr_ridge, hpr_ar1, hpr_curvature)
% performs grid search over all combinations of the given prior
% hyperparameters (arrays of values). Ridge controls the magnitude of the
% weights (a bias towards zero, i.e. ridgre regression). AR1 is a prior on
% the first temporal derivative of the weights, encouraging them to be
% close to constant. Curvature is a prior on the second temporal derivative
% of the weights, encouraging them to be smooth. For example,
% PSYCHOPHYSICALKERNEL(D, R, [0 0.5 1], [0], [10 100]) will fit using all
% of [0 0 10], [0 0 100], [.5 0 10], [.5 0 100], [1 0 10], [1 0 100] for
% the ridge/ar1/curvature hyperparameters respectively. Three additional
% return values contain the MAP estimate of each hyperparameter..
% [ w, p, e, map_ridge, map_ar1, map_curvature ] = PSYCHOPHYSICALKERNEL(...)
%
% ... = PSYCHOPHYSICALKERNEL(..., standardize) sets standardization rule.
% When standardize == 0, no preprocessing is done. When standardize == 1,
% the data are standardized assuming iid and zero-mean. When standardize ==
% 2, the data are standardized using the built-in zscore() function.
% Defaults to 0.


% Standardize each regressor.
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

% Add a column of ones to data for a bias term.
% convert boolean to float type
responses = 1.0 * responses(:);

% Grid search will be done over all combinations (i.e. the cartesian
% product) of given hyperparameters.
grid_size = [length(hpr_ridge)];
n_gridpts = prod(grid_size);

% Each entry in 'results' will itself be a cell array containing the return
% values (i.e. {weights, postval, errors, ...})
results = cell(n_gridpts, 1);

for i=1:n_gridpts
    % Determine which hyperparameters to use this iteration by treating i as a 1d index into the 3d grid of ridge/ar1/curvature values.
    [idx_ridge] = ind2sub(grid_size, i);
    
    [weights, postVal, errors] = do_fit(data, responses, oval_coordinates, time_stamps, num_frames, hpr_ridge(idx_ridge));
    
    % Record all results for this set of hyperparameters.
    results{i} = {weights, postVal, errors, hpr_ridge(idx_ridge)};
end

% Find and return the MAP result.
postVals = cellfun(@(result) result{2}, results);
[~, idx_map] = max(postVals);
[weights, postVal, errors, map_ridge] = results{idx_map}{:};

end

function [optim_weights, postVal, errors] = do_fit(data, responses, oval_coordinates, time_stamps, num_frames, hpr_ridge)
% Perform fitting at a single setting of hyperparameters / at a single
% point in the grid search.

p = 8;
prior_matrix = hpr_ridge * eye(p);

    function [nlp] = neg_log_posterior(params)
        % COPIED FROM CustomRegression.PsychophysicalKernel
        neg_log_prior = 0.5 * params' * prior_matrix * params;
        nlp = neg_log_prior + neg_bernoulli_log_likelihood(data, responses, oval_coordinates, time_stamps, num_frames,params);
    end

options = optimoptions(@fminunc, 'Display', 'off', 'MaxFunctionEvaluations', 1e5);
try
    [optim_weights, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_log_posterior, rand(8,1), options);
catch
    [optim_weights, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_log_posterior, rand(8,1), options);
end
postVal = -negPostVal;
% attempt to invert the hessian for standard error estimate - this sometimes fails silently returning NaN.
errors = sqrt(diag(inv(hessian)));
end

function LL = neg_bernoulli_log_likelihood(data, responses, oval_coordinates, time_stamps, num_frames, params)
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
log_bernoulli = -1*(log(0.5*lapse+(1-lapse)*sigmoid(logits(:))).*responses)-1*(log(1-0.5*lapse-(1-lapse)*sigmoid(logits(:))).*(1-responses));
LL = sum(log_bernoulli);
end

