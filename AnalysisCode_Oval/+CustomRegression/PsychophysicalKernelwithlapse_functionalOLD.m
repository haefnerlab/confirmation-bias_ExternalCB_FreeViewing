function [optim_params] = PsychophysicalKernelwithlapse_functionalOLD(data, responses, oval_coordinates, num_frames)
responses = 1.0 * responses(:);
% num_im = fix(size(oval_coordinates,3)/num_frames);
% for i=1:num_frames
%     data_temp(:,i) = data(:,(i-1)*num_im + 1);
% end
% glm_weights = glmfit(data_temp, responses, 'binomial');
% glm_weights = glm_weights(2:end)/max(glm_weights(2:end));
% ab = CustomRegression.expFit(glm_weights(2:end));
% init_vals = [rand(4,1); ab(1); glm_weights(1); rand];% ab(2)];

init_vals = [randi([2,10],1); rand(5,1); 0.001];
options = optimoptions(@fminunc, 'Display', 'off', 'MaxFunctionEvaluations', 1e5);

% % init_vals = rand(8,1);
% % options = optimoptions(@fminunc,'Display', 'off', ...
% %     'Algorithm', 'trust-region', ...
% %     'SpecifyObjectiveGradient', true, ...
% %     'HessianFcn', 'objective','MaxFunctionEvaluations', 1e5);
% % options = optimset('Display','off','MaxIter', 10000000,'MaxFunEvals', 100000000,'TolX', 1e-12,'TolFun', 1e-12);

try
    [optim_params] = fminunc(@neg_log_posterior, init_vals, options);
catch
    [optim_params] = fminunc(@neg_log_posterior, init_vals, options);
end

    function [nlp] = neg_log_posterior(params)
        nlp = neg_bernoulli_log_likelihood(data, responses, oval_coordinates, num_frames, params);
    end
end

function LL = neg_bernoulli_log_likelihood(data, responses, oval_coordinates, num_frames, params)
scaling = params(1)^2;
dist_param = exp(params(2));
mu = params(3);
sigma = params(4)^2;
temporal = (params(5));
bias = params(6);
% lapse = 1e-4+(1-1e-4)*sigmoid(params(7));%
lapse = params(7)^2;
trials = size(data,1);
num_im = fix(size(oval_coordinates,3)/num_frames);

for fx=1:num_frames
    temporal_weights(((fx-1) * num_im + 1):fx * num_im) = ones(1,num_im) * (exp(temporal*(fx)));
end
for tr=1:trials
    elliptical_pixel_dist = squeeze(sqrt(oval_coordinates(tr,1,:).^2 + dist_param * oval_coordinates(tr,2,:).^2));
    spatial_weights(tr,:) = (exp((-1*(elliptical_pixel_dist-mu).^2)./(sigma^2))) * scaling;
    logits(tr) = dot(data(tr,:),spatial_weights(tr,:) .* temporal_weights) + bias;
end
log_bernoulli = -1*(log(0.5*lapse+(1-lapse)*sigmoid(logits(:))).*responses)-1*(log(1-0.5*lapse-(1-lapse)*sigmoid(logits(:))).*(1-responses));
LL = sum(log_bernoulli);
end
