function [sobl, postVal] = LinearPK_with_lapse_duration(data, responses, standardize)
%LinearPK Regress PK as a linear function.
%
% sob = LinearPK(data, responses) returns sob = [slope offset bias].

if nargin < 3, standardize = 0; end

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

% convert boolean to float type
% assert(islogical(responses));
responses = 1.0 * responses(:);

[~, frames] = size(data);

    function NLL = neg_bernoulli_log_likelihood(sobl)
        weights = sobl(2) + (0:frames-1) * sobl(1);
        logits = data * weights(:) + sobl(3);
        lapse = sobl(4)^2;
        neg_log_bernoulli = -log(exp(responses .* logits(:)) + lapse*(responses - 0.5).*(1-exp(logits(:)))) + log(1 + exp(logits(:)));
        %         neg_log_bernoulli = -log(exp(logits(:)).*(responses*(1-lapse)+0.5*lapse) + ((1-responses)*(1-lapse)+0.5*lapse)) + log(1 + exp(logits(:)));
        
        %         neg_log_bernoulli = -responses .* logits(:) + log(1 + exp(logits(:)));
        NLL = sum(neg_log_bernoulli);
    end

compute_error = nargout > 2;
options=optimset('MaxIter', 1000000, 'MaxFunEvals', 1000000,'Display','off');
lb = [-inf*ones(1,3) -1.0];
ub = [inf*ones(1,3) 1.0];
% Fit weights using 'fminunc', only computing the hessian (which is
% slow) if errors are requested.
% if compute_error
%     [sobl, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_bernoulli_log_likelihood, zeros(4,1),options);
%     [sobl, negPostVal, ~, ~, ~, hessian] = fmincon(@neg_bernoulli_log_likelihood, zeros(1,4),[],[],[],[],lb,ub,[],options);
%
%
% attempt to invert the hessian for standard error estimate - this
% sometimes fails silently, returning NaN.
%     errors = sqrt(diag(abs(inv(-hessian))));
% else
%     [sobl, negPostVal] = fminunc(@neg_bernoulli_log_likelihood, zeros(1,4),options);
[sobl, negPostVal] = fmincon(@neg_bernoulli_log_likelihood, zeros(1,4),[],[],[],[],lb,ub,[],options);
% end
postVal = -negPostVal;
end
