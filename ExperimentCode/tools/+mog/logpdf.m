function L = logpdf(x, mog)
%MOG.LOGPDF compute log pdf of 1d mixture-of-gaussians more stably than
%log(mog.pdf())
%
%See mog.pdf for more information.
modes = size(mog, 1);
log_mode_prob = zeros(numel(x), modes);
for m=1:modes
    mu = mog(m, 1);
    sigma2 = mog(m, 2)^2;
    logpi = log(mog(m, 3));
    for i=1:numel(x)
        log_mode_prob(i, m) = logpi - ((x(i)-mu)^2 / sigma2 + log(2*pi*sigma2)) / 2;
    end
end
L = logsumexp(log_mode_prob, 2);
L = reshape(L, size(x));
end

function s = logsumexp(a, dim)
% Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
% Default is dim = 1 (columns).
% logsumexp(a, 2) will sum across rows instead of columns.
% Unlike matlab's "sum", it will not switch the summing direction
% if you provide a row vector.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

if nargin < 2
  dim = 1;
end

% subtract the largest in each column
[y, ~] = max(a,[],dim);
dims = ones(1,ndims(a));
dims(dim) = size(a,dim);
a = a - repmat(y, dims);
s = y + log(sum(exp(a),dim));
i = find(~isfinite(y));
if ~isempty(i)
  s(i) = y(i);
end
end