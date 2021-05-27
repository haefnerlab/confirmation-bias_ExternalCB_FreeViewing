function x = sample(mog, n)
%MOG.SAMPLE sample from a 1d mixture of gaussians.
%
%See MOG.PDF for format.

if nargin < 2, n = 1; end

% First choose a mode.
for i=n:-1:1
	mode(i) = find(rand <= cumsum(mog(:,3)), 1);
end

% Then sample from each chosen mode.
x = mog(mode, 1) + randn(n, 1) .* mog(mode, 2);
end