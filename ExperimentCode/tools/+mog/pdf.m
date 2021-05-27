function P = pdf(x, mog, discretize)
%MOG.PDF compute pdf of 1d mixture-of-gaussians.
%
%A MoG (in 1d only) is specified by a mean, standard deviation, and weight
%at each mode. A distribution with N modes is represented by an Nx3 matrix:
%
%   mog = [mu_1 sigma_1 pi_1; 
%          ...;
%          mu_n sigma_n, pi_n]
flat_x = x(:)';
modes = size(mog, 1);
if all(mog(:,2) == 0)
    l_each_mode = 1e4 * mog(:,3) .* (flat_x == mog(:, 1));
else
    l_each_mode = normpdf(flat_x .* ones(modes, 1), mog(:,1), mog(:,2));
end
P = l_each_mode' * mog(:,3);
if nargin > 2 && discretize, P = P / sum(P); end
P = reshape(P, size(x));
end