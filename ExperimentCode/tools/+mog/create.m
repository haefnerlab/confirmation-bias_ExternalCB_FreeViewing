function mog = create(mus, sigmas, pis)
%MOG.CREATE create a mixture-of-gaussians.
%
%A MoG (in 1d only) is specified by a mean, standard deviation, and weight
%at each mode. A distribution with N modes is represented by an Nx3 matrix:
%
%   mog = [mu_1 sigma_1 pi_1; 
%          ...;
%          mu_n sigma_n, pi_n]
mog = [mus(:) sigmas(:), pis(:)];
end