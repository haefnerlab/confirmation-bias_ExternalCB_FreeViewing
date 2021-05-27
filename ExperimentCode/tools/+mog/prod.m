function prod = prod(d1, d2)
%MOG.PROD compute product of two 1d mixture-of-gaussian distributions.
%
%See MOG.PDF for format.
modes1 = size(d1, 1);
modes2 = size(d2, 1);
modes_out = modes1 * modes2;

prod = zeros(modes_out, 3);
k = 1;
for i=1:modes1
    for j=1:modes2
        mu_i = d1(i,1); var_i = d1(i,2)^2; pi_i = d1(i,3);
        mu_j = d2(j,1); var_j = d2(j,2)^2; pi_j = d2(j,3);
        mu_k = (var_i * mu_j + var_j * mu_i) / (var_i + var_j);
        var_k = (var_i * var_j) / (var_i + var_j);
        pi_k = pi_i * pi_j * normpdf(mu_i, mu_j, sqrt(var_i + var_j));
        prod(k, :) = [mu_k, sqrt(var_k), pi_k];
        k = k+1;
    end
end
prod(:,3) = prod(:,3) / sum(prod(:,3));
end