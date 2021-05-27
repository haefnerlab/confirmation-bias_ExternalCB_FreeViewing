function [ cmap ] = fadecolors( color1, color2, n )
%fadecolors Create a colormap with n points that fades from color1 to
%color2.

cmap = zeros(n, 3);

cmap(:, 1) = linspace(color1(1), color2(1), n)';
cmap(:, 2) = linspace(color1(2), color2(2), n)';
cmap(:, 3) = linspace(color1(3), color2(3), n)';

end