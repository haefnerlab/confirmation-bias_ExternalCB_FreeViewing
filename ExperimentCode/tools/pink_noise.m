function noise = pink_noise(size)

I = randn(size);

IF = fftshift(fft2(I));
xs = linspace(-1, 1, size);

[xx, yy] = meshgrid(xs, xs);
rr = sqrt(xx.^2 + yy.^2);

windowF = exp(-rr * 5);

noise = 3 * real(ifft2(ifftshift(windowF .* IF)));

end