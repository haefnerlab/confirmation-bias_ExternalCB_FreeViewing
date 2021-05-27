function [img, c, aperture] = grating_no_annulus(rotation, contrast, frames)
% rotation: a vector with length=frames, each value is clockwise in radian
% phase: number between 0 and 2*pi. passed in as a vector of length=frames
% contrast: w of the dominant direction. 0.5 < w <= 1
% frames: number of frames; set to 10
sf = 18; % spatial freq in cycles per image
width = 121;%481;
annulusPix = 0;%200; %25
phase = 0;
rotation = -rotation;

% Create separate [-1, 1] range meshgrid for pixel-space filters.
[px, py] = meshgrid(linspace(-1, 1, width));
pr = sqrt(px.^2 + py.^2);

% cut out annulus
aperture = exp(-4 * pr.^2);
aperture = aperture .* (1 + erf(10 * (pr - annulusPix / width)));

% blur = 20;
% [xx, yy] = meshgrid(linspace(-2,2,blur));
% kernel = normpdf(sqrt(xx.^2 + yy.^2));

threshold = 0.6;
if contrast<=threshold
    flip = 2;
    delta = 0.01;
else
    flip = 0;
    delta = 0.01;

end
img = zeros(frames, width, width);
cat = ones(1,frames-flip);
cat = [cat -1*ones(1,flip)]; 
cat = Shuffle(cat);
c = zeros(1,frames);
for f=1:frames
    % Z is one frame
    Z1 = (1+sign(sin(sf*(px+py*rotation(f))+phase(f))))/2; 
    Z2 = (1+sign(sin(sf*(px-py*rotation(f))+phase(f))))/2; 
%     blurred_noise = conv2(white_noise, kernel, 'same'); 
%     Z = Z + noise*blurred_noise;
%     Z = reshape(Z,[width,width]);
    if frames==1
        c(f)=contrast;
    else
        c(f) = contrast + cat(f)*delta;
    end
    img(f, :, :) = squeeze(aperture .* max(c(f)*Z1,(1-c(f))*Z2));
    %img(f, :, :) = squeeze(aperture .* (contrast*Z1+(1-contrast)*Z2));
end

if frames>1
    disp('c:');
    disp(c);
end

% [img, ~] = grating([-1,], 0, 1, 1);
% imagesc(squeeze(img));colormap(gray); axis image; axis('off');


end