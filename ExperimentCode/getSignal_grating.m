function sig = getSignal_grating(im, noise, template)
%BPG.GETSIGNAL(bpg_im, oriDEG, oriKappa) compute the amount of energy there
%is in the given image(s) at the given orientation band. Spatial Frequency
%is integrated out.
%
%BPG.GETSIGNAL(bpg_im, oriDEG, oriKappa, spFreqCPP, spFreqStdCPP) also
%includes a spatial frequency filter.
%
% bpg_im must be an image (a matrix) or a [frames x height x width] array
% of images.
%
% Returns a [frames x 1] vector of signal levels.
[frames, width, ~] = size(im);
px = width*width;
% % % % persistent K;
% % % % 
% % % % if ismatrix(im)
% % % %     im = double(reshape(im, [1 size(im)]));
% % % % end
% % % % [frames, width, ~] = size(im);
% % % % px = width*width;
% % % % 
% % % % %TODO - if the blur ever changes, this can't stay simply 'persistent'
% % % % %(since the same Cinv is always used regardless of function inputs)
% % % % if isempty(K)
% % % %     blur = 11;
% % % %     [xx, yy] = meshgrid(linspace(-2,2,blur));
% % % %     kernel = normpdf(sqrt(xx.^2 + yy.^2));
% % % %     
% % % %     K = zeros(px, px);
% % % %     for p=1:px
% % % %         impulse = zeros(width);
% % % %         impulse(p) = 1;
% % % %         response = conv2(impulse, kernel, 'same');
% % % %         K(:, p) = response(:);
% % % %     end
% % % % end
% % % % 
% % % % scaledK = sqrt(noise) * K;

%% Get signal of each frame.
sig = zeros(frames, 1);
for f=1:frames
    I = squeeze(reshape(im(f,:,:), [px, 1]));
%     sig(f) = (scaledK \ double(I(:)))' * (scaledK \ double(template(:)));
    sig(f) = (1.0/noise*noise) * double(I(:))' * double(template(:));
end
end