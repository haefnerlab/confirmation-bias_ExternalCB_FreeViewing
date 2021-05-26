function [xDPS, yDPS] = estimateEyeVelocity(tMS, xDVA, yDVA, h)
%ESTIMATEEYEVELOCITY get x and y velocity of eye with center finite difference estimate. That is,
%speed dx/dt is estimated using (x(t+h)-x(t-h))/2h, where h is the 'half width'.
%
%[xDPS, yDPS] = ESTIMATEEYEVELOCITY(t, x, y, h) uses time points t (in milliseconds), xy positions
%given by xDVA and yDVA, and half-width h. Returns x and y speeds in DPS (degrees per second)
%
%xDVA and yDVA may be size [time points x trials]

eyeTimeS = tMS(:) / 1000;

dH = xDVA(1+h:end, :) - xDVA(1:end-h, :);
dV = yDVA(1+h:end, :) - yDVA(1:end-h, :);
dt = eyeTimeS(1+h:end, :) - eyeTimeS(1:end-h, :);

xDPS = dH ./ dt;
yDPS = dV ./ dt;

% Pad beginning/end with repeats of first/last value so that output is same size as input
xDPS = vertcat(repmat(xDPS(1, :), ceil(h/2), 1), xDPS, repmat(xDPS(end, :), floor(h/2), 1));
yDPS = vertcat(repmat(yDPS(1, :), ceil(h/2), 1), yDPS, repmat(yDPS(end, :), floor(h/2), 1));
end