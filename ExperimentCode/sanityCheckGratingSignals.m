T1 = grating(+pi/4, 0, 20, 1);
T2 = grating(-pi/4, 0, 20, 1);

contrast = 20;
noise_levels = [0.01, 0.1, 1, 10];
repeats = 1000;

sigsPlus = zeros(length(noise_levels), repeats);
sigsMinus = zeros(length(noise_levels), repeats);

for idx=1:length(noise_levels)
    n = noise_levels(idx);
    g = grating(+pi / 4 * ones(1, repeats), n, contrast, repeats);
    sigsPlus(idx, :) = getSignal_grating(g, n, T1) - getSignal_grating(g, n, T2);
    g = grating(-pi / 4 * ones(1, repeats), n, contrast, repeats);
    sigsMinus(idx, :) = getSignal_grating(g, n, T1) - getSignal_grating(g, n, T2);
end

%% Plot

figure;
hold on;
for idx=1:length(noise_levels)
    allsigs = horzcat(sigsPlus(idx, :), sigsMinus(idx, :));
    ksdensity(zscore(allsigs));
end
legend(arrayfun(@(n) sprintf('noise=%.1e', n), noise_levels, 'UniformOutput', false));

%% Same as above but for contrast

noise = 10;
contrasts = linspace(0, 64, 5);
repeats = 1000;

sigsPlus = zeros(length(contrasts), repeats);
sigsMinus = zeros(length(contrasts), repeats);

for idx=1:length(contrasts)
    c = contrasts(idx);
    T1 = grating(+pi/4, 0, c, 1);
    T2 = grating(-pi/4, 0, c, 1);
    g = grating(+pi / 4 * ones(1, repeats), noise, c, repeats);
    sigsPlus(idx, :) = getSignal_grating(g, noise, T1) - getSignal_grating(g, noise, T2);
    g = grating(-pi / 4 * ones(1, repeats), noise, c, repeats);
    sigsMinus(idx, :) = getSignal_grating(g, noise, T1) - getSignal_grating(g, noise, T2);
end

%% Plot

figure;
hold on;
for idx=1:length(contrasts)
    allsigs = horzcat(sigsPlus(idx, :), sigsMinus(idx, :));
    ksdensity(zscore(allsigs));
end
legend(arrayfun(@(c) sprintf('contrast=%.1f', c), contrasts, 'UniformOutput', false));