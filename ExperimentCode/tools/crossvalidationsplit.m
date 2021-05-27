function [training_data, held_out_data] = crossvalidationsplit(data, nfold)
%CROSSVALIDATIONSPLIT split data into 'training' and 'held out' for
%cross-validation.
%
%[training_data, held_out_data] = CROSSVALIDATIONSPLIT(data, nfold) where
%data is an n-dimensional array (size [d1 x d2 x ... x dn]), returns two
%[nfold x 1] cell arrays. The disjoint sets 'training_data{k}' and
%'held_out_data{k}' together reconstruct the full dataset. The size of
%'held_out_data{k}' is approximately d1/nfold

sz = size(data, 1);
assert(nfold <= sz, sprintf('Cannot split data of size %d into %d pieces!', sz, nfold));

% Compute sizes of each held-out component
subset_sizes = diff(round(linspace(0, sz, nfold+1)))';

% Shuffle data along 1st dimension
shuffled_data = reshape(data(randperm(sz), :), size(data));

% Split 'data' into 'nfold' components of approximately equal size. These
% will be held-out data.
sz_cell = num2cell(sz);
sz_cell{1} = subset_sizes;
held_out_data = mat2cell(shuffled_data, sz_cell{:});

% For each subset k, all data 1:k-1 and k+1:end constitute the training
% data.
training_data = cell(size(held_out_data));
for k=1:nfold
    training_data{k} = vertcat(held_out_data{[1:k-1 k+1:nfold]});
end

end

