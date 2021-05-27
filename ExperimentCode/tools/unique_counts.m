function [unq, counts] = unique_counts(A)
%UNIQUE_COUNTS get unique values and counts of each instance in an array.
%
%Where u = unique(A) only gives unique values, [u, c] = UNIQUE_COUNTS(A)
%also includes a vector c the same size as u, where c(i) is the count of
%instances of u(i) in A.
%
% see https://www.mathworks.com/matlabcentral/newsreader/view_thread/309839

sortA = sort(A(:));
% If sortA is [5 5 5 6 7 12 12 12 12] then diff(sortA) will be [0 0 1 1 5 0
% 0 0], and length(A) is 9. So find([9 0 0 1 1 5 0 0 0 9]) gives p=[1 4 5 6
% 10]. Note that sortA(p(1:end-1)) gives unique values and diff(p) gives
% the counts of them!
p = find([length(A); diff(sortA); length(A)]);
unq = sortA(p(1:end-1));
counts = diff(p);
end