function res = selectNth(x, N, direction)
% Select the Nth ordered element from x, largest or smallest by 'direction'.

if ~exist('direction', 'var'), direction = 'descend'; end

[~, ind] = sort(x(:), direction);

res = x(ind(N));

% Maybe I should implement a Median of Medians based selecting alg.
% since it's \theta(n)

end

