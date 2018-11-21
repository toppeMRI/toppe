
%function [N] = numtrailzeros(g,geps)
%
%	Function looks at the gradient and returns
%	the number of points at the end with magnitude
%	less than geps.
%
%	INPUT:
%		g = (Ng x D) gradient samples, D=#dimensions.
%		geps = (1x1) tolerance.
%
%	OUTPUT:
%		N = number of points at end that are within geps of zero.

function [N] = numtrailzeros(g,geps)

sg = size(g);
Ng = sg(1);
g = g.';
normg = sqrt(sum(g.*g,1));

f = find(normg > geps);
N = Ng-max(f);


