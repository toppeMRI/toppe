function g = makeevenlength(g)
% If length not even, pad with zeroes at end
% function g = makeevenlength(g)

if mod(length(g),2)
	[nrows,ncols] = size(g);
	if nrows > ncols
		g = [g; zeros(2-mod(length(g),2),ncols)];
	else
		g = [g  zeros(nrows,2-mod(length(g),2))];
	end
end

% EOF
