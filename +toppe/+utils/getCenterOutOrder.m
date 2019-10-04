function vco = getCenterOutOrder(v)

n = length(v);

if mod(n,2)
	error('length(v) must be even');
end

for i = 1:n/2
	vco(2*(i-1)+1) = n/2 + i;
	vco(2*(i-1)+2) = n/2 - i + 1;
end

% Or:
%I = 1:n/2;
%vco(2*(I-1)+1) = n/2 + I;
%vco(2*(I-1)+2) = n/2 - I + 1;


