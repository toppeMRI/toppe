function f = fermi2d(sz,radius,width,fltgeom)
%  usage ...  fermi2d(sz,radius,width)
%  sz - matrix size, radius - radius of window, width - transition width
%  sz         matrix size
%  radius     radius of window
%  width      transition width
%  fltgeom    'circ' (default) or 'rect'
%  negative width will give "standard width" = 20*radius/128
%
% From Doug Noll.
% jfn added rectangular filter output

if ~exist('fltgeom', 'var')
	fltgeom = 'circ';
end

%
% [x,y]= meshdom(-64:1:63,  63:-1:-64);
% f= 1 ./ (1 + exp((sqrt(x.^2 + y.^2) - radius) ./ (10*steepness/256)));
if width < 0,
   width = 20*radius/128;
end
cent = sz/2 + 1;
x= (1:sz);
y= (1:sz)';
X= ones(size(y))*x;
Y= y*ones(size(x));

X = X-cent;
Y = Y-cent;

if strcmp(fltgeom, 'circ')
	R = abs( X + Y.*1i );
	f = 1 ./ (1 + exp( (R - radius) ./ width ));
else
	fx = 1 ./ (1 + exp( (abs(X) - radius) ./ width ));
	fy = 1 ./ (1 + exp( (abs(Y) - radius) ./ width ));
	f = fx.*fy;
end

return;
