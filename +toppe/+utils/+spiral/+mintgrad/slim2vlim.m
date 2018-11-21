%
%	Function [Imax,Vmax,Rcoil,Lcoil,eta] = slim2vlim(gmax,smax)
%
%	Function takes the constant-slew-rate-limited model and 
%	returns parameters so that voltage limited functions
%	will return the constant-slew-rate-limit answers.
%
%       INPUT:
%               gmax  = maximum gradient (G/cm)
%               smax  = maximum slew (G/cm/s)
%
%       OUTPUT:
%               Imax  = maximum amplifier current (A)
%               Vmax  = maximum amplifier voltage (V)
%               Rcoil = coil resistance (ohms)
%               Lcoil = coil inductance (H)
%               eta   = coil efficiency (G/cm/A)
%
%
%       B. Hargreaves, April 2003.



function [Imax,Vmax,Rcoil,Lcoil,eta] = slim2vlim(gmax,smax)

if (nargin < 2)
	error('Not enough arguments');
end;

%	The translation here is fairly arbitrary... but once 
%	eta and Lcoil are chosen, we are done.  
%
%	Here we choose eta=1 so that it looks like amps=G/cm,
%	and Lcoil = 1 so that volts look like G/cm/s.
	
eta = ones(size(gmax));		% Arbitrary.  Hopefully this is reasonable.
Lcoil = ones(size(gmax));	% Arbitrary.  Again, hope this is reasonable.

Rcoil = zeros(size(gmax));	% Contant-slew does not take into account R.
Imax = gmax./eta; 		% Max current limit.
Vmax = Lcoil.*smax./eta;	% Max voltage limit.


